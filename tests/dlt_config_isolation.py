"""Shared dlt config isolation helpers for tests.

These wrap dlt's own test-isolation technique (see `dlt-hub/dlt`'s `tests/utils.py`,
`inject_providers`/`reset_providers`): swap the config provider chain held by the
`PluggableRunContext` in dlt's global `Container` for disposable in-memory providers,
then restore the original chain on exit. `dlt.config`/`dlt.secrets` remain the real
dlt accessor objects throughout -- only the underlying provider chain changes -- so
code under test (including a real `dlt.pipeline()`/`dlt.destination()`) behaves exactly
as it would in production, without ever touching the process's real config.toml/
secrets.toml/env-derived provider chain.
"""

from collections.abc import Generator
from contextlib import contextmanager
from copy import deepcopy
from typing import Any

import dlt
import pytest
from dlt.common.configuration.container import Container
from dlt.common.configuration.providers import EnvironProvider
from dlt.common.configuration.providers.doc import BaseDocProvider
from dlt.common.configuration.specs.config_providers_context import ConfigProvidersContainer
from dlt.common.configuration.specs.pluggable_run_context import PluggableRunContext
from dlt.common.pipeline import PipelineContext


class InMemoryConfigProvider(BaseDocProvider):
    """Writable, non-secret, in-memory config provider.

    _ConfigAccessor.writable_provider (in dlt/common/configuration/accessors.py) requires
    "the first writable provider that does not support secrets", so this class must not support secrets.
    """

    def __init__(self, initial: dict[str, Any] | None = None) -> None:
        """Instantiation."""
        super().__init__(deepcopy(initial or {}))

    def reseed(self, initial: dict[str, Any] | None = None) -> None:
        """Replace this provider's config doc with a deep copy of initial."""
        self._config_doc = deepcopy(initial) if initial else {}

    @property
    def name(self) -> str:
        """Name."""
        return "InMemoryConfigProvider"

    @property
    def supports_secrets(self) -> bool:
        """Supports secrets."""
        return False

    @property
    def supports_sections(self) -> bool:
        """Supports sections."""
        return True

    @property
    def is_writable(self) -> bool:
        """Is writable."""
        return True


class InMemorySecretsProvider(InMemoryConfigProvider):
    """Writable, secrets-capable, in-memory provider. The in-memory analogue of SecretsTomlProvider."""

    @property
    def name(self) -> str:
        """Name."""
        return "InMemorySecretsProvider"

    @property
    def supports_secrets(self) -> bool:
        """Supports secrets."""
        return True


@contextmanager
def dlt_config_unset() -> Generator[None]:
    """Set dlt.config to None. Note: this is not possible using dlt and only exists for testing purposes."""
    with pytest.MonkeyPatch.context() as m:
        m.setattr(dlt, "config", None)
        yield


@contextmanager
def deactivated_pipeline() -> Generator[None]:
    """Deactivate any dlt pipeline active before or created during the wrapped block.

    Mirrors dlt's own `deactivate_pipeline` test fixture: a `dlt.pipeline()` created inside
    the block must not leak forward as the "active pipeline" once the block exits.
    """
    container = Container()
    if container[PipelineContext].is_active():
        container[PipelineContext].deactivate()
    try:
        yield
    finally:
        if container[PipelineContext].is_active():
            container[PipelineContext].deactivate()


@contextmanager
def isolated_dlt_config(
    config: dict[str, Any] | None = None,
    secrets: dict[str, Any] | None = None,
) -> Generator[InMemoryConfigProvider]:
    """Swap dlt's config provider chain for isolated in-memory providers, restored on exit.

    Keeps `EnvironProvider` in the chain, ahead of the in-memory providers, so environment
    variable overrides (`monkeypatch.setenv`, moto/CEPH AWS credentials, etc.) still resolve
    exactly as they do against dlt's real provider chain.
    """
    run_ctx: PluggableRunContext = Container()[PluggableRunContext]
    old_providers: ConfigProvidersContainer = run_ctx.providers
    config_provider = InMemoryConfigProvider(config)
    secrets_provider = InMemorySecretsProvider(secrets)
    run_ctx.providers = ConfigProvidersContainer(
        initial_providers=[EnvironProvider(), config_provider, secrets_provider]
    )
    try:
        with deactivated_pipeline():
            yield config_provider
    finally:
        run_ctx.providers = old_providers
