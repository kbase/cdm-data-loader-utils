"""Shared model-registry and scenario-directory fixtures for the JSONL pipeline tests."""

import sys
import types
from collections.abc import Callable, Generator
from pathlib import Path
from uuid import uuid4

import pytest
from pydantic import BaseModel, Field


class Widget(BaseModel):
    """Test model. widget_id must be non-empty. count must be zero or more."""

    widget_id: str = Field(min_length=1)
    count: int = Field(ge=0)


@pytest.fixture
def widget_model() -> type[Widget]:
    """Return the Widget model."""
    return Widget


@pytest.fixture
def entity_models_module_factory() -> Generator[Callable[[dict[str, type[BaseModel]]], str]]:
    """Return a factory that registers a module holding an ENTITY_MODELS dict.

    Each call creates a module with a unique name and adds it to
    sys.modules. Removes each created module after the test.
    """
    created_names: list[str] = []

    def _factory(entity_models: dict[str, type[BaseModel]]) -> str:
        module_name = f"_test_entity_models_{uuid4().hex}"
        module = types.ModuleType(module_name)
        module.ENTITY_MODELS = entity_models  # type: ignore[attr-defined]
        sys.modules[module_name] = module
        created_names.append(module_name)
        return module_name

    yield _factory

    for name in created_names:
        sys.modules.pop(name, None)


@pytest.fixture
def entity_models_module(
    entity_models_module_factory: Callable[[dict[str, type[BaseModel]]], str], widget_model: type[Widget]
) -> str:
    """Return a module name that maps 'widget' to the Widget model."""
    return entity_models_module_factory({"widget": widget_model})


@pytest.fixture
def scenario_input_dir(test_data_dir: Path) -> Callable[[str], str]:
    """Return a function that maps a scenario name to its data directory."""

    def _resolve(scenario: str) -> str:
        path = test_data_dir / "jsonlines" / scenario
        assert path.is_dir(), f"missing test fixture directory: {path}"
        return str(path)

    return _resolve
