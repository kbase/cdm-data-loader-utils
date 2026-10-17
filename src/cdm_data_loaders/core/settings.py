"""Common defaults for running CDM data loading pipelines."""

import logging
from typing import TYPE_CHECKING, Any, Final, Self

from frozendict import frozendict
from pydantic import AliasChoices, AliasPath, computed_field, field_validator, model_validator
from pydantic.fields import FieldInfo
from pydantic_settings import BaseSettings, SettingsConfigDict, SettingsError

from cdm_data_loaders.core.destination import PIPELINE_METADATA_NOT_LOCAL, is_local, normalise_dir
from cdm_data_loaders.core.fields import (
    BUFFER_SIZE,
    DEFAULTS,
    DLT_DEV_MODE,
    INPUT_DIR,
    LOADER_FILE_FORMAT,
    LOG_CONFIG_FILE,
    LOG_INTERVAL,
    OUTPUT_DIR,
    USE_DESTINATION,
    USE_OUTPUT_DIR_FOR_PIPELINE_METADATA,
    BufferSize,
    DatasetName,
    DltDevMode,
    InputDir,
    LoaderFileFormat,
    LogConfigFile,
    LogInterval,
    OutputDir,
    UseDestination,
    UseOutputDirForPipelineMetadata,
)

if TYPE_CHECKING:
    from collections.abc import Mapping

logger = logging.getLogger(__name__)


DEFAULT_CTS_SETTINGS = frozendict(
    {
        k: DEFAULTS[k]
        for k in [
            BUFFER_SIZE,
            DLT_DEV_MODE,
            INPUT_DIR,
            LOG_CONFIG_FILE,
            LOG_INTERVAL,
            OUTPUT_DIR,
            USE_DESTINATION,
            USE_OUTPUT_DIR_FOR_PIPELINE_METADATA,
        ]
    }
)


DEFAULT_SETTINGS_CONFIG_DICT = frozendict(
    {
        "cli_exit_on_error": False,
        "cli_ignore_unknown_args": True,
        "cli_kebab_case": True,
        "env_prefix": "CDL_",
        "str_strip_whitespace": True,
    }
)

CLI_SHORTCUTS: frozendict[str, str] = frozendict(
    {
        INPUT_DIR: "i",
        OUTPUT_DIR: "o",
        USE_DESTINATION: "d",
        USE_OUTPUT_DIR_FOR_PIPELINE_METADATA: "p",
    }
)

# argparse registers -h/--help on every parser; claiming either makes the parser fail to build
RESERVED_CLI_NAMES: Final[frozenset[str]] = frozenset({"h", "help"})


def _alias_names(alias: str | AliasChoices | AliasPath | None) -> set[str]:
    """Top-level names introduced by a field alias or validation alias."""
    if alias is None:
        return set()
    if isinstance(alias, str):
        return {alias}
    if isinstance(alias, AliasPath):
        return {alias.path[0]} if isinstance(alias.path[0], str) else set()
    return set().union(*(_alias_names(choice) for choice in alias.choices))


def default_settings_with_shortcuts(
    *,
    cli_prog_name: str | None = None,
    cli_shortcuts: dict[str, str | list[str]] | frozendict[str, str | list[str]] | None = None,
) -> SettingsConfigDict:
    """Build a SettingsConfigDict with the default settings config and the given CLI shortcuts."""
    extra_args = {}
    if cli_prog_name is not None:
        extra_args["cli_prog_name"] = cli_prog_name
    if cli_shortcuts is not None:
        extra_args["cli_shortcuts"] = {
            field_name.replace("_", "-"): value for field_name, value in cli_shortcuts.items()
        }
    return SettingsConfigDict(**{**DEFAULT_SETTINGS_CONFIG_DICT, **extra_args})


class CdmDataLoadersBase(BaseSettings):
    """Base for all CDM Data Loaders settings classes.

    - Sets up some basic CLI parsing defaults
    - Adds a class method to check for conflicting aliases
    """

    model_config = SettingsConfigDict(**DEFAULT_SETTINGS_CONFIG_DICT)

    @classmethod
    def __pydantic_init_subclass__(cls, **kwargs: Any) -> None:
        """Validate the subclass's cli shortcuts and aliases."""
        super().__pydantic_init_subclass__(**kwargs)
        cls.check_aliases()

    @classmethod
    def _cli_name(cls, name: str) -> str:
        """Get the CLI argument name pydantic-settings registers for a field or alias name."""
        return name.replace("_", "-") if cls.model_config.get("cli_kebab_case") else name

    @classmethod
    def _field_cli_names(cls, field_name: str, field_info: FieldInfo) -> set[str]:
        """CLI argument names a field may register: its name plus any aliases.

        Deliberately over-inclusive -- reserving a name that is never registered can only cause
        a false positive, never a missed collision.
        """
        names = {field_name} | _alias_names(field_info.alias) | _alias_names(field_info.validation_alias)
        return {cls._cli_name(name) for name in names}

    @classmethod
    def check_aliases(cls) -> None:
        """Raise `SettingsError` if the class's CLI argument names or shortcuts are invalid or collide.

        pydantic-settings silently drops shortcuts whose target does not match a registered argument
        name, and silently drops any name already registered by an earlier argument, so both are
        checked here. Names differing only in '_' vs '-' are treated as the same name.
        """
        owners: dict[str, str] = {}
        errors = []

        def claim(name: str, owner: str, errors: list[str]) -> None:
            key = name.replace("_", "-")
            if key in RESERVED_CLI_NAMES:
                err_msg = f"{cls.__name__}: CLI name {name!r} for {owner!r} is reserved by argparse"
                errors.append(err_msg)
                return
            other = owners.setdefault(key, owner)
            if other == owner and key in claimed:
                err_msg = f"{cls.__name__}: CLI name {name!r} is claimed twice by {owner!r}"
                errors.append(err_msg)
                return
            if other != owner:
                both = " and ".join(sorted(repr(o) for o in (other, owner)))
                err_msg = f"{cls.__name__}: CLI name {name!r} is claimed by both {both}"
                errors.append(err_msg)
                return
            claimed.add(key)

        claimed: set[str] = set()
        registered: dict[str, str] = {}
        for field_name, field_info in cls.model_fields.items():
            cli_names = cls._field_cli_names(field_name, field_info)
            # claim each normalised name once per field, so a field's own name and alias can coincide
            for cli_name in {name.replace("_", "-"): name for name in cli_names}.values():
                claim(cli_name, field_name, errors)
            registered.update(dict.fromkeys(cli_names, field_name))

        cli_shortcuts: Mapping[str, str | list[str]] = cls.model_config.get("cli_shortcuts") or {}
        for target, shortcuts in cli_shortcuts.items():
            owner = registered.get(target)
            if owner is None:
                suggestion = cls._cli_name(target)
                hint = f" (did you mean {suggestion!r}?)" if suggestion in registered else ""
                err_msg = f"{cls.__name__}: cli_shortcuts target {target!r} does not match any CLI argument{hint}"
                errors.append(err_msg)
                continue
            for shortcut in [shortcuts] if isinstance(shortcuts, str) else shortcuts:
                if not shortcut or shortcut.startswith("-"):
                    err_msg = f"{cls.__name__}: invalid shortcut {shortcut!r} for {owner!r}; omit the leading dashes"
                    errors.append(err_msg)
                    continue
                claim(shortcut, owner, errors)
        if errors:
            raise SettingsError("\n".join(errors))


class LoggerSettings(CdmDataLoadersBase):
    """Configuration for a class with a logger config."""

    log_config_file: LogConfigFile


class InputOutputSettings(LoggerSettings):
    """Configuration with basic input and output settings.

    ``output_dir`` is None when not given. A blank value (e.g. ``CDL_OUTPUT_DIR=""``) also counts as not given.
    """

    model_config = default_settings_with_shortcuts(cli_shortcuts={INPUT_DIR: "i", OUTPUT_DIR: "o"})

    input_dir: InputDir
    output_dir: OutputDir

    @field_validator(OUTPUT_DIR, mode="before")
    @classmethod
    def blank_output_dir_is_unset(cls, value: Any) -> Any:  # noqa: ANN401
        """Treat a blank output_dir as unset, so it can fall back to the destination config."""
        if isinstance(value, str) and not value.strip():
            return None
        return value

    @field_validator(INPUT_DIR, OUTPUT_DIR, mode="after")
    @classmethod
    def validate_dir_path(cls, value: str | None) -> str | None:
        """Remove any trailing slashes from directory paths."""
        return None if value is None else normalise_dir(value)


class CtsSettings(InputOutputSettings):
    """Configuration for running a basic DLT pipeline.

    Building this class only parses the CLI and env vars; it never reads dlt.config. If
    ``output_dir`` is not given, it stays None until
    :func:`cdm_data_loaders.pipelines.core.resolve_cts_settings` fills it in from
    ``destination.<use_destination>.bucket_url``. ``run_cli`` does this automatically.
    """

    model_config: SettingsConfigDict = default_settings_with_shortcuts(cli_shortcuts=CLI_SHORTCUTS)

    buffer_size: BufferSize
    dlt_dev_mode: DltDevMode
    log_interval: LogInterval
    use_destination: UseDestination
    use_output_dir_for_pipeline_metadata: UseOutputDirForPipelineMetadata

    @model_validator(mode="after")
    def check_pipeline_metadata_location(self) -> Self:
        """Fail early if an explicit output_dir is remote and pipeline metadata is meant to go there.

        Output locations taken from the dlt config are checked during resolution instead.
        """
        if self.use_output_dir_for_pipeline_metadata and self.output_dir and not is_local(self.output_dir):
            raise ValueError(PIPELINE_METADATA_NOT_LOCAL.format(url=self.output_dir))
        return self

    def _output_subdir(self, name: str) -> str | None:
        if self.output_dir is None:
            return None
        separator = "" if self.output_dir.endswith("/") else "/"
        return f"{self.output_dir}{separator}{name}"

    @computed_field
    @property
    def output_is_local(self) -> bool | None:
        """Whether output_dir is on the local filesystem; None until output_dir is resolved."""
        return None if self.output_dir is None else is_local(self.output_dir)

    @computed_field
    @property
    def raw_data_dir(self) -> str | None:
        """Directory to save downloaded raw data files to: ``<output_dir>/raw_data``.

        None until output_dir is resolved.
        """
        return self._output_subdir("raw_data")

    @computed_field
    @property
    def pipeline_dir(self) -> str | None:
        """Custom directory to save pipeline metadata to.

        If use_output_dir_for_pipeline_metadata is true, this is ``<output_dir>/.dlt_conf``;
        otherwise None, and dlt uses its default.
        """
        if self.use_output_dir_for_pipeline_metadata:
            return self._output_subdir(".dlt_conf")
        return None


class FileFormatIngestSettings(CtsSettings):
    """Configuration for DLT pipelines that import a particular file format."""

    model_config: SettingsConfigDict = default_settings_with_shortcuts(
        cli_shortcuts={**CLI_SHORTCUTS, LOADER_FILE_FORMAT: "f"}
    )

    dataset_name: DatasetName
    loader_file_format: LoaderFileFormat
