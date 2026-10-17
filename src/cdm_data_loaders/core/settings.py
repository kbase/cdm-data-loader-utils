"""Settings classes for data loading and wrangling pipelines."""

from collections import defaultdict
from pathlib import Path
from typing import TYPE_CHECKING, Any, Final

from frozendict import frozendict
from pydantic import AliasChoices, AliasPath, computed_field
from pydantic.fields import FieldInfo
from pydantic_settings import BaseSettings, SettingsConfigDict, SettingsError

from cdm_data_loaders.core.fields import (
    INPUT_DIR,
    LOADER_FILE_FORMAT,
    OUTPUT_DIR,
    USE_DESTINATION,
    USE_OUTPUT_DIR_FOR_PIPELINE_METADATA,
    BufferSize,
    DatasetName,
    DltDevMode,
    LoaderFileFormat,
    LogConfigFile,
    LogInterval,
    UseDestination,
    UseOutputDirForPipelineMetadata,
)
from cdm_data_loaders.core.path_fields import LocalInputDir, LocalOutputDir

if TYPE_CHECKING:
    from collections.abc import Mapping


DEFAULT_SETTINGS_CONFIG_DICT = frozendict(
    {
        "cli_exit_on_error": False,
        "cli_ignore_unknown_args": False,
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
# dir where dlt will put run metadata, schemas, etc.
PIPELINE_DIR_NAME: Final[str] = ".dlt_conf"


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
        owners = defaultdict(set)
        for field_name, value in cli_shortcuts.items():
            value_list = [value] if not isinstance(value, list) else value
            for v in value_list:
                owners[v].add(field_name.replace("_", "-"))
        violations = {w2: w1s for w2, w1s in owners.items() if len(w1s) > 1}
        if violations:
            err_msg = "Invalid shortcut/name pairs found"
            raise ValueError(err_msg)

        extra_args["cli_shortcuts"] = {
            field_name.replace("_", "-"): value for field_name, value in cli_shortcuts.items()
        }
    return SettingsConfigDict(**{**DEFAULT_SETTINGS_CONFIG_DICT, **extra_args})


def _alias_names(alias: str | AliasChoices | AliasPath | None) -> set[str]:
    """Top-level names introduced by a field alias or validation alias."""
    if alias is None:
        return set()
    if isinstance(alias, str):
        return {alias}
    if isinstance(alias, AliasPath):
        return {alias.path[0]} if isinstance(alias.path[0], str) else set()
    return set().union(*(_alias_names(choice) for choice in alias.choices))


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

    # TODO: simplify this bs
    @classmethod
    def check_aliases(cls) -> None:
        """Raise `SettingsError` if the class's CLI argument names or shortcuts are invalid or collide.

        pydantic-settings silently drops shortcuts whose target does not match a registered argument
        name, and silently drops any name already registered by an earlier argument, so both are
        checked here. Names differing only in '_' vs '-' are treated as the same name.
        """
        # owners = defaultdict(set)
        # for w1, w2 in set(pairs):  # set() drops exact duplicates
        #     owners[w2].add(w1)

        # violations = {w2: w1s for w2, w1s in owners.items() if len(w1s) > 1}

        owners: dict[str, str] = {}
        errors = []
        claimed: set[str] = set()
        registered: dict[str, str] = {}
        dd = defaultdict(set)
        pairs: set[tuple[str, str]] = set()

        def claim(name: str, owner: str, errors: list[str]) -> None:
            key = name.replace("_", "-")
            dd[owner].add(name)
            #     if key in RESERVED_CLI_NAMES:
            #         err_msg = f"{cls.__name__}: CLI name {name!r} for {owner!r} is reserved by argparse"
            #         errors.append(err_msg)
            #         return
            other = owners.setdefault(key, owner)

        #     if other == owner and key in claimed:
        #         err_msg = f"{cls.__name__}: CLI name {name!r} is claimed twice by {owner!r}"
        #         errors.append(err_msg)
        #         return
        #     if other != owner:
        #         both = " and ".join(sorted(repr(o) for o in (other, owner)))
        #         err_msg = f"{cls.__name__}: CLI name {name!r} is claimed by both {both}"
        #         errors.append(err_msg)
        #         return
        #     claimed.add(key)

        for field_name, field_info in cls.model_fields.items():
            cli_names = cls._field_cli_names(field_name, field_info)
            # claim each normalised name once per field, so a field's own name and alias can coincide
            for cli_name in {name.replace("_", "-"): name for name in cli_names}.values():
                pairs.add((cli_name, field_name))
                claim(cli_name, field_name, errors)
            registered.update(dict.fromkeys(cli_names, field_name))

        cli_shortcuts: Mapping[str, str | list[str]] = cls.model_config.get("cli_shortcuts") or {}
        for target, shortcuts in cli_shortcuts.items():
            owner = registered.get(target)
            # if owner is None:
            #     suggestion = cls._cli_name(target)
            #     hint = f" (did you mean {suggestion!r}?)" if suggestion in registered else ""
            #     err_msg = f"{cls.__name__}: cli_shortcuts target {target!r} does not match any CLI argument{hint}"
            #     errors.append(err_msg)
            #     continue
            for shortcut in [shortcuts] if isinstance(shortcuts, str) else shortcuts:
                pairs.add((shortcut, target))
                if not shortcut or shortcut.startswith("-"):
                    err_msg = f"{cls.__name__}: invalid shortcut {shortcut!r} for {owner!r}; omit the leading dashes"
                    errors.append(err_msg)
                    continue
                claim(shortcut, (owner or target), errors)
        violations = {w2: w1s for w2, w1s in dd.items() if len(w1s) > 1}
        if violations:
            print(violations)
        if errors:
            raise SettingsError("\n".join(errors))


class LoggerSettings(CdmDataLoadersBase):
    """Configuration for a class with a logger config."""

    log_config_file: LogConfigFile


class LocalInputOutputSettings(LoggerSettings):
    """Configuration for a DLT pipeline where both input and output directories are local."""

    model_config = default_settings_with_shortcuts(cli_shortcuts={INPUT_DIR: "i", OUTPUT_DIR: "o"})

    input_dir: LocalInputDir
    output_dir: LocalOutputDir


class CtsSettings(LocalInputOutputSettings):
    """Configuration for running a basic DLT pipeline."""

    model_config: SettingsConfigDict = default_settings_with_shortcuts(cli_shortcuts=CLI_SHORTCUTS)

    buffer_size: BufferSize
    dlt_dev_mode: DltDevMode
    log_interval: LogInterval
    use_destination: UseDestination
    use_output_dir_for_pipeline_metadata: UseOutputDirForPipelineMetadata

    @computed_field
    @property
    def pipeline_dir(self) -> Path | None:
        """Custom directory to save pipeline metadata to.

        If use_output_dir_for_pipeline_metadata is true, this is ``<output_dir>/.dlt_conf``;
        otherwise None, and dlt uses its default.
        """
        if self.use_output_dir_for_pipeline_metadata:
            return self.output_dir / PIPELINE_DIR_NAME
        return None


class FileFormatIngestSettings(CtsSettings):
    """Configuration for DLT pipelines that import a particular file format."""

    model_config: SettingsConfigDict = default_settings_with_shortcuts(
        cli_shortcuts={**CLI_SHORTCUTS, LOADER_FILE_FORMAT: "f"}
    )

    dataset_name: DatasetName
    loader_file_format: LoaderFileFormat
