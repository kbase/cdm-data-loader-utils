"""Unit tests for the module-level settings accessors and the entity model registry loader."""

import sys
import types
from collections.abc import Callable
from typing import Any

import pytest
from pydantic import BaseModel

from cdm_data_loaders.pipelines.jsonlines.base import load_registry_dict
from cdm_data_loaders.pipelines.jsonlines.extract_pydantic_validate_pipeline import load_entity_models


def test_load_entity_models_pass_returns_registered_mapping(
    entity_models_module_factory: Callable[[dict[str, type[BaseModel]]], str],
    widget_model: type[BaseModel],
) -> None:
    """load_entity_models returns exactly the ENTITY_MODELS dict a module defines."""
    module_name = entity_models_module_factory({"widget": widget_model})
    assert load_entity_models(module_name) == {"widget": widget_model}


def test_load_entity_models_fail_module_not_found() -> None:
    """load_entity_models propagates ImportError for a nonexistent dotted path."""
    with pytest.raises(ImportError):
        load_entity_models("this.module.does.not.exist")


@pytest.mark.parametrize(
    ("attribute_present", "value", "match"),
    [
        pytest.param(False, None, "ENTITY_MODELS", id="attribute_missing_entirely"),
        pytest.param(True, {}, "ENTITY_MODELS", id="empty_dict"),
        pytest.param(True, ["widget"], "ENTITY_MODELS", id="wrong_type_list_instead_of_dict"),
    ],
)
def test_load_entity_models_fail_invalid_registry(
    attribute_present: bool,
    value: object,
    match: str,
) -> None:
    """load_entity_models raises ValueError for a missing, empty, or wrongly-typed ENTITY_MODELS attribute."""
    module_name = f"_test_invalid_registry_{id(value)}"
    module = types.ModuleType(module_name)
    if attribute_present:
        module.ENTITY_MODELS = value  # type: ignore[attr-defined]
    sys.modules[module_name] = module
    try:
        with pytest.raises(ValueError, match=match):
            load_entity_models(module_name)
    finally:
        sys.modules.pop(module_name, None)


def test_load_registry_dict_pass_returns_mapping(
    schema_module_factory: Callable[[dict[str, dict[str, Any]]], str],
    widget_schema: dict[str, Any],
) -> None:
    """load_registry_dict returns a non-empty dict attribute from the imported module."""
    expected = {"widget": widget_schema}
    module_name = schema_module_factory(expected)
    assert load_registry_dict(module_name, "ENTITY_SCHEMAS") == expected


@pytest.mark.parametrize(
    ("attribute_present", "value", "type_hint", "expected_message"),
    [
        pytest.param(
            False,
            None,
            "",
            "_test_load_registry_dict_invalid does not define a non-empty CUSTOM_REG",
            id="missing_no_type_hint",
        ),
        pytest.param(
            True,
            {},
            " dict[str, int]",
            "_test_load_registry_dict_invalid does not define a non-empty CUSTOM_REG dict[str, int]",
            id="empty_with_type_hint",
        ),
        pytest.param(
            True,
            ("not", "a", "dict"),
            " dict[str, Any]",
            "_test_load_registry_dict_invalid does not define a non-empty CUSTOM_REG dict[str, Any]",
            id="wrong_type_with_type_hint",
        ),
    ],
)
def test_load_registry_dict_fail_invalid_attribute(
    monkeypatch: pytest.MonkeyPatch,
    attribute_present: bool,
    value: object,
    type_hint: str,
    expected_message: str,
) -> None:
    """load_registry_dict raises ValueError with optional type_hint when attribute is missing, empty, or not a dict."""
    module_name = "_test_load_registry_dict_invalid"
    module = types.ModuleType(module_name)
    if attribute_present:
        module.CUSTOM_REG = value  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, module_name, module)
    with pytest.raises(ValueError, match=r"^" + expected_message.replace("[", r"\[").replace("]", r"\]") + r"$"):
        load_registry_dict(module_name, "CUSTOM_REG", type_hint=type_hint)
