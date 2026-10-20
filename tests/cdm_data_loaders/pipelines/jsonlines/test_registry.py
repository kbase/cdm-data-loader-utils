"""Unit tests for the module-level settings accessors and the entity model registry loader."""

import sys
import types
from collections.abc import Callable

import pytest
from pydantic import BaseModel

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
