"""Unit tests for jsonschema_to_pyspark.dereferencing."""

from typing import Any

import pytest
from frozendict import frozendict
from referencing import Registry, Resource, Specification
from referencing.jsonschema import DRAFT7, DRAFT201909, DRAFT202012

from cdm_data_loaders.converters.jsonschema_to_pyspark.dereferencer import (
    DereferenceContext,
    DereferencingError,
    RefSiblingMode,
    _dereference,
    _merge_all_of,
    _resolve_ref_and_all_of,
    _should_apply_ref_siblings,
    dereference_schema,
)


def _make_ctx(
    schema: dict[str, Any],
    spec: Specification = DRAFT202012,
    ref_sibling_mode: RefSiblingMode = "auto",
    max_ref_depth: int = 50,
    additional_resources: dict[str, dict[str, Any]] | None = None,
) -> DereferenceContext:
    """Build a DereferenceContext wired up with a resolver rooted at `schema`."""
    root_resource = Resource.from_contents(schema, default_specification=spec)
    resources = [
        (uri, Resource.from_contents(doc, default_specification=spec))
        for uri, doc in (additional_resources or {}).items()
    ]
    registry = Registry().with_resources(resources)
    resolver = registry.resolver_with_root(root_resource)
    return DereferenceContext(
        resolver=resolver, spec=spec, max_ref_depth=max_ref_depth, ref_sibling_mode=ref_sibling_mode
    )


"""dereference_schema"""


def test_dereference_schema_pass_resolves_local_ref() -> None:
    """dereference_schema() inlines a local $ref pointing into '$defs'."""
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "properties": {"addr": {"$ref": "#/$defs/Address"}},
        "$defs": {"Address": {"type": "string"}},
    }
    result = dereference_schema(schema)
    assert result["properties"]["addr"] == {"type": "string"}
    assert "$ref" not in str(result)


def test_dereference_schema_pass_resolves_external_ref() -> None:
    """dereference_schema() resolves a $ref into a document supplied via additional_resources."""
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "properties": {"addr": {"$ref": "https://example.com/common.json#/Address"}},
    }
    resources = {"https://example.com/common.json": {"Address": {"type": "string"}}}
    result = dereference_schema(schema, additional_resources=resources)
    assert result["properties"]["addr"] == {"type": "string"}


def test_dereference_schema_pass_merges_all_of() -> None:
    """dereference_schema() merges 'allOf' branches' properties/required into a single schema."""
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "allOf": [
            {"properties": {"a": {"type": "string"}}, "required": ["a"]},
            {"properties": {"b": {"type": "integer"}}, "required": ["b"]},
        ],
    }
    result = dereference_schema(schema)
    assert "allOf" not in result
    assert set(result["properties"]) == {"a", "b"}
    assert set(result["required"]) == {"a", "b"}


def test_dereference_schema_pass_removes_definition_containers() -> None:
    """dereference_schema() strips out now-unreachable '$defs'/'definitions' containers."""
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "properties": {"x": {"$ref": "#/$defs/X"}},
        "$defs": {"X": {"type": "string"}},
    }
    result = dereference_schema(schema)
    assert "$defs" not in result


@pytest.mark.parametrize(
    ("path", "expected"),
    [
        (["properties", "a"], {"type": "string"}),
        (["items"], {"type": "string"}),
    ],
)
def test_dereference_schema_pass_resolves_refs_in_nested_positions(path: list[str], expected: dict[str, Any]) -> None:
    """dereference_schema() recursively resolves $refs nested inside properties/items."""
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "properties": {"a": {"$ref": "#/$defs/S"}},
        "items": {"$ref": "#/$defs/S"},
        "$defs": {"S": {"type": "string"}},
    }
    result = dereference_schema(schema)
    node: Any = result
    for key in path:
        node = node[key]
    assert node == expected


def test_dereference_schema_pass_resolves_refs_within_tuple_style_items_list() -> None:
    """dereference_schema() resolves $refs nested inside a Draft-07-style tuple 'items' list."""
    schema = {
        "$schema": "http://json-schema.org/draft-07/schema#",
        "type": "array",
        "items": [{"$ref": "#/definitions/S"}, {"type": "integer"}],
        "definitions": {"S": {"type": "string"}},
    }
    result = dereference_schema(schema)
    assert result["items"] == [{"type": "string"}, {"type": "integer"}]


def test_dereference_schema_pass_ref_siblings_ignored_pre_2019_09(caplog: pytest.LogCaptureFixture) -> None:
    """dereference_schema() ignores keywords declared alongside $ref under Draft-07 semantics by default."""
    schema = {
        "$schema": "http://json-schema.org/draft-07/schema#",
        "type": "object",
        "properties": {"a": {"$ref": "#/definitions/S", "description": "ignored"}},
        "definitions": {"S": {"type": "string"}},
    }
    result = dereference_schema(schema)
    assert "description" not in result["properties"]["a"]
    assert "Ignoring keyword" in caplog.text


def test_dereference_schema_pass_ref_siblings_applied_2019_09_plus() -> None:
    """dereference_schema() applies sibling keywords next to $ref under Draft 2019-09+/'auto' mode."""
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "properties": {"a": {"$ref": "#/$defs/S", "description": "kept"}},
        "$defs": {"S": {"type": "string"}},
    }
    result = dereference_schema(schema)
    assert result["properties"]["a"]["description"] == "kept"


def test_dereference_schema_pass_ref_sibling_mode_always_forces_merge() -> None:
    """dereference_schema(ref_sibling_mode='always') merges siblings even under Draft-07."""
    schema = {
        "$schema": "http://json-schema.org/draft-07/schema#",
        "type": "object",
        "properties": {"a": {"$ref": "#/definitions/S", "description": "kept"}},
        "definitions": {"S": {"type": "string"}},
    }
    result = dereference_schema(schema, ref_sibling_mode="always")
    assert result["properties"]["a"]["description"] == "kept"


def test_dereference_schema_pass_ref_sibling_mode_never_ignores_regardless_of_draft() -> None:
    """dereference_schema(ref_sibling_mode='never') ignores siblings even under Draft 2020-12."""
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "properties": {"a": {"$ref": "#/$defs/S", "description": "ignored"}},
        "$defs": {"S": {"type": "string"}},
    }
    result = dereference_schema(schema, ref_sibling_mode="never")
    assert "description" not in result["properties"]["a"]


@pytest.mark.parametrize(
    "keyword",
    [
        "additionalProperties",
        "additionalItems",
        "contains",
        "propertyNames",
        "contentSchema",
        "not",
        "if",
        "then",
        "else",
        "unevaluatedItems",
        "unevaluatedProperties",
    ],
)
def test_dereference_schema_pass_resolves_ref_within_each_single_schema_keyword(keyword: str) -> None:
    """dereference_schema() recursively resolves $refs nested inside each single-schema keyword."""
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        keyword: {"$ref": "#/$defs/S"},
        "$defs": {"S": {"type": "string"}},
    }
    result = dereference_schema(schema)
    assert result[keyword] == {"type": "string"}


@pytest.mark.parametrize("value", [True, False])
def test_dereference_schema_pass_boolean_valued_single_schema_keyword_passthrough(value: bool) -> None:
    """dereference_schema() passes a boolean (true/false) value of a single-schema keyword through unchanged."""
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "additionalProperties": value,
    }
    result = dereference_schema(schema)
    assert result["additionalProperties"] is value


def test_dereference_schema_pass_resolves_ref_within_dependent_schemas() -> None:
    """dereference_schema() recursively resolves $refs nested inside 'dependentSchemas' values."""
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "dependentSchemas": {"a": {"$ref": "#/$defs/S"}},
        "$defs": {"S": {"type": "string"}},
    }
    result = dereference_schema(schema)
    assert result["dependentSchemas"]["a"] == {"type": "string"}


def test_dereference_schema_pass_non_schema_keywords_pass_through_unchanged() -> None:
    """dereference_schema() leaves non-schema-valued keywords (annotations, assertions, dependentRequired) untouched."""
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "title": "My Schema",
        "pattern": "^a",
        "minLength": 1,
        "dependentRequired": {"a": ["b"]},
    }
    result = dereference_schema(schema)
    assert result["title"] == "My Schema"
    assert result["pattern"] == "^a"
    assert result["minLength"] == 1
    assert result["dependentRequired"] == {"a": ["b"]}


def test_dereference_schema_pass_unreferenced_defs_entry_dropped_without_resolution_even_if_broken() -> None:
    """dereference_schema() drops an unreferenced $defs entry entirely without erroring, even if it's internally broken."""
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "properties": {"a": {"type": "string"}},
        "$defs": {"Unused": {"$ref": "#/$defs/DoesNotExist"}},
    }
    result = dereference_schema(schema)  # must not raise
    assert "$defs" not in result


def test_dereference_schema_pass_all_of_branch_is_a_ref() -> None:
    """dereference_schema() resolves a $ref used as one of allOf's branches ('extends a base schema' pattern)."""
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "allOf": [
            {"$ref": "#/$defs/Base"},
            {"properties": {"extra": {"type": "string"}}},
        ],
        "$defs": {"Base": {"properties": {"id": {"type": "integer"}}, "required": ["id"]}},
    }
    result = dereference_schema(schema)
    assert set(result["properties"]) == {"id", "extra"}
    assert result["required"] == ["id"]


def test_dereference_schema_pass_nested_all_of_branches_merge_recursively() -> None:
    """dereference_schema() flattens an allOf branch that itself contains a nested allOf."""
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "allOf": [
            {"allOf": [{"properties": {"a": {"type": "string"}}}]},
            {"properties": {"b": {"type": "integer"}}},
        ],
    }
    result = dereference_schema(schema)
    assert "allOf" not in result
    assert set(result["properties"]) == {"a", "b"}


def test_merge_all_of_pass_later_branch_scalar_keyword_wins_over_earlier_branch() -> None:
    """_merge_all_of() lets a later allOf branch's scalar keyword value win over an earlier branch's."""
    schema = {"allOf": [{"title": "from first branch"}, {"title": "from second branch"}]}
    ctx = _make_ctx(schema)
    result = _merge_all_of(schema, ctx)
    assert result["title"] == "from second branch"


def test_dereference_schema_pass_empty_all_of_list_removes_keyword_with_no_merged_content() -> None:
    """dereference_schema() removes an empty 'allOf' list, leaving the rest of the schema untouched."""
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "allOf": [],
        "title": "Empty AllOf",
    }
    result = dereference_schema(schema)
    assert "allOf" not in result
    assert result["title"] == "Empty AllOf"


def test_dereference_schema_pass_repeated_sibling_ref_is_not_treated_as_circular() -> None:
    """dereference_schema() allows the same $ref to be used by multiple sibling properties without a false circular-ref error."""
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "properties": {
            "a": {"$ref": "#/$defs/S"},
            "b": {"$ref": "#/$defs/S"},
        },
        "$defs": {"S": {"type": "string"}},
    }
    result = dereference_schema(schema)
    assert result["properties"]["a"] == {"type": "string"}
    assert result["properties"]["b"] == {"type": "string"}


def test_dereference_schema_pass_ref_chain_exactly_at_max_depth_succeeds() -> None:
    """dereference_schema() succeeds when a $ref chain's length is exactly max_ref_depth."""
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "properties": {"a": {"$ref": "#/$defs/A"}},
        "$defs": {"A": {"$ref": "#/$defs/B"}, "B": {"type": "string"}},
    }
    result = dereference_schema(schema, max_ref_depth=2)
    assert result["properties"]["a"] == {"type": "string"}


def test_dereference_schema_fail_ref_chain_one_over_max_depth() -> None:
    """dereference_schema() raises when a $ref chain's length is exactly one more than max_ref_depth."""
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "properties": {"a": {"$ref": "#/$defs/A"}},
        "$defs": {"A": {"$ref": "#/$defs/B"}, "B": {"$ref": "#/$defs/C"}, "C": {"type": "string"}},
    }
    with pytest.raises(DereferencingError, match="exceeded max depth"):
        dereference_schema(schema, max_ref_depth=2)


def test_dereference_schema_pass_multi_hop_external_ref_chain() -> None:
    """dereference_schema() follows a $ref chain that hops from the root into one external resource, then into another."""
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "properties": {"a": {"$ref": "https://example.com/first.json"}},
    }
    resources = {
        "https://example.com/first.json": {
            "$schema": "https://json-schema.org/draft/2020-12/schema",
            "$ref": "https://example.com/second.json",
        },
        "https://example.com/second.json": {
            "$schema": "https://json-schema.org/draft/2020-12/schema",
            "type": "string",
        },
    }
    result = dereference_schema(schema, additional_resources=resources)
    assert result["properties"]["a"] == {"type": "string"}


def test_dereference_schema_pass_resolver_context_does_not_leak_between_sibling_external_refs() -> None:
    """dereference_schema() correctly restores resolver context between sibling properties referencing different external resources, and does not leak either resource's own '$schema' into the merged output."""
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "properties": {
            "a": {"$ref": "https://example.com/a.json"},
            "b": {"$ref": "https://example.com/b.json"},
        },
    }
    resources = {
        "https://example.com/a.json": {"$schema": "https://json-schema.org/draft/2020-12/schema", "type": "string"},
        "https://example.com/b.json": {"$schema": "https://json-schema.org/draft/2020-12/schema", "type": "integer"},
    }
    result = dereference_schema(schema, additional_resources=resources)
    assert result["properties"]["a"] == {"type": "string"}
    assert result["properties"]["b"] == {"type": "integer"}


def test_dereference_schema_pass_root_schema_keyword_preserved_after_identity_stripping_fix() -> None:
    """dereference_schema() still preserves the root document's own top-level '$schema' keyword."""
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "properties": {"a": {"$ref": "#/$defs/S"}},
        "$defs": {"S": {"type": "string"}},
    }
    result = dereference_schema(schema)
    assert result["$schema"] == "https://json-schema.org/draft/2020-12/schema"


def test_dereference_schema_pass_whole_document_ref_does_not_leak_target_schema_keyword() -> None:
    """dereference_schema() strips the target document's own '$schema' when a $ref points at an external document's root."""
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "properties": {"customer": {"$ref": "https://example.com/customer.json"}},
    }
    resources = {
        "https://example.com/customer.json": {
            "$schema": "https://json-schema.org/draft/2020-12/schema",
            "type": "object",
            "properties": {"name": {"type": "string"}},
        }
    }
    result = dereference_schema(schema, additional_resources=resources)
    assert "$schema" not in result["properties"]["customer"]


def test_dereference_schema_fail_circular_ref_through_all_of_branch() -> None:
    """dereference_schema() detects a circular $ref even when the cycle passes through an allOf branch."""
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "properties": {"node": {"$ref": "#/$defs/Node"}},
        "$defs": {"Node": {"allOf": [{"$ref": "#/$defs/Node"}]}},
    }
    with pytest.raises(DereferencingError, match="Circular \\$ref"):
        dereference_schema(schema)


def test_dereference_schema_fail_circular_reference() -> None:
    """dereference_schema() raises DereferencingError for a self-referencing $ref chain."""
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "properties": {"child": {"$ref": "#/$defs/Node"}},
        "$defs": {"Node": {"$ref": "#/$defs/Node"}},
    }
    with pytest.raises(DereferencingError, match="Circular \\$ref"):
        dereference_schema(schema)


def test_dereference_schema_fail_ref_to_boolean_schema() -> None:
    """dereference_schema() raises a clear DereferencingError (not a raw TypeError) when a $ref resolves to a boolean schema."""
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "properties": {"anything": {"$ref": "#/$defs/Anything"}},
        "$defs": {"Anything": True},
    }
    with pytest.raises(DereferencingError, match="boolean schema"):
        dereference_schema(schema)


def test_dereference_schema_fail_max_ref_depth_exceeded() -> None:
    """dereference_schema() raises DereferencingError when a $ref chain exceeds max_ref_depth."""
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "properties": {"a": {"$ref": "#/$defs/A"}},
        "$defs": {
            "A": {"$ref": "#/$defs/B"},
            "B": {"$ref": "#/$defs/C"},
            "C": {"type": "string"},
        },
    }
    with pytest.raises(DereferencingError, match="exceeded max depth"):
        dereference_schema(schema, max_ref_depth=1)


def test_dereference_schema_fail_unresolvable_ref() -> None:
    """dereference_schema() raises DereferencingError for a $ref pointing nowhere resolvable."""
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "properties": {"a": {"$ref": "#/$defs/DoesNotExist"}},
    }
    with pytest.raises(DereferencingError, match="Unable to resolve"):
        dereference_schema(schema)


def test_dereference_schema_pass_missing_schema_keyword_defaults_to_draft_2020_12_behavior() -> None:
    """dereference_schema() defaults to Draft 2020-12 semantics (e.g. ref-siblings applied) when '$schema' is absent."""
    schema = {
        "type": "object",
        "properties": {"a": {"$ref": "#/$defs/S", "description": "kept"}},
        "$defs": {"S": {"type": "string"}},
    }
    result = dereference_schema(schema)
    assert result["properties"]["a"]["description"] == "kept"


def test_dereference_schema_pass_resolver_dgaf_about_invalid_root_specification() -> None:
    """dereference_schema() soldiers on regardless if the root schema's dialect can't be determined."""
    schema = frozendict({"$schema": "not-a-real-dialect-uri", "type": "object"})
    deref_schema = dereference_schema({**schema})
    assert deref_schema == schema


def test_dereference_schema_pass_resolver_dgaf_about_additional_resource_specification() -> None:
    """dereference_schema() doesn't care if an external resource's dialect can't be determined."""
    schema = frozendict(
        {
            "$schema": "https://json-schema.org/draft/2020-12/schema",
            "type": "object",
        }
    )
    resources = frozendict({"https://example.com/bad.json": {"$schema": "not-a-real-dialect-uri"}})
    deref_schema = dereference_schema({**schema}, additional_resources={**resources})
    assert deref_schema == schema


"""_should_apply_ref_siblings"""


@pytest.mark.parametrize(
    ("mode", "has_siblings", "spec", "expected"),
    [
        ("never", True, DRAFT202012, False),
        ("never", False, DRAFT202012, False),
        ("always", True, DRAFT7, True),
        ("always", False, DRAFT7, False),
        ("auto", True, DRAFT202012, True),
        ("auto", True, DRAFT7, False),
        ("auto", False, DRAFT202012, False),
    ],
)
def test_should_apply_ref_siblings_pass_mode_and_spec_combinations(
    mode: RefSiblingMode, has_siblings: bool, spec: Specification, expected: bool
) -> None:
    """_should_apply_ref_siblings() honours ref_sibling_mode and falls back to per-draft rules under 'auto'."""
    ctx = _make_ctx({}, spec=spec, ref_sibling_mode=mode)
    assert _should_apply_ref_siblings(ctx, has_siblings) is expected


@pytest.mark.parametrize(
    ("mode", "has_siblings", "spec", "expected"),
    [
        ("auto", False, DRAFT7, False),
        ("auto", True, DRAFT201909, True),
    ],
)
def test_should_apply_ref_siblings_pass_additional_mode_and_spec_combinations(
    mode: RefSiblingMode,
    has_siblings: bool,
    spec: Specification,
    expected: bool,
) -> None:
    """_should_apply_ref_siblings() covers additional auto-mode combinations, including Draft 2019-09."""
    ctx = _make_ctx({}, spec=spec, ref_sibling_mode=mode)
    assert _should_apply_ref_siblings(ctx, has_siblings) is expected


"""_resolve_ref_and_all_of"""


def test_resolve_ref_and_all_of_pass_returns_schema_unchanged_when_no_ref_or_all_of() -> None:
    """_resolve_ref_and_all_of() is a no-op for a schema without $ref/allOf."""
    schema = {"type": "string", "minLength": 1}
    ctx = _make_ctx(schema)
    assert _resolve_ref_and_all_of(schema, ctx) == schema


def test_resolve_ref_and_all_of_fail_circular_ref_detected_directly() -> None:
    """_resolve_ref_and_all_of() raises DereferencingError as soon as a $ref reappears in active_refs."""
    root = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$defs": {"Node": {"$ref": "#/$defs/Node"}},
    }
    ctx = _make_ctx(root)
    with pytest.raises(DereferencingError, match="Circular \\$ref"):
        _resolve_ref_and_all_of({"$ref": "#/$defs/Node"}, ctx)


def test_resolve_ref_and_all_of_pass_strips_schema_keyword_from_resolved_ref_target() -> None:
    """_resolve_ref_and_all_of() strips a stray '$schema' declaration carried by the resolved $ref target."""
    root = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$defs": {"S": {"$schema": "https://json-schema.org/draft/2020-12/schema", "type": "string"}},
    }
    ctx = _make_ctx(root)
    result = _resolve_ref_and_all_of({"$ref": "#/$defs/S"}, ctx)
    assert result == {"type": "string"}


@pytest.mark.parametrize(
    ("identity_keyword", "value"),
    [
        ("$schema", "https://json-schema.org/draft/2020-12/schema"),
        ("$id", "https://example.com/foo"),
        ("$anchor", "Foo"),
    ],
)
def test_resolve_ref_and_all_of_pass_strips_various_identity_keywords_from_resolved_ref_target(
    identity_keyword: str, value: str
) -> None:
    """_resolve_ref_and_all_of() strips resource identity/dialect keywords ($schema/$id/$anchor) from $ref targets."""
    root = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$defs": {"S": {identity_keyword: value, "type": "string"}},
    }
    ctx = _make_ctx(root)
    result = _resolve_ref_and_all_of({"$ref": "#/$defs/S"}, ctx)
    assert identity_keyword not in result
    assert result["type"] == "string"


def test_resolve_ref_and_all_of_fail_ref_resolves_to_boolean_schema() -> None:
    """_resolve_ref_and_all_of() raises DereferencingError (not a raw TypeError/AttributeError) for a $ref to a boolean schema."""
    root = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$defs": {"Anything": True},
    }
    ctx = _make_ctx(root)
    with pytest.raises(DereferencingError, match="boolean schema"):
        _resolve_ref_and_all_of({"$ref": "#/$defs/Anything"}, ctx)


"""_merge_all_of"""


def test_merge_all_of_pass_merges_properties_required_and_pattern_properties() -> None:
    """_merge_all_of() unions 'properties'/'patternProperties' and de-duplicates 'required' across branches."""
    schema = {
        "allOf": [
            {"properties": {"a": {"type": "string"}}, "required": ["a"], "patternProperties": {"^x-": {}}},
            {"properties": {"b": {"type": "integer"}}, "required": ["a", "b"]},
        ]
    }
    ctx = _make_ctx(schema)
    result = _merge_all_of(schema, ctx)
    assert set(result["properties"]) == {"a", "b"}
    assert result["required"] == ["a", "b"]
    assert "^x-" in result["patternProperties"]


def test_merge_all_of_pass_sibling_keywords_take_final_precedence() -> None:
    """_merge_all_of() lets keywords declared alongside 'allOf' itself win over branch values."""
    schema = {
        "allOf": [{"title": "from branch"}],
        "title": "from sibling",
    }
    ctx = _make_ctx(schema)
    result = _merge_all_of(schema, ctx)
    assert result["title"] == "from sibling"


"""_dereference"""


@pytest.mark.parametrize("bool_schema", [True, False])
def test_dereference_pass_boolean_schema_passthrough(bool_schema: bool) -> None:
    """_dereference() passes boolean schemas (true/false) through unchanged."""
    ctx = _make_ctx({})
    assert _dereference(bool_schema, ctx) == bool_schema


def test_dereference_pass_resolves_ref_within_list_of_schemas_keyword() -> None:
    """_dereference() resolves $refs nested inside list-of-schemas keywords like 'anyOf'."""
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "anyOf": [{"$ref": "#/$defs/S"}],
        "$defs": {"S": {"type": "string"}},
    }
    ctx = _make_ctx(schema)
    result: dict[str, Any] = _dereference(schema, ctx)  # pyright: ignore[reportAssignmentType]
    assert result["anyOf"] == [{"type": "string"}]


def test_dereference_pass_dependencies_keyword_handles_mixed_schema_and_list_values() -> None:
    """_dereference() resolves schema-valued 'dependencies' entries but leaves list-valued entries untouched."""
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "dependencies": {
            "a": {"$ref": "#/$defs/S"},
            "b": ["c", "d"],
        },
        "$defs": {"S": {"type": "string"}},
    }
    ctx = _make_ctx(schema)
    result: dict[str, Any] = _dereference(schema, ctx)  # pyright: ignore[reportAssignmentType]
    assert result["dependencies"]["a"] == {"type": "string"}
    assert result["dependencies"]["b"] == ["c", "d"]
