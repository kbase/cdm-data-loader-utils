"""Resolve all $ref references in a schema to create a fully standalone JSONSchema document.

Fully resolves `$ref` — including references into external JSON Schema
documents supplied via `additional_resources` — and merges `allOf`
throughout a JSON Schema document, producing an equivalent, self-contained
schema with no remaining `$ref`/`allOf` keywords anywhere in the tree.

`jsonschema_to_pyspark.converter.JSONSchemaToPySpark.convert()` requires its
input to already be in this dereferenced form; run it through
`dereference_schema()` first if your schema uses `$ref` (local or external)
or `allOf`.
"""

import logging
from collections.abc import Mapping
from dataclasses import dataclass
from dataclasses import field as dataclass_field
from typing import Any, Final, Literal

from referencing import Registry, Resource
from referencing._core import Resolver
from referencing.exceptions import CannotDetermineSpecification, Unresolvable
from referencing.jsonschema import DRAFT201909, DRAFT202012, specification_with

logger = logging.getLogger(__name__)

RefSiblingMode = Literal["auto", "always", "never"]

# Drafts under which keywords declared alongside `$ref` are meaningful and applied
REF_SIBLING_AWARE_SPECS = (DRAFT202012, DRAFT201909)

# Keywords whose values are {name: schema} maps.
DICT_OF_SCHEMAS_KEYWORDS: frozenset[str] = frozenset({"properties", "patternProperties", "dependentSchemas"})
# Keywords whose values are lists of schemas.
LIST_OF_SCHEMAS_KEYWORDS: frozenset[str] = frozenset({"anyOf", "oneOf", "prefixItems"})
# Keywords whose value is a single schema (or boolean schema).
SINGLE_SCHEMA_KEYWORDS: frozenset[str] = frozenset(
    {
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
    }
)
# Pure definition containers: only ever reached via `$ref`, which is already
# resolved/inlined by the time we get here, so they don't need to be carried
# forward into the dereferenced output.
DECLARATION_ONLY_KEYWORDS: frozenset[str] = frozenset({"$defs", "definitions"})

# Resource identity/dialect keywords. These describe a schema *resource*
# (which document it is, which draft it's written against) rather than a
# validation constraint on instance data.
REF_TARGET_IDENTITY_KEYWORDS: Final[frozenset[str]] = frozenset(
    {"$anchor", "$dynamicAnchor", "$id", "$schema", "$vocabulary", "id"}
)


class DereferencingError(ValueError):
    """Raised when a JSON Schema document cannot be fully dereferenced."""


@dataclass
class DereferenceContext:
    """Mutable state for a single `dereference_schema()` call.

    :param resolver: the `referencing` resolver rooted at the schema currently being dereferenced
    :type resolver: Resolver
    :param spec: the JSON Schema specification/dialect in effect
    :type spec: Any
    :param max_ref_depth: maximum depth of a `$ref` chain before giving up
    :type max_ref_depth: int
    :param ref_sibling_mode: how to handle keywords declared alongside `$ref`
    :type ref_sibling_mode: RefSiblingMode
    :param active_refs: the set of `$ref` values currently being resolved, used for circular-reference detection
    :type active_refs: set[str]
    """

    resolver: Resolver
    spec: Any
    max_ref_depth: int
    ref_sibling_mode: RefSiblingMode
    active_refs: set[str] = dataclass_field(default_factory=set)


def dereference_schema(
    schema: dict[str, Any],
    additional_resources: Mapping[str, dict[str, Any]] | None = None,
    *,
    max_ref_depth: int = 50,
    ref_sibling_mode: RefSiblingMode = "auto",
) -> dict[str, Any]:
    """Resolve every `$ref` in a schema using `additional_resources` where needed.

    Merges any `allOf` tags to return a fully self-contained schema.

    :param schema: the JSON Schema document to dereference
    :type schema: dict[str, Any]
    :param additional_resources: mapping of URI to schema dict, used to resolve external `$ref`s
    :type additional_resources: Mapping[str, dict[str, Any]] | None
    :param max_ref_depth: maximum depth of a `$ref` chain before giving up
    :type max_ref_depth: int
    :param ref_sibling_mode: how to handle keywords declared alongside `$ref`:
        'auto' (spec-compliant per detected draft), 'always' (permissive), or
        'never' (strict pre-2019-09 semantics)
    :type ref_sibling_mode: RefSiblingMode
    :return: an equivalent, `$ref`/`allOf`-free JSON Schema document
    :rtype: dict[str, Any]
    :raises DereferencingError: if a `$ref` is unresolvable, circular, or exceeds `max_ref_depth`,
        or if a dialect cannot be determined for the root schema or an external resource
    """
    additional_resources = additional_resources or {}

    dialect_id = schema.get("$schema")
    spec = specification_with(dialect_id, default=DRAFT202012) if dialect_id else DRAFT202012

    try:
        root_resource = Resource.from_contents(schema, default_specification=spec)
    except CannotDetermineSpecification as exc:
        err_msg = f"Could not determine a JSON Schema dialect for the root schema: {exc}"
        raise DereferencingError(err_msg) from exc

    try:
        resources = [
            (uri, Resource.from_contents(doc, default_specification=spec)) for uri, doc in additional_resources.items()
        ]
    except CannotDetermineSpecification as exc:
        err_msg = f"Could not determine a JSON Schema dialect for an external resource in `additional_resources`: {exc}"
        raise DereferencingError(err_msg) from exc

    registry: Registry = Registry().with_resources(resources)
    resolver: Resolver[Any] = registry.resolver_with_root(root_resource)

    ctx = DereferenceContext(
        resolver=resolver,
        spec=spec,
        max_ref_depth=max_ref_depth,
        ref_sibling_mode=ref_sibling_mode,
    )
    return _dereference(schema, ctx)


def _should_apply_ref_siblings(ctx: DereferenceContext, has_siblings: bool) -> bool:  # noqa: FBT001
    """Decide whether keywords declared alongside a `$ref` should be applied.

    :param ctx: the active dereference context
    :type ctx: DereferenceContext
    :param has_siblings: whether the schema declares any keywords alongside `$ref`
    :type has_siblings: bool
    :return: True if sibling keywords should be merged into the resolved schema
    :rtype: bool
    """
    if not has_siblings or ctx.ref_sibling_mode == "never":
        return False
    if ctx.ref_sibling_mode == "always":
        return True
    # ref_sibling_mode is "auto"
    return ctx.spec in REF_SIBLING_AWARE_SPECS


def _resolve_ref_and_all_of(schema: dict[str, Any], ctx: DereferenceContext) -> dict[str, Any]:
    """Resolve `allOf` and `$ref` references at a single schema level.

    Does *not* descend into nested sub-schemas (e.g. `properties` values) — see `_dereference` instead.

    :param schema: the schema fragment to resolve at this level
    :type schema: dict[str, Any]
    :param ctx: the active dereference context
    :type ctx: DereferenceContext
    :return: the schema with `$ref`/`allOf` resolved at this level (not recursively)
    :rtype: dict[str, Any]
    :raises DereferencingError: if a `$ref` is unresolvable, circular, or the chain
        exceeds `ctx.max_ref_depth`
    """
    if "$ref" in schema:
        ref = schema["$ref"]

        if len(ctx.active_refs) >= ctx.max_ref_depth:
            msg = f"$ref chain exceeded max depth ({ctx.max_ref_depth}) near {ref!r}."
            raise DereferencingError(msg)

        if ref in ctx.active_refs:
            err_msg = (
                f"Circular $ref detected at {ref!r}. PySpark schemas must be finite; JSON Schemas with "
                "self-referential structures (e.g. recursive tree/linked-list) cannot be dereferenced."
            )
            raise DereferencingError(err_msg)

        try:
            resolved = ctx.resolver.lookup(ref)
        except Unresolvable as exc:
            err_msg = f"Unable to resolve $ref {ref!r}: {exc}"
            raise DereferencingError(err_msg) from exc

        if isinstance(resolved.contents, bool):
            err_msg = (
                f"$ref {ref!r} resolves to a boolean schema ({resolved.contents!r}). Boolean schemas "
                "reached via $ref are not currently supported by the dereferencer; inline the boolean "
                "schema directly at the reference site instead."
            )
            raise DereferencingError(err_msg)

        sibling_keys = {k: v for k, v in schema.items() if k != "$ref"}
        apply_siblings = _should_apply_ref_siblings(ctx, bool(sibling_keys))

        # Strip resource identity/dialect keywords from the *resolved* content before
        # merging it in -- see `REF_TARGET_IDENTITY_KEYWORDS` docstring above.
        merged = {k: v for k, v in resolved.contents.items() if k not in REF_TARGET_IDENTITY_KEYWORDS}
        if apply_siblings:
            merged.update(sibling_keys)
        elif sibling_keys:
            logger.warning(
                "Ignoring keyword(s) %s declared alongside $ref %r, per JSON Schema draft semantics "
                "(siblings of $ref are not applied in this draft). Pass ref_sibling_mode='always' to "
                "override.",
                sorted(sibling_keys),
                ref,
            )

        previous_resolver = ctx.resolver
        ctx.resolver = resolved.resolver
        ctx.active_refs.add(ref)
        try:
            return _resolve_ref_and_all_of(merged, ctx)
        finally:
            ctx.active_refs.discard(ref)
            ctx.resolver = previous_resolver

    if "allOf" in schema:
        return _merge_all_of(schema, ctx)

    return schema


def _merge_all_of(schema: dict[str, Any], ctx: DereferenceContext) -> dict[str, Any]:
    """Merge an `allOf` keyword's branches into a single flattened schema.

    `properties`/`required`/`patternProperties` are unioned across all `allOf`
    branches and any sibling declaration of the same keywords; other sibling
    keywords take final precedence over their `allOf`-branch counterparts.

    :param schema: the schema fragment containing the `allOf` keyword to merge
    :type schema: dict[str, Any]
    :param ctx: the active dereference context
    :type ctx: DereferenceContext
    :return: the schema with `allOf` merged away
    :rtype: dict[str, Any]
    """
    merged: dict[str, Any] = {}
    merged_properties: dict[str, Any] = {}
    merged_required: list[str] = []
    merged_pattern_properties: dict[str, Any] = {}

    for sub_schema in schema["allOf"]:
        resolved_sub = _resolve_ref_and_all_of(sub_schema, ctx)
        merged_properties.update(resolved_sub.get("properties", {}))
        merged_required.extend(resolved_sub.get("required", []))
        merged_pattern_properties.update(resolved_sub.get("patternProperties", {}))
        merged.update(resolved_sub)  # last-write-wins for scalar keywords

    # Keywords declared alongside `allOf` itself take final precedence, but
    # properties/required/patternProperties still union with the allOf branches
    # rather than replacing them outright.
    sibling_schema = {k: v for k, v in schema.items() if k != "allOf"}
    merged_properties.update(sibling_schema.get("properties", {}))
    merged_required.extend(sibling_schema.get("required", []))
    merged_pattern_properties.update(sibling_schema.get("patternProperties", {}))
    merged.update(sibling_schema)

    if merged_properties:
        merged["properties"] = merged_properties
    if merged_required:
        merged["required"] = list(dict.fromkeys(merged_required))  # de-dup, keep order
    if merged_pattern_properties:
        merged["patternProperties"] = merged_pattern_properties

    return merged


def _dereference(schema: dict[str, Any] | bool, ctx: DereferenceContext) -> dict[str, Any] | bool:  # noqa: FBT001
    """Recursively resolve `$ref`/`allOf` at every level of `schema`.

    :param schema: the schema fragment (or boolean schema) to dereference
    :type schema: dict[str, Any] | bool
    :param ctx: the active dereference context
    :type ctx: DereferenceContext
    :return: the fully dereferenced schema fragment, with nested schemas (in `properties`,
        `items`, `anyOf`, etc.) recursively dereferenced as well
    :rtype: dict[str, Any] | bool
    :raises DereferencingError: if any nested `$ref` is unresolvable, circular, or exceeds
        `ctx.max_ref_depth`
    """
    if isinstance(schema, bool):
        return schema

    resolved = _resolve_ref_and_all_of(schema, ctx)

    result: dict[str, Any] = {}
    for key, value in resolved.items():
        if key in DECLARATION_ONLY_KEYWORDS:
            continue
        if key in DICT_OF_SCHEMAS_KEYWORDS and isinstance(value, dict):
            result[key] = {k: _dereference(v, ctx) for k, v in value.items()}
        elif key in LIST_OF_SCHEMAS_KEYWORDS and isinstance(value, list):
            result[key] = [_dereference(v, ctx) for v in value]
        elif key == "items" and isinstance(value, list):
            # Draft-07-style tuple validation: list of per-position schemas.
            result[key] = [_dereference(v, ctx) for v in value]
        elif key == "dependencies" and isinstance(value, dict):
            # Legacy keyword: values are either a schema or a list of required property names.
            result[key] = {k: (_dereference(v, ctx) if isinstance(v, (dict, bool)) else v) for k, v in value.items()}
        elif key in SINGLE_SCHEMA_KEYWORDS or key == "items":
            result[key] = _dereference(value, ctx) if isinstance(value, (dict, bool)) else value
        else:
            result[key] = value
    return result
