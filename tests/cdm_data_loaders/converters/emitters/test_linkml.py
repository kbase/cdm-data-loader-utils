"""LinkML schema emission from typed converter documents."""

from pathlib import Path

import pytest
from linkml.linter.linter import Linter
from linkml_runtime.dumpers import yaml_dumper
from linkml_runtime.linkml_model.meta import ClassDefinition, EnumDefinition, SchemaDefinition, SlotDefinition
from linkml_runtime.loaders import yaml_loader

from cdm_data_loaders.converters.core.ir import Field, SchemaDocument, TypedNode
from cdm_data_loaders.converters.emitters.linkml import LinkMLEmitter, LinkMLEmitterError
from cdm_data_loaders.converters.readers.json_schema import JsonSchemaReader


def test_emit_pass_nested_enum_array_constraints() -> None:
    """Objects, enums, arrays, descriptions, and scalar constraints emit LinkML definitions."""
    document = JsonSchemaReader().read(
        {
            "$schema": "https://json-schema.org/draft/2020-12/schema",
            "$id": "https://example.org/record",
            "title": "record",
            "description": "A record.",
            "type": "object",
            "properties": {
                "id": {"type": "integer", "minimum": 1, "maximum": 10},
                "status": {"enum": ["active", "retired"]},
                "tags": {"type": "array", "items": {"type": "string", "pattern": "^[a-z]+$"}},
                "address": {
                    "type": "object",
                    "description": "Mailing address.",
                    "properties": {"city": {"type": "string"}},
                    "required": ["city"],
                },
            },
            "required": ["id", "status"],
        }
    )

    result = LinkMLEmitter().emit(document)
    assert isinstance(result, SchemaDefinition)
    assert isinstance(result.classes["record"], ClassDefinition)
    assert isinstance(result.classes["record"].attributes["id"], SlotDefinition)
    assert isinstance(result.enums["record_status_enum"], EnumDefinition)
    assert result == SchemaDefinition(
        id="https://example.org/record",
        name="record",
        prefixes={"record": "https://example.org/record#", "linkml": "https://w3id.org/linkml/"},
        default_prefix="record",
        imports=["linkml:types"],
        description="A record.",
        classes={
            "record": {
                "attributes": {
                    "id": {"range": "integer", "minimum_value": 1, "maximum_value": 10, "required": True},
                    "status": {"range": "record_status_enum", "required": True},
                    "tags": {"range": "string", "pattern": "^[a-z]+$", "multivalued": True},
                    "address": {"range": "record_address", "description": "Mailing address."},
                }
            },
            "record_address": {
                "attributes": {"city": {"range": "string", "required": True}},
                "description": "Mailing address.",
            },
        },
        enums={"record_status_enum": {"permissible_values": {"active": None, "retired": None}}},
    )
    assert yaml_loader.loads(yaml_dumper.dumps(result), target_class=SchemaDefinition) == result


def test_emit_pass_linkml_metamodel_validation(tmp_path: Path) -> None:
    """Emitted schemas conform to the LinkML metamodel."""
    document = JsonSchemaReader().read(
        {
            "$schema": "https://json-schema.org/draft/2020-12/schema",
            "type": "object",
            "properties": {
                "id": {"type": "integer", "minimum": 1},
                "status": {"enum": ["active", "retired"]},
            },
            "required": ["id"],
        }
    )
    schema_path = tmp_path / "record.yaml"
    schema_path.write_text(yaml_dumper.dumps(LinkMLEmitter().emit(document)), encoding="utf-8")

    assert list(Linter.validate_schema(str(schema_path))) == []


@pytest.mark.parametrize(
    "property_schema",
    [
        {"type": ["string", "integer"]},
        {"anyOf": [{"type": "string"}, {"type": "integer"}]},
        {"oneOf": [{"type": "string"}, {"type": "integer"}]},
        {"type": "array", "items": [{"type": "string"}, {"type": "integer"}]},
        {"enum": [1, 2]},
        {"type": "null"},
        {"type": "string", "minLength": 1},
    ],
    ids=["union-type", "any-of", "one-of", "tuple-array", "numeric-enum", "null", "unsupported-constraint"],
)
def test_emit_fail_unsupported_node(property_schema: dict[str, object]) -> None:
    """Unsupported JSON Schema constructs fail rather than silently changing their meaning."""
    document = JsonSchemaReader().read(
        {
            "$schema": "https://json-schema.org/draft/2020-12/schema",
            "type": "object",
            "properties": {"value": property_schema},
        }
    )

    with pytest.raises(LinkMLEmitterError):
        LinkMLEmitter().emit(document)


def test_emit_fail_non_object_root() -> None:
    """Reject documents without an object root."""
    with pytest.raises(LinkMLEmitterError, match="requires an object root"):
        LinkMLEmitter().emit(SchemaDocument(root=TypedNode(type="string")))


@pytest.mark.parametrize(
    "node",
    [
        TypedNode(type="map"),
        TypedNode(type="map", key_type=TypedNode(type="string")),
        TypedNode(type="map", value_type=TypedNode(type="integer")),
        TypedNode(type="string", constraints={"enum": "active"}),
    ],
    ids=["missing-map-types", "missing-map-value", "missing-map-key", "non-sequence-enum"],
)
def test_emit_node_fail_incomplete_schema(node: TypedNode) -> None:
    """Reject incomplete map definitions and non-sequence enum constraints."""
    with pytest.raises(LinkMLEmitterError):
        LinkMLEmitter().emit_node(node)


@pytest.mark.parametrize(
    ("title", "override", "expected_name", "expected_class"),
    [
        (None, None, "schema", "Root"),
        ("123 report", None, "_123_report", "_123_report"),
        ("!!!", None, "schema", "schema"),
        ("record", "custom schema", "custom_schema", "record"),
    ],
    ids=["unnamed", "leading-digit", "punctuation-only", "schema-name-override"],
)
def test_emit_pass_names(title: str | None, override: str | None, expected_name: str, expected_class: str) -> None:
    """Normalize schema and class names while retaining valid empty schemas."""
    document = SchemaDocument(root=TypedNode(type="object"), annotations={"title": title})
    result = LinkMLEmitter(schema_name=override).emit(document)

    assert result.name == expected_name
    assert result.id == f"urn:linkml:{expected_name}"
    assert result.default_prefix == expected_name
    assert result.classes == {expected_class: ClassDefinition(name=expected_class)}
    assert result.enums == {}


def test_emit_pass_unique_names_and_reuse() -> None:
    """Repeated emission keeps previous models intact and resets generated names."""
    document = SchemaDocument(
        name="record",
        root=TypedNode(
            type="object",
            properties=(
                Field("a-b", TypedNode(type="object")),
                Field("a_b", TypedNode(type="object")),
                Field("status", TypedNode(type="string", constraints={"enum": ("active", "retired")})),
            ),
        ),
    )
    emitter = LinkMLEmitter()
    first = emitter.emit(document)
    original_yaml = yaml_dumper.dumps(first)
    second = emitter.emit(SchemaDocument(name="other", root=TypedNode(type="object")))

    assert set(first.classes) == {"record", "record_a_b", "record_a_b_2"}
    assert first.classes["record"].attributes["a-b"].range == "record_a_b"
    assert first.classes["record"].attributes["a_b"].range == "record_a_b_2"
    assert set(first.enums) == {"record_status_enum"}
    assert second.classes == {"other": ClassDefinition(name="other")}
    assert second.enums == {}
    assert yaml_dumper.dumps(first) == original_yaml
    assert emitter.emit(document) == first


def test_emit_field_pass_description_precedence() -> None:
    """Field descriptions override node descriptions in official slot models."""
    result = LinkMLEmitter().emit_field(
        Field(
            "label",
            TypedNode(type="string", annotations={"description": "Node description."}),
            annotations={"description": "Field description."},
        )
    )
    assert result == SlotDefinition(name="label", range="string", description="Field description.")


def test_emit_into_pass_unique_names_and_descriptions() -> None:
    """Appending documents preserves classes, enums, descriptions, and collision-free references."""
    document = SchemaDocument(
        name="record",
        annotations={"description": "Table description."},
        root=TypedNode(
            type="object",
            properties=(
                Field("child", TypedNode(type="object")),
                Field("status", TypedNode(type="string", constraints={"enum": ("active",)})),
            ),
        ),
    )
    emitter = LinkMLEmitter()
    schema = SchemaDefinition(id="urn:test:namespace", name="namespace")
    assert emitter.emit_into(document, schema) is schema
    assert emitter.emit_into(document, schema) is schema
    assert schema.classes == {
        name: ClassDefinition(
            name=name,
            description="Table description.",
            attributes=[
                SlotDefinition(name="child", range=f"{name}_child"),
                SlotDefinition(name="status", range=f"{name}_status_enum"),
            ],
        )
        for name in ("record", "record_2")
    } | {name: ClassDefinition(name=name) for name in ("record_child", "record_2_child")}
    assert schema.enums == {
        name: EnumDefinition(name=name, permissible_values={"active": None})
        for name in ("record_status_enum", "record_2_status_enum")
    }


@pytest.mark.parametrize(
    "root",
    [
        TypedNode(type="string"),
        TypedNode(
            type="object",
            properties=(
                Field("child", TypedNode(type="object")),
                Field("status", TypedNode(type="string", constraints={"enum": ("active",)})),
                Field("unsupported", TypedNode(type="null")),
            ),
        ),
    ],
    ids=["non-object-root", "failure-after-nested-definitions"],
)
def test_emit_into_fail_preserves_destination(root: TypedNode) -> None:
    """Failed emission leaves the destination unchanged and permits a subsequent conversion."""
    emitter = LinkMLEmitter()
    schema = emitter.emit(SchemaDocument(name="existing", root=TypedNode(type="object")))
    original_yaml = yaml_dumper.dumps(schema)
    with pytest.raises(LinkMLEmitterError):
        emitter.emit_into(SchemaDocument(name="broken", root=root), schema)
    assert yaml_dumper.dumps(schema) == original_yaml
    emitter.emit_into(SchemaDocument(name="good", root=TypedNode(type="object")), schema)
    assert schema.classes == {
        "existing": ClassDefinition(name="existing"),
        "good": ClassDefinition(name="good"),
    }
    assert schema.enums == {}
