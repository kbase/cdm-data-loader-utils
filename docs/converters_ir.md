# Converter IR and Reader Contract

The converters package translates schema definitions between formats. Every conversion is a
reader (source -> IR) followed by an emitter (IR -> target), composed by a thin facade that
owns the public entry points, input guards and target policy defaults. The intermediate
representation (IR) is a frozen node tree that preserves source facts without choosing target
policy; emitters make every target-specific decision. Captured legacy outputs in the golden
file are the compatibility reference for the four original routes, not the definition of
source truth.

## Package Layout

All paths below are under `src/cdm_data_loaders/converters/`.

```
core/                          shared, target-independent foundation
  errors.py                    ConversionError base class
  ir.py                        NodeType, Source, Provenance, NodeHints, TypedNode,
                               Field, SchemaDocument, Reader/Emitter protocols
  ir_values.py                 Value, freeze_value, freeze_mapping, mutable_value
  extensions.py                ExtensionError, ExtensionSpec, ExtensionRegistry,
                               Extensions, DEFAULT_EXTENSIONS, validated_extensions
  guards.py                    require_schema_keyword, require_object_root,
                               reject_unresolved_references
  inference.py                 infer_implicit_type, json_type_from_enum, is_unconstrained_node,
                               decimal_places, decimal_pattern, resolve_union_type
  metadata.py                  ConversionContext, metadata_keys_for, get_known_jsonschema_keywords,
                               DEFAULT_METADATA_KEYWORDS, REF_AND_IDENTITY_KEYWORDS,
                               STRUCTURAL_OR_COMPOSITIONAL_KEYWORDS
  io.py                        load_schema_text, load_schema_file
  paths.py                     NestedPathError, set_nested
readers/                       source -> IR
  json_schema.py               JsonSchemaReader
  dlt.py                       DltReader
  iceberg.py                   IcebergReader, iceberg_value
emitters/                      IR -> target
  json_schema.py               JsonSchemaEmitter
  dlt.py                       DltEmitter
  pyspark.py                   PySparkEmitter, format-map helpers
  pyarrow.py                   PyArrowEmitter, format-map helpers
  linkml.py                    LinkMLEmitter
jsonschema_to_dlt.py           facade: JSON Schema -> dlt stored schema
jsonschema_to_pyspark.py       facade: JSON Schema -> PySpark StructType
jsonschema_to_pyarrow.py       facade: JSON Schema -> PyArrow Schema / empty Table
dlt_to_jsonschema.py           facade: dlt -> JSON Schema documents
pyiceberg_to_jsonschema.py     functions over IcebergReader + JsonSchemaEmitter
pyiceberg_to_linkml.py         functions over IcebergReader + LinkMLEmitter
pyiceberg_catalog_to_jsonschema.py  CLI: dump catalog tables as JSON Schema
pyiceberg_catalog_to_linkml.py      CLI: dump catalog tables as LinkML
```

`core/` imports no dlt, PyIceberg or PySpark modules. The only cross-package import is the
`x-xsv-config` contract, which reuses `X_XSV_CONFIG_SCHEMA` from
`cdm_data_loaders.readers.jsonschema_xsv.xsv_validator.custom_metaschema`. dlt table
normalization lives outside the package in `cdm_data_loaders.utils.dlt.dlt_normalization`
(`DltTableNode`, `unflatten_tables`, `flatten_nodes`, `tables_parents_first`, `child_key`,
`child_table_name`).

| Module                    | Public API                                                                                                         |
| ------------------------- | ------------------------------------------------------------------------------------------------------------------ |
| `core/ir.py`              | `NodeType`, `Source`, `Provenance`, `NodeHints`, `TypedNode`, `Field`, `SchemaDocument`, `Reader`, `Emitter`       |
| `core/ir_values.py`       | `Value`, `freeze_value`, `freeze_mapping`, `mutable_value`                                                         |
| `core/extensions.py`      | `ExtensionError`, `ExtensionSpec`, `ExtensionRegistry`, `Extensions`, `DEFAULT_EXTENSIONS`, `validated_extensions` |
| `core/guards.py`          | `require_schema_keyword`, `require_object_root`, `reject_unresolved_references`                                    |
| `core/inference.py`       | `infer_implicit_type`, `json_type_from_enum`, `decimal_places`, `decimal_pattern`, `resolve_union_type`            |
| `core/io.py`              | `load_schema_text`, `load_schema_file`                                                                             |
| `core/paths.py`           | `set_nested`                                                                                                       |
| `readers/json_schema.py`  | `JsonSchemaReader`, including `read_node`                                                                          |
| `readers/dlt.py`          | `DltReader`                                                                                                        |
| `readers/iceberg.py`      | `IcebergReader`, including `read_node` and `read_field`; `iceberg_value`                                           |
| `emitters/json_schema.py` | `JsonSchemaEmitter`, including `emit_node` and `emit_field`                                                        |
| `emitters/dlt.py`         | `DltEmitter`                                                                                                       |
| `core/metadata.py`        | `ConversionContext`, `metadata_keys_for`, `get_known_jsonschema_keywords`, metadata keyword constants              |
| `emitters/pyspark.py`     | `PySparkEmitter`, `merge_format_map`                                                                               |
| `emitters/pyarrow.py`     | `PyArrowEmitter`, `merge_format_map`, `infer_type_from_enum`                                                       |
| `emitters/linkml.py`      | `LinkMLEmitter`, including `emit_node` and `emit_field`                                                            |

The IR, value helpers and extension registry import no dlt, PyIceberg or PySpark modules.
The XSV contract imports only the existing custom metaschema module, not Spark readers.
Models are frozen, slotted dataclasses. Constructors validate and recursively copy mappings
into `frozendict` and sequences into tuples. There is no unchecked construction path.

### LinkML Output

`LinkMLEmitter.emit()` and `table_to_linkml()` return official `linkml-runtime`
`SchemaDefinition` objects, not dictionaries. `emit_node()`, `emit_field()`,
`convert_type()`, `convert_field()`, and `convert_struct()` return `SlotDefinition`
objects. Access model attributes directly, for example `schema.classes["record"].attributes["id"].range`.

Serialize these models using `linkml_runtime.dumpers.yaml_dumper.dumps(schema)`.
The catalog CLI writes UTF-8 YAML to `*.schema.yaml`, with one UTC generation
timestamp per run in `annotations.generated_at.value`. It does not emit JSON or
the nonstandard top-level `x-iceberg` property.

By default, the CLI writes one `<namespace>.<table>.schema.yaml` file per table.
Set `group_by_namespace=True` (CLI: `--group-by-namespace true`, environment:
`CDL_GROUP_BY_NAMESPACE=true`) to write one `<namespace>.schema.yaml` file instead.
Each table becomes a class with local attributes; nested structs and maps retain
their supporting classes. Table comments become class descriptions. Name collisions
receive numeric suffixes, with tables processed in sorted identifier order.
Empty namespaces and namespaces where every table fails produce no file.
Failed tables are logged and skipped without adding partial definitions.

For in-memory grouping, `table_to_linkml(table, identifier, schema=destination)`
appends to an existing `SchemaDefinition` and returns it. `LinkMLEmitter.emit_into()`
provides the same operation for typed schema documents.

### PyArrow Output

`PyArrowEmitter.emit()` returns a `pyarrow.Schema`; `emit_table()` returns an empty `pyarrow.Table`.
`JSONSchemaToPyArrow` is the facade: `convert()` and `convert_from_string/file()` return the schema,
`convert_to_table()` returns the empty table. The root must resolve to a struct.

JSON dispatch follows the PySpark emitter: declared or inferred types precede combiners, `oneOf`
then `anyOf` approximate their first branch, and `treat_unknown_as_string` controls the string
fallback. Differences from Spark:

- Integers: `int32` when both bounds fit, otherwise `int64`.
- Numbers: `decimal128(38, scale)` from `multipleOf`, otherwise `float64`.
- `date` is `date32`; `date-time` is `timestamp[us, tz=UTC]`; `time` stays a string.
- Arrays are `list<item>` and dynamic objects are `map<string, value>`. List elements and map values
  are nullable, except where a typed source (Iceberg, dlt) declares them required.
- Typed sources use their hints: bit widths select `int8` to `int64`, `float` is `float32`,
  decimals above precision 38 (and dlt `wei`) are `decimal256`, `timestamp` precision selects the
  unit (0: `s`, 1-3: `ms`, 4-6: `us`, 7-9: `ns`), `time` is `time64[us]`, and `fixed` is
  `fixed_size_binary(length)`. A missing timestamp timezone is treated as UTC.
- Field metadata uses Arrow string pairs: `jsonschema` (JSON object of selected keywords, decimals
  encoded as strings), `comment` (description) and `original_union` (JSON array).

## Models

The following are complete constructor fields. `Mapping` preserves insertion order.
All metadata mappings default to empty immutable mappings; all fields not marked optional
are required constructor arguments. Metadata values retain absent keys separately from null.

```python
type Value = None | bool | int | float | Decimal | str | tuple[Value, ...] | Mapping[str, Value]
type Source = Literal["json-schema", "dlt", "iceberg"]
type NodeType = Literal[
    "object", "map", "array", "string", "integer", "number", "boolean",
    "null", "any", "never", "unknown",
]

Provenance(
    source: Source,
    path: tuple[str | int, ...] = (),
    metadata: Mapping[str, Value] = {},
)

NodeHints(
    logical_type: str | None = None,
    precision: int | None = None,
    scale: int | None = None,
    bit_width: int | None = None,
    timezone: bool | None = None,
    length: int | None = None,
)

TypedNode(
    type: NodeType,
    nullable: bool | None = None,
    declared_type: str | tuple[str, ...] | None = None,
    inferred_type: str | None = None,
    hints: NodeHints = NodeHints(),
    properties: tuple[Field, ...] = (),
    required_names: tuple[str, ...] | None = None,
    items: TypedNode | tuple[TypedNode, ...] | None = None,
    prefix_items: tuple[TypedNode, ...] | None = None,
    pattern_properties: Mapping[str, TypedNode] = {},
    additional_properties: TypedNode | None = None,
    any_of: tuple[TypedNode, ...] | None = None,
    one_of: tuple[TypedNode, ...] | None = None,
    schema_keywords: Mapping[str, TypedNode] = {},
    schema_maps: Mapping[str, Mapping[str, TypedNode]] = {},
    key_type: TypedNode | None = None,
    value_type: TypedNode | None = None,
    constraints: Mapping[str, Value] = {},
    annotations: Mapping[str, Value] = {},
    extensions: Mapping[str, Value] = {},
    provenance: Provenance | None = None,
    source_keywords: tuple[str, ...] = (),
)

Field(
    name: str,
    node: TypedNode,
    required: bool = False,
    annotations: Mapping[str, Value] = {},
    extensions: Mapping[str, Value] = {},
    provenance: Provenance | None = None,
)

SchemaDocument(
    root: TypedNode,
    name: str | None = None,
    dialect: str | None = None,
    identifier: str | None = None,
    annotations: Mapping[str, Value] = {},
    extensions: Mapping[str, Value] = {},
    provenance: Provenance | None = None,
    jsonschema_validator_cls: type | None = None,
)
```

`JsonSchemaReader.read()` detects the draft and records its validator class in
`jsonschema_validator_cls`, preserving the original URI in `dialect`. The PySpark emitter
reuses this class for metadata selection. Documents without a detected validator and
standalone nodes without an explicit conversion context use Draft 2020-12 metadata rules.

`Field.required` describes property presence, not whether its value accepts null.
`TypedNode.nullable=None` means unresolved. JSON unions are not flattened or selected.
`declared_type=None` means absent; a scalar string and a tuple union declaration remain
distinct, including a one-element union. `inferred_type` records the shared keyword inference
without replacing the declaration or choosing a combiner branch. Enum/const values remain
in `constraints`; inference from them is an emitter choice. `unknown` is not `string`.
Boolean schema `true` is `any`; `false` is `never`; explicit `type: null` is `null`.

`required_names` retains original order and names without declared properties. For dlt it
also retains legacy required-list order, including duplicate names when a child overwrites
a same-named column. `source_keywords` retains JSON keyword presence and order, distinguishing
an omitted `properties` or `type` from an explicit empty declaration. Boolean schemas have
no keywords. Target dispatch must consider these facts, not just `node.type`.

Schema structure is separate from literal values:

- `properties`, `pattern_properties`, `items`, `prefix_items`, `additional_properties`,
  `any_of` and `one_of` hold typed children, preserving declaration order.
- `schema_keywords` holds `additionalItems`, `unevaluatedItems`, `contains`,
  `unevaluatedProperties`, `propertyNames`, `not`, `if`, `then`, `else`, `contentSchema`.
- `schema_maps` holds `$defs`, `definitions`, `dependentSchemas`, and schema-valued
  `dependencies`. Array-valued `dependencies` stays in `constraints`.
- `constraints` retains all remaining non-extension keywords, including exact `multipleOf`,
  min/max and exclusive bounds, enum, const, pattern, format, content encoding/media type,
  lengths, item counts, uniqueness, and unknown readable keywords.
- `annotations` holds `$schema`, `$id`, `$anchor`, `$dynamicAnchor`, `$comment`, `$vocabulary`,
  title, description, default, examples, readOnly, writeOnly and deprecated.

Decimals remain finite `Decimal` objects, with no conversion through binary floats.
`freeze_value(value: object) -> Value` and
`freeze_mapping(value: Mapping[str, object]) -> Mapping[str, Value]` reject arbitrary objects,
non-string mapping keys and non-finite numbers. `mutable_value(value: Value) -> object`
returns independent dict/list containers and preserves Decimal. This is not a JSON encoder:
future serializers must encode Decimal deliberately and must not convert it to float.

## Extension Contracts

```python
ExtensionSpec(namespace: str, schema: Mapping[str, Value] | bool)
ExtensionRegistry(specs: tuple[ExtensionSpec, ...] = ())
ExtensionRegistry.register(spec: ExtensionSpec, *, replace: bool = False) -> ExtensionRegistry
ExtensionRegistry.extend_payload(namespace: str, properties: Mapping[str, object]) -> ExtensionRegistry
ExtensionRegistry.validate(values: Mapping[str, object]) -> Mapping[str, Value]
Extensions(payload: Mapping[str, Value] = {}, registry: ExtensionRegistry = DEFAULT_EXTENSIONS)
validated_extensions(values: Mapping[str, Value]) -> Extensions
```

Only top-level namespace names require a nonempty `x-` prefix. Internal names such as
`data_type` and `field_id` are ordinary keys. Key spelling is never rewritten with `replace`.
Each namespace has an explicit draft 2020-12 validation schema. Arrays and null are valid
when that schema permits them. Normal IR construction wraps plain mappings in `Extensions`
using `DEFAULT_EXTENSIONS`; a prevalidated `Extensions` instance retains its custom registry.
Schemas, additional property contracts and payloads are recursively copied and immutable.
Registration returns a new registry. Duplicate namespaces need `replace=True`; extending
a payload rejects key collisions. Ordinary payload objects reject unknown keys.

Builtins:

- `x-dlt`: name, data_type, nullable, description, precision, scale, primary_key, unique,
  foreign_key, sort, cluster, partition, merge_key, row_key, root_key, variant and timezone.
  Boolean hints permit null to preserve observed metadata. Precision/scale are nonnegative
  integers or null. Data type names are strings or null so unknown source types remain readable.
  Custom hints, including internal `x-*` hints, require `extend_payload("x-dlt", {...})`.
- `x-iceberg`: field_id, required, initial_default, write_default, logical_type, precision,
  scale, length, element_id, element_required, key_id, value_id, value_required, identifier,
  schema_id, format_version, current_snapshot_id, location, properties, partition_spec,
  identifier_field_ids and generated_at. Defaults allow arbitrary IR values; other keys
  have explicit scalar, array or closed-object contracts.
- `x-xsv-config`: exactly `X_XSV_CONFIG_SCHEMA` from the existing custom metaschema,
  including its nonempty recognized-key requirement and closed property set. This does
  not impose a nonempty root JSON Schema `required` array.
- `x-file-glob`, `x-dlt-prefix`, `x-dlt-split`: strings; `x-delimiter`: one character;
  `x-pii`: boolean. These match observed repository tests.

Unregistered extensions raise `ExtensionError(ConversionError)`. This is intentionally
stricter than the pre-phase-5 facades: requesting `extra_metadata_keywords={"x-custom"}`
selects metadata but does not register its schema. All facade constructors accept an actual
immutable `extension_registry` instance; arbitrary mappings are not coerced into registries.
The Iceberg functions accept the registry as an optional keyword-only argument.
Forward and reverse dlt facade errors retain their public exception classes and chain
reader/emitter failures, including `ExtensionError` and its underlying validation error.
`preserve_unknown_hints=False` filters emitted hints, never bypasses reader validation.

Explicit registration example:

```python
from cdm_data_loaders.converters.core.extensions import DEFAULT_EXTENSIONS, ExtensionSpec
from cdm_data_loaders.converters.jsonschema_to_pyspark.converter import JSONSchemaToPySpark

registry = DEFAULT_EXTENSIONS.register(
    ExtensionSpec(
        "x-vendor",
        {
            "type": "object",
            "properties": {"classification": {"enum": ["public", "internal"]}},
            "required": ["classification"],
            "additionalProperties": False,
        },
    )
)
converter = JSONSchemaToPySpark(
    extension_registry=registry,
    extra_metadata_keywords=frozenset({"x-vendor"}),
)
schema = converter.convert(
    {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "properties": {
            "name": {
                "type": "string",
                "x-vendor": {"classification": "public"},
            }
        },
    }
)
```

Custom column hints require `registry.extend_payload("x-dlt", {"x-owner":
{"type": "string", "enum": ["lab"]}})`. Unknown payload keys and invalid values still fail.
Registration returns a new registry and never mutates the builtin contracts. `x-pii` remains
boolean by default. A legacy metadata-selection test used `"some-value"` for `x-pii`; that
test now explicitly replaces the namespace contract with a string schema, retaining its
original output assertion. This is a documented validation tightening, not inferred permission.

## Readers

```python
JsonSchemaReader(
    extension_registry: ExtensionRegistry = DEFAULT_EXTENSIONS,
    error_type: type[ConversionError] = ConversionError,
    converter_name: str = "JsonSchemaReader",
    dereference_function: str = "dereference_schema",
)
JsonSchemaReader.read(source: Mapping[str, Any]) -> SchemaDocument
JsonSchemaReader.read_node(source: Mapping[str, Any] | bool) -> TypedNode

DltReader(
    include_dlt_columns: bool = False,
    include_variant_columns: bool = False,
    child_table_mode: Literal["object", "array"] = "object",
    extension_registry: ExtensionRegistry = DEFAULT_EXTENSIONS,
)
DltReader.read(source: Mapping[str, Any]) -> dict[str, SchemaDocument]

IcebergReader(extension_registry: ExtensionRegistry = DEFAULT_EXTENSIONS, mapping_target: str = "IR")
IcebergReader.read(source: Schema) -> SchemaDocument
IcebergReader.read_table(table: Table, identifier: tuple[str, ...]) -> SchemaDocument
IcebergReader.read_node(source: IcebergType) -> TypedNode
IcebergReader.read_field(source: NestedField) -> Field
iceberg_value(value: object) -> Value

JsonSchemaEmitter.emit_node(node: TypedNode) -> dict[str, Any] | bool
JsonSchemaEmitter.emit_field(field: Field) -> dict[str, Any] | bool
PySparkEmitter.emit_node(node: TypedNode, ctx: ConversionContext | None = None) -> DataType
PyArrowEmitter.emit_node(node: TypedNode, ctx: ConversionContext | None = None) -> pa.DataType

class Reader[Input, Output](Protocol):
    def read(self, source: Input) -> Output: ...

class Emitter[Output](Protocol):
    def emit(self, document: SchemaDocument) -> Output: ...
```

### JSON Schema

Input is already validated and dereferenced. The reader requires `$schema` and an
object-compatible root: absent type, `object`, or a type list containing `object`.
It rejects `$ref`, unmerged `allOf`, `$dynamicRef`, and `$recursiveRef` only at schema
positions. Literal enum/const/default/example and extension data are never scanned as
schemas. It does not run a generic root metaschema validation or perform dereferencing.
Nested boolean schemas and both tuple forms are retained. No `treat_unknown_as_string`
reader setting exists. Root annotations/extensions are mirrored into the document envelope;
emitters emit them once, while retaining root type constraints.

Facades apply their original stricter root guard first: declared root type lists, including
`["object", "null"]`, remain rejected with the original converter-specific message.
Missing `$schema` raises the facade's own `InvalidJSONSchemaError`. The reader's diagnostic
options preserve `$ref`/`allOf` messages while checking every schema-valued branch, including
ones that old target dispatch ignored. Newly rejected nested references are intentional
defensive validation; literal values and extension payloads are not scanned as schemas.

### dlt

`unflatten_tables` owns graph validation and reconstruction, including multiple roots,
explicit parents, cycles and ordering. Empty tables are rejected. Each root becomes a
document; document provenance contains the stored-schema envelope and reader policy.
Table and column provenance contain immutable source metadata, including filtered columns.
Emitters must traverse actual nodes, not re-run an old converter over this metadata.

Filtering uses `_dlt_` and `__v_` names, matching the legacy converter, before scalar
classification and required propagation. A child with exactly one visible `value` column
and no children becomes an array of that typed value. Other children use the configured
object/array mode, recorded in document provenance. A child is required when at least one
direct visible column is nonnullable; grandchildren alone do not make it required.

Incomplete and `json` columns become `any`, retaining their nullable fact. Unknown type names
become `unknown` and log, never explicit strings. Decimal is numeric with precision/scale;
decimal-as-string and decimal patterns belong to the JSON emitter. `bigint.precision`
is retained as bit width. A missing timestamp timezone hint remains unknown, not naive.
All source column keys are validated under `x-dlt`, not silently discarded.

### Iceberg

Primitive widths, logical types, decimal precision/scale, fixed length, field IDs and
list/map element IDs are retained. Maps have distinct `key_type` and `value_type` nodes;
they are not JSON dynamic objects or prematurely rendered key/value arrays.
Timestamp and timestamptz are `timestamp-without-tz`/`False` and
`timestamp-with-tz`/`True`. Binary is readable even though the old JSON converter rejects it.
Unmapped types, including current nanosecond/unknown types, raise `NotImplementedError`.

Field descriptions/defaults/IDs live on `Field`; type facts live on `TypedNode`; table
metadata lives on `SchemaDocument`. This separation supports legacy nullable wrapper
placement. `read_table` reads only loaded table accessors, capturing snapshot ID, partition
spec, properties, location, format version, schema ID and identifier IDs without I/O.

`iceberg_value` explicitly converts UUID to text, temporal objects to ISO text and bytes to
base64; finite Decimal stays exact. Other unsupported Python objects fail. PyIceberg's
`NestedField.__init__` currently populates both defaults and doc even when omitted. The
reader preserves those observable null values and cannot recover presence already lost
by PyIceberg. JSON and dlt metadata preserve absent versus null directly.

## Phase 4 Emitter Rules

Consume `document.root`, fields, logical hints, constraints and scoped extensions directly.
Do not reconstruct source schemas as an intermediate call into old monoliths.

1. JSON-to-dlt: implement dlt's existing dispatch order from declared/inferred/enum facts,
   with combiners overriding declared types where the legacy converter does so. Both
   forward converters prefer `oneOf` over `anyOf` and approximate the first branch.
   Keep nesting limits, scalar flattening and JSON fallback in the emitter.
2. JSON-to-PySpark: recognized declared or inferred types precede combiners. Retain field
   presence separately from nullable values; apply strict/lenient unknown handling there.
   Use constraints and annotations for formats, numeric widths and field metadata.
3. dlt-to-JSON: apply existing logical-type templates, decimal string pattern, precision/
   scale and hint-selection policy. Put column descriptions and `x-dlt` inside nullable
   wrappers. Empty `any` nodes without emitted metadata must stay `{}`. Use `required_names`
   for exact legacy required-list order. Root titles/IDs are derived from document names;
   table descriptions are retained on root/child nodes.
4. Iceberg-to-JSON: derive integer bounds, formats, decimal string patterns and map-as-pairs
   encoding from logical hints. Type metadata goes inside optional `anyOf`; field metadata
   and field descriptions go outside it. Nullable list elements and map values use their
   own wrappers. Legacy mode omits null defaults, suppresses element IDs it never emitted,
   and preserves its unsupported-type errors even where IR can represent more.

Identity claims are limited to supported behavior. JSON constraints preserved by the reader
need not survive a legacy target conversion. dlt child-table ambiguity is irrecoverable.
An IR-aware future mode may improve output, but compatibility mode must match captured
structured outputs, not claim arbitrary JSON/YAML byte identity.

## Facades and Composition

The four original converter modules contain no recursive conversion traversal or target
type dispatch. `JSONSchemaToDlt.to_schema` and `to_yaml` delegate to `DltEmitter`; both forward
`convert_from_string` methods remain JSON-only. Reverse dlt conversion still accepts JSON
or YAML text and maps root documents through `JsonSchemaEmitter`. Unknown dlt type errors
retain `Column 'name' has unknown dlt data_type 'type'.` and the old catch identity.

Iceberg `convert_type`, `convert_field`, `convert_struct`, and `table_to_json_schema` call
fragment/document APIs directly. `decimal_pattern` is an alias to the JSON emitter function.
The catalog converter still stamps `x-iceberg.generated_at` after conversion. PySpark's
tested private methods are short reader/emitter bridges; metadata constants, context and
helpers are owned by the emitter and exported at the original paths. The public dlt format
constants also remain aliases derived from the emitter's single format map. Dead untested private
traversals and the Iceberg dispatch map were removed, not retained as a superclass.

Direct dlt-to-PySpark composition needs no JSON intermediate or Spark session:

```python
from cdm_data_loaders.converters.readers.dlt import DltReader
from cdm_data_loaders.converters.emitters.pyspark import PySparkEmitter

documents = DltReader().read(
    {
        "records": {
            "columns": {
                "amount": {"data_type": "decimal", "precision": 20, "scale": 8, "nullable": False},
            }
        }
    }
)
spark_schema = PySparkEmitter().emit(documents["records"])
```

This produces `decimal(20,8)` directly. A dlt-to-JSON route intentionally produces a decimal
string pattern instead. JSON -> dlt -> JSON tests assert the supported required-scalar
subset; they do not claim general lossless conversion of constraints, array origins or metadata.

Logging is configured through standard module loggers, without a callback:

```python
import logging

logging.getLogger("cdm_data_loaders.converters.emitters.dlt").setLevel(logging.WARNING)
logging.getLogger("cdm_data_loaders.converters.emitters.pyspark").setLevel(logging.WARNING)
logging.getLogger("cdm_data_loaders.converters.readers.dlt").setLevel(logging.WARNING)
```

Warnings now belong to the reader/emitter that makes the decision, not the facade module.
Logger-name assertions were updated accordingly; type/output assertions were preserved.

## Tests and Limits

`tests/data/converters/ir/legacy_outputs.json` contains seven fixed cases with inputs,
options and mechanically captured outputs from all four pre-rewrite converter routes.
Spark outputs use `StructType.jsonValue()` without a Spark session. Other outputs are
stored dlt dictionaries or JSON Schema documents. Cases cover nested trees, scalar arrays,
unions, precedence, dynamic objects, metadata, multi-root dlt, and Iceberg maps/defaults.
The existing precise converter tests remain primary; the corpus is not exhaustive.

Tests cover readers, emitters, the existing facade entry points, the core IR modules
(`core/ir.py`, `core/ir_values.py`, `core/extensions.py`), and the per-facade and
end-to-end test modules. Emitter comparisons that became tautologies were replaced
with independent explicit expectations. The captured golden file is unchanged. No internal
mocks, external services or Spark sessions are used by these tests.

The validation commands are:

```sh
uv run pytest tests/cdm_data_loaders/converters -m 'not requires_spark and not requires_ceph and not external_request' -o addopts='' -o log_cli=false -q
```

Coverage adds `--cov=cdm_data_loaders.converters.<module>` for the four facade modules,
`readers` and `emitters`, with `--cov-branch --cov-report=term-missing`. The terminal-only
report leaves the repository's existing coverage XML unchanged.

Readers recurse over typed children; unlike the iterative dlt normalizer, they do not
promise support beyond Python's recursion limit. JSON input validation is upstream.
Extension validation uses registered JSON Schemas; caller-provided schemas are trusted
configuration and should be self-contained. Emitters retain exact Decimal objects; a public
Decimal JSON serializer is not provided. These are implementation boundaries, not
output-equivalence claims.
