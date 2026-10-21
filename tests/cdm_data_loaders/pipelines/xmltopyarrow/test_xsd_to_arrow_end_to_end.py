"""End-to-end specification tests: XSD + XML -> process_xml_file_to_pyarrow -> PyArrow Table batches."""

import datetime as dt
import json
from collections.abc import Callable
from pathlib import Path

import pyarrow as pa
import pytest

from cdm_data_loaders.pipelines.xmltopyarrow.settings import XmlToPyArrowSettings
from cdm_data_loaders.readers.xml import process_xml_file_to_pyarrow
from cdm_data_loaders.readers.xsd import generate_pyarrow_schema, load_schema


def _json_default(obj: object) -> str:
    """Serialize Python date and bytes values from Arrow tables for comparison with expected JSON."""
    if isinstance(obj, dt.date):
        return obj.isoformat()
    if isinstance(obj, bytes):
        return obj.decode("utf-8")
    msg = f"Object of type {type(obj)} is not JSON serializable"
    raise TypeError(msg)


def _strip_none(value: object) -> object:
    """Recursively remove keys with None values from dicts to compare Arrow-materialized rows with sparse expected.json."""
    if isinstance(value, dict):
        return {k: _strip_none(v) for k, v in value.items() if v is not None}
    if isinstance(value, list):
        return [_strip_none(item) for item in value]
    return value


@pytest.mark.parametrize(
    ("xsd_rel_path", "xml_rel_path", "expected_rel_path", "xml_tag"),
    [
        pytest.param(
            "xsd/spec/01_scalars_and_simple_types/schema.xsd",
            "xsd/spec/01_scalars_and_simple_types/data.xml",
            "xsd/spec/01_scalars_and_simple_types/expected.json",
            "record",
            id="spec-01-scalars-and-simple-types",
        ),
        pytest.param(
            "xsd/spec/02_namespaces_and_nillable/schema.xsd",
            "xsd/spec/02_namespaces_and_nillable/data.xml",
            "xsd/spec/02_namespaces_and_nillable/expected.json",
            "entry",
            id="spec-02-namespaces-and-nillable",
        ),
        pytest.param(
            "xsd/spec/03_mixed_content_raw_text/schema.xsd",
            "xsd/spec/03_mixed_content_raw_text/data.xml",
            "xsd/spec/03_mixed_content_raw_text/expected.json",
            "article",
            id="spec-03-mixed-content-raw-text",
        ),
        pytest.param(
            "xsd/spec/04_recursion_depth/schema.xsd",
            "xsd/spec/04_recursion_depth/data.xml",
            "xsd/spec/04_recursion_depth/expected.json",
            "record",
            id="spec-04-recursion-depth",
        ),
        pytest.param(
            "xsd/spec/05_groups_and_compositors/schema.xsd",
            "xsd/spec/05_groups_and_compositors/data.xml",
            "xsd/spec/05_groups_and_compositors/expected.json",
            "item",
            id="spec-05-groups-and-compositors",
        ),
        pytest.param(
            "xsd/spec/06_context_dependent_paths/schema.xsd",
            "xsd/spec/06_context_dependent_paths/data.xml",
            "xsd/spec/06_context_dependent_paths/expected.json",
            "root",
            id="spec-06-context-dependent-paths",
        ),
        pytest.param(
            "xsd/spec/07_multi_root_shared_type/schema.xsd",
            "xsd/spec/07_multi_root_shared_type/data.xml",
            "xsd/spec/07_multi_root_shared_type/expected.json",
            "entry",
            id="spec-07-multi-root-shared-type",
        ),
        pytest.param(
            "xsd/spec/08_nullability_and_empty_tags/schema.xsd",
            "xsd/spec/08_nullability_and_empty_tags/data.xml",
            "xsd/spec/08_nullability_and_empty_tags/expected.json",
            "entry",
            id="spec-08-nullability-and-empty-tags",
        ),
        pytest.param(
            "xsd/spec/09_empty_marker_elements/schema.xsd",
            "xsd/spec/09_empty_marker_elements/data.xml",
            "xsd/spec/09_empty_marker_elements/expected.json",
            "book",
            id="spec-09-empty-marker-elements",
        ),
        pytest.param(
            "xsd/spec/10_combinations/schema.xsd",
            "xsd/spec/10_combinations/data.xml",
            "xsd/spec/10_combinations/expected.json",
            "record",
            id="spec-10-combinations",
        ),
        pytest.param(
            "uniprot/uniref/uniref.xsd",
            "xsd/examples/uniref/minimal.xml",
            "xsd/examples/uniref/minimal_expected.json",
            "entry",
            id="example-uniref-minimal",
        ),
        pytest.param(
            "uniprot/uniref/uniref.xsd",
            "xsd/examples/uniref/maximal.xml",
            "xsd/examples/uniref/maximal_expected.json",
            "entry",
            id="example-uniref-maximal",
        ),
        pytest.param(
            "uniprot/uniprot_kb/uniprot.xsd",
            "xsd/examples/uniprot/minimal.xml",
            "xsd/examples/uniprot/minimal_expected.json",
            "entry",
            id="example-uniprot-minimal",
        ),
        pytest.param(
            "uniprot/uniprot_kb/uniprot.xsd",
            "xsd/examples/uniprot/maximal.xml",
            "xsd/examples/uniprot/maximal_expected.json",
            "entry",
            id="example-uniprot-maximal",
        ),
        pytest.param(
            "pubmed/pubmed.xsd",
            "xsd/examples/pubmed_article/minimal.xml",
            "xsd/examples/pubmed_article/minimal_expected.json",
            "PubmedArticle",
            id="example-pubmed-article-minimal",
        ),
        pytest.param(
            "pubmed/pubmed.xsd",
            "xsd/examples/pubmed_article/maximal.xml",
            "xsd/examples/pubmed_article/maximal_expected.json",
            "PubmedArticle",
            id="example-pubmed-article-maximal",
        ),
        pytest.param(
            "pubmed/pubmed.xsd",
            "xsd/examples/pubmed_book_article/minimal.xml",
            "xsd/examples/pubmed_book_article/minimal_expected.json",
            "PubmedBookArticle",
            id="example-pubmed-book-article-minimal",
        ),
        pytest.param(
            "pubmed/pubmed.xsd",
            "xsd/examples/pubmed_book_article/maximal.xml",
            "xsd/examples/pubmed_book_article/maximal_expected.json",
            "PubmedBookArticle",
            id="example-pubmed-book-article-maximal",
        ),
    ],
)
def test_process_xml_file_to_pyarrow_pass_yields_arrow_tables_for_spec_and_examples(
    settings_factory: Callable[..., XmlToPyArrowSettings],
    test_data_dir: Path,
    *,
    xsd_rel_path: str,
    xml_rel_path: str,
    expected_rel_path: str,
    xml_tag: str,
) -> None:
    """process_xml_file_to_pyarrow yields pa.Table batches matching the XSD PyArrow schema and expected JSON."""
    xsd_path = test_data_dir / xsd_rel_path
    xml_path = test_data_dir / xml_rel_path
    expected = json.loads((test_data_dir / expected_rel_path).read_text(encoding="utf-8"))

    settings = settings_factory(
        table_name=xml_tag,
        xml_tag=xml_tag,
        xsd_file=str(xsd_path),
    )

    items = list(process_xml_file_to_pyarrow(settings, file_path=xml_path))
    assert len(items) == 1
    assert items[0].meta.table_name == xml_tag
    assert isinstance(items[0].data, pa.Table)

    table: pa.Table = items[0].data
    expected_pa_schema = generate_pyarrow_schema(load_schema(xsd_path), xml_tag, max_recursion_depth=2)
    assert table.schema.equals(expected_pa_schema, check_metadata=True)
    assert table.num_rows == len(expected)

    actual_rows = json.loads(json.dumps(_strip_none(table.to_pylist()), default=_json_default))
    assert actual_rows == expected
