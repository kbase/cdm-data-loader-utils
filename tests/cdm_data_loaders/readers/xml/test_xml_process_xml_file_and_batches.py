"""Tests for `process_xml_file`, `process_xml_file_items`, `build_xml_file_resource`, and `process_xml_file_to_dict` from cdm_data_loaders.readers.xml."""

import gzip
import logging
from collections.abc import Callable, Generator, Iterable, Iterator
from pathlib import Path
from typing import Any, Final
from unittest.mock import MagicMock

import pytest
import xmltodict
from dlt.common.storages.fsspec_filesystem import FileItemDict
from dlt.extract.items import DataItemWithMeta, TableNameMeta
from dlt.sources.filesystem import filesystem
from lxml.etree import Element, XMLSyntaxError, tostring

import cdm_data_loaders.readers.xml as xml_module
import cdm_data_loaders.utils.buffer as buffer_module
from cdm_data_loaders.core.fields import BUFFER_SIZE, DEFAULT_XML_FILE_GLOB, DEFAULTS, LOG_INTERVAL
from cdm_data_loaders.core.settings import CtsSettings
from cdm_data_loaders.pipelines.xmltodict.settings import XmlToDictSettings
from cdm_data_loaders.readers.xml import (
    DEFAULT_XMLTODICT_ARGS,
    process_xml_file,
    process_xml_file_batches,
    process_xml_file_items,
    process_xml_file_to_dict,
)

ParseFn = Callable[..., dict[str, list[dict[str, Any]]]]

UNIPROT_NS: Final[str] = "http://uniprot.org/uniprot"

# Arbitrary non-default value used to prove settings are forwarded to process_xml_file.
CUSTOM_LOG_INTERVAL: Final[int] = 7
CUSTOM_BATCH_LOG_INTERVAL: Final[int] = 50

PEOPLE_XML_2: Final[str] = """<?xml version="1.0" encoding="UTF-8"?>
<people>
    <person id="1">
        <name>Anne Example</name>
        <email>a@example.com</email>
    </person>
    <person id="2">
        <name>Belinda Carlisle</name>
        <email>b@example.com</email>
    </person>
</people>
"""

PEOPLE_XML_2_ONE_MISSING_EMAIL: Final[str] = """<?xml version="1.0" encoding="UTF-8"?>
<people>
    <person id="1">
        <name>Anne Example</name>
        <email>a@example.com</email>
    </person>
    <person id="2">
        <name>Belinda Carlisle</name>
    </person>
</people>
"""

PEOPLE_XML_5: Final[str] = """<?xml version="1.0" encoding="UTF-8"?>
<people>
    <person id="1"><name>Anne Example</name><email>a@example.com</email></person>
    <person id="2"><name>Belinda Carlisle</name><email>b@example.com</email></person>
    <person id="3"><name>Carol Singer</name><email>c@example.com</email></person>
    <person id="4"><name>Debbie Downer</name><email>d@example.com</email></person>
    <person id="5"><name>Ex Ample</name><email>e@example.com</email></person>
</people>
"""

PEOPLE_PARSED = [
    {"id": "1", "name": "Anne Example", "email": "a@example.com"},
    {"id": "2", "name": "Belinda Carlisle", "email": "b@example.com"},
    {"id": "3", "name": "Carol Singer", "email": "c@example.com"},
    {"id": "4", "name": "Debbie Downer", "email": "d@example.com"},
    {"id": "5", "name": "Ex Ample", "email": "e@example.com"},
]

# The same PEOPLE_XML_5 entries in the raw xmltodict shape process_xml_file_to_dict produces: one
# {xml_tag: {...}} dict per entry, with attr_prefix="_" turning the "id" attribute into "_id".
PEOPLE_XML_TO_DICT_PARSED = [
    {"person": {"_id": p["id"], "name": p["name"], "email": p["email"]}} for p in PEOPLE_PARSED
]


PEOPLE_XML_EMPTY: Final[str] = """<?xml version="1.0" encoding="UTF-8"?>
<people>
</people>
"""

PEOPLE_XML_MALFORMED: Final[str] = """<?xml version="1.0" encoding="UTF-8"?>
<people>
    <person id="1"><name>Anne Example</name></person>
    <person id="2"><name>Belinda Carlisle</name></persoo>
</people>
"""

UNIPROT_XML_2: Final[str] = f"""<?xml version="1.0" encoding="UTF-8"?>
<uniprot xmlns="{UNIPROT_NS}">
    <entry dataset="Swiss-Prot">
        <accession>P12345</accession>
        <name>TEST_HUMAN</name>
        <organism>
            <name type="scientific">Homo sapiens</name>
        </organism>
    </entry>
    <entry dataset="Swiss-Prot">
        <accession>Q67890</accession>
        <name>TEST2_HUMAN</name>
        <organism>
            <name type="scientific">Homo sapiens</name>
        </organism>
    </entry>
</uniprot>
"""


def _write_xml(directory: Path, filename: str, content: str, *, gzip_compress: bool = False) -> Path:
    """Write an XML string to `directory / filename`, optionally gzip-compressed."""
    file_path = directory / filename
    data = content.encode("utf-8")
    if gzip_compress:
        with gzip.open(file_path, "wb") as f:
            f.write(data)
    else:
        file_path.write_bytes(data)
    return file_path


def parse_person_entry(entry: Element, file_path: Path) -> dict[str, list[dict[str, Any]]]:
    """Reference parse_fn: turn a <person> element into a single "people" row."""
    return {
        "people": [
            {
                "id": entry.get("id"),
                "name": entry.findtext("name"),
                "email": entry.findtext("email"),
                "source_file": str(file_path),
            }
        ]
    }


def parse_person_multi_table(
    entry: Element,
    file_path: Path,  # noqa: ARG001
) -> dict[str, list[dict[str, Any]]]:
    """Reference parse_fn: emit both a "people" row and a related "emails" row per entry."""
    person_id = entry.get("id")
    return {
        "people": [{"id": person_id, "name": entry.findtext("name")}],
        "emails": [{"person_id": person_id, "email": entry.findtext("email")}],
    }


def parse_person_skip_missing_email(
    entry: Element,
    file_path: Path,  # noqa: ARG001
) -> dict[str, list[dict[str, Any]]]:
    """Reference parse_fn: return an empty dict (no rows) for entries missing an <email>."""
    email = entry.findtext("email")
    if email is None:
        return {}
    return {"people": [{"id": entry.get("id"), "email": email}]}


def parse_uniprot_entry(
    entry: Element,
    file_path: Path,  # noqa: ARG001
) -> dict[str, list[dict[str, Any]]]:
    """Reference parse_fn: extract accession/name from a namespaced UniProt <entry>."""
    ns = f"{{{UNIPROT_NS}}}"
    return {
        "proteins": [
            {
                "accession": entry.findtext(f"{ns}accession"),
                "name": entry.findtext(f"{ns}name"),
            }
        ]
    }


def parse_returns_non_dict(entry: Element, file_path: Path) -> list[str]:  # noqa: ARG001
    """Malformed parse_fn: returns a list instead of the required dict[str, rows]."""
    return [entry.get("id") or ""]


def make_failing_parse_fn(fail_on_call: int) -> ParseFn:
    """Build a parse_fn that raises ValueError on the N-th invocation (1-indexed)."""
    state = {"n": 0}

    def _parse_fn(entry: Element, file_path: Path) -> dict[str, list[dict[str, Any]]]:  # noqa: ARG001
        state["n"] += 1
        if state["n"] == fail_on_call:
            msg = f"boom on call {fail_on_call}"
            raise ValueError(msg)
        return {"people": [{"id": entry.get("id")}]}

    return _parse_fn


def _table_and_data(items: Iterable[DataItemWithMeta]) -> list[tuple[str, Any]]:
    """Flatten a sequence of DataItemWithMeta into (table_name, data) pairs for easy assertion."""
    result: list[tuple[str, Any]] = []
    for item in items:
        assert isinstance(item, DataItemWithMeta)
        assert isinstance(item.meta, TableNameMeta)
        result.append((item.meta.table_name, item.data))
    return result


def fake_settings(
    buffer_size: int = DEFAULTS[BUFFER_SIZE],
    log_interval: int = DEFAULTS[LOG_INTERVAL],
    settings_class: type = CtsSettings,
    input_dir: str | Path = ".",
) -> CtsSettings | XmlToDictSettings:
    """Build a MagicMock stand-in for CtsSettings, XmlToDictSettings, etc.

    ``process_xml_file`` reads ``buffer_size`` and ``log_interval`` off the settings
    object, so these are exposed as real integers (rather than nested mocks) to make
    the arithmetic in the function under test behave correctly. ``process_xml_file_items``
    also reads ``input_dir``.
    """
    settings = MagicMock(spec=settings_class)
    settings.buffer_size = buffer_size
    settings.log_interval = log_interval
    settings.input_dir = Path(input_dir)
    return settings


def fake_xmltodict_settings(
    table_name: str = "people",
    xml_tag: str = "person",
    buffer_size: int = DEFAULTS[BUFFER_SIZE],
    log_interval: int = DEFAULTS[LOG_INTERVAL],
) -> XmlToDictSettings:
    """Build a MagicMock XmlToDictSettings with table_name/xml_tag set, for process_xml_file_to_dict tests."""
    settings: XmlToDictSettings = fake_settings(  # pyright: ignore[reportAssignmentType]
        buffer_size=buffer_size, log_interval=log_interval, settings_class=XmlToDictSettings
    )
    settings.table_name = table_name
    settings.xml_tag = xml_tag
    return settings


def test_process_xml_file_pass_single_table_single_entry(tmp_path: Path) -> None:
    """Verify a single matching element is parsed into one correctly table-tagged row."""
    file_path = _write_xml(tmp_path, "one.xml", PEOPLE_XML_2)
    items = list(process_xml_file(fake_settings(), "person", parse_person_entry, file_path=file_path))  # pyright: ignore[reportArgumentType]
    assert len(items) == 1
    parsed = _table_and_data(items)
    assert parsed[0][0] == "people"
    assert parsed[0][1] == [{"source_file": str(file_path), **p_data} for p_data in PEOPLE_PARSED[:2]]


def test_process_xml_file_pass_multiple_entries(tmp_path: Path) -> None:
    """Verify all entries in a multi-entry XML file are streamed and parsed in document order."""
    file_path = _write_xml(tmp_path, "five.xml", PEOPLE_XML_5)
    items = list(process_xml_file(fake_settings(), "person", parse_person_entry, file_path=file_path))  # pyright: ignore[reportArgumentType]
    parsed = _table_and_data(items)
    assert parsed[0][0] == "people"
    assert parsed[0][1] == [{"source_file": str(file_path), **p_data} for p_data in PEOPLE_PARSED]


@pytest.mark.parametrize(
    ("buffer_size", "expected_batch_sizes"),
    [
        pytest.param(1, [1, 1, 1, 1, 1], id="one-row-batches"),
        pytest.param(2, [2, 2, 1], id="two-row-batches-with-remainder"),
        pytest.param(3, [3, 2], id="three-row-batches-with-remainder"),
        pytest.param(10, [5], id="larger-than-input"),
    ],
)
def test_process_xml_file_pass_buffer_size_controls_batching(
    tmp_path: Path, buffer_size: int, expected_batch_sizes: list[int]
) -> None:
    """Verify rows are yielded in batches of at most `buffer_size`, with a final partial batch.

    ``buffer_size`` is now supplied via the settings object and constrained to be a
    positive integer by ``CtsSettings``.
    """
    file_path = _write_xml(tmp_path, "five.xml", PEOPLE_XML_5)

    items = _table_and_data(
        process_xml_file(
            fake_settings(buffer_size=buffer_size),  # pyright: ignore[reportArgumentType]
            "person",
            parse_person_entry,
            file_path=file_path,
        )
    )

    assert [len(rows) for _, rows in items] == expected_batch_sizes
    assert [row["id"] for _, rows in items for row in rows] == [str(index) for index in range(1, 6)]


def test_process_xml_file_pass_multiple_tables_per_entry(tmp_path: Path) -> None:
    """Verify a parse_fn returning multiple tables yields one buffered item per table."""
    file_path = _write_xml(tmp_path, "two.xml", PEOPLE_XML_2)
    items = list(process_xml_file(fake_settings(), "person", parse_person_multi_table, file_path=file_path))  # pyright: ignore[reportArgumentType]
    tagged = _table_and_data(items)
    assert [table for table, _ in tagged] == ["people", "emails"]
    assert [row["id"] for row in tagged[0][1]] == ["1", "2"]
    assert [row["email"] for row in tagged[1][1]] == ["a@example.com", "b@example.com"]


def test_process_xml_file_pass_gzip_file(tmp_path: Path) -> None:
    """Verify gzip-compressed XML input (.gz suffix) is transparently decompressed and parsed."""
    plain_path = _write_xml(tmp_path, "plain.xml", PEOPLE_XML_2)
    gz_path = _write_xml(tmp_path, "compressed.xml.gz", PEOPLE_XML_2, gzip_compress=True)

    plain_items = _table_and_data(process_xml_file(fake_settings(), "person", parse_person_entry, file_path=plain_path))  # pyright: ignore[reportArgumentType]
    gz_items = _table_and_data(process_xml_file(fake_settings(), "person", parse_person_entry, file_path=gz_path))  # pyright: ignore[reportArgumentType]

    plain_ids = [row["id"] for _, rows in plain_items for row in rows]
    gz_ids = [row["id"] for _, rows in gz_items for row in rows]
    assert gz_ids == plain_ids == ["1", "2"]


def test_process_xml_file_pass_namespaced_tag(tmp_path: Path) -> None:
    """Verify namespace-qualified tags (e.g. UniProt-style {ns}entry) are matched and parsed."""
    file_path = _write_xml(tmp_path, "uniprot.xml", UNIPROT_XML_2)
    xml_tag = f"{{{UNIPROT_NS}}}entry"
    items = _table_and_data(process_xml_file(fake_settings(), xml_tag, parse_uniprot_entry, file_path=file_path))  # pyright: ignore[reportArgumentType]
    accessions = [row["accession"] for _, rows in items for row in rows]
    assert accessions == ["P12345", "Q67890"]


def test_process_xml_file_pass_partial_empty_parse_results(tmp_path: Path) -> None:
    """Verify entries for which parse_fn returns an empty dict contribute nothing, and surviving rows are exact."""
    file_path = _write_xml(tmp_path, "missing_email.xml", PEOPLE_XML_2_ONE_MISSING_EMAIL)
    items = _table_and_data(
        process_xml_file(fake_settings(), "person", parse_person_skip_missing_email, file_path=file_path)  # pyright: ignore[reportArgumentType]
    )
    assert items == [("people", [{"id": "1", "email": "a@example.com"}])]


def test_process_xml_file_pass_no_matching_elements(tmp_path: Path) -> None:
    """Verify a well-formed XML file with no elements matching xml_tag yields nothing, without error."""
    file_path = _write_xml(tmp_path, "empty.xml", PEOPLE_XML_EMPTY)
    items = list(process_xml_file(fake_settings(), "person", parse_person_entry, file_path=file_path))  # pyright: ignore[reportArgumentType]
    assert items == []


def test_process_xml_file_pass_file_path_passed_through(tmp_path: Path) -> None:
    """Verify the correct file_path is forwarded to parse_fn for every entry."""
    file_path = _write_xml(tmp_path, "two.xml", PEOPLE_XML_2)
    items = _table_and_data(process_xml_file(fake_settings(), "person", parse_person_entry, file_path=file_path))  # pyright: ignore[reportArgumentType]
    for _, rows in items:
        assert rows[0]["source_file"] == str(file_path)


def test_process_xml_file_pass_logs_reading_info_message_once(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    """Verify the INFO-level 'Reading from <path>' message is logged exactly once with the correct path."""
    caplog.set_level(logging.INFO)
    file_path = _write_xml(tmp_path, "two.xml", PEOPLE_XML_2)
    list(process_xml_file(fake_settings(), "person", parse_person_entry, file_path=file_path))  # pyright: ignore[reportArgumentType]

    info_records = [r for r in caplog.records if r.levelno == logging.INFO]
    assert len(info_records) == 1
    assert info_records[0].msg == "Reading from %s"
    assert info_records[0].args == (str(file_path),)


def test_process_xml_file_pass_no_matching_elements_logs_no_progress(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """Verify no elements matching xml_tag also means no 'Processed' progress log is emitted (final-if is False)."""
    caplog.set_level(logging.DEBUG)
    file_path = _write_xml(tmp_path, "empty.xml", PEOPLE_XML_EMPTY)
    items = list(process_xml_file(fake_settings(), "person", parse_person_entry, file_path=file_path))  # pyright: ignore[reportArgumentType]

    assert items == []
    progress_records = [r for r in caplog.records if r.msg.startswith("Processed")]
    assert progress_records == []


@pytest.mark.parametrize(
    ("log_interval", "expected_interim_counts", "expected_final_count"),
    [
        pytest.param(1, [1, 2, 3, 4, 5], None, id="divides_evenly_every_entry"),
        pytest.param(5, [5], None, id="divides_evenly_equal_to_count"),
        pytest.param(2, [2, 4], 5, id="remainder_leaves_trailing_entries"),
        pytest.param(10, [], 5, id="interval_larger_than_entry_count"),
    ],
)
def test_process_xml_file_pass_log_interval_boundary(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
    log_interval: int,
    expected_interim_counts: list[int],
    expected_final_count: int | None,
) -> None:
    """Verify the exact entry counts logged at each interim checkpoint and in the trailing summary, not just their number.

    ``log_interval`` is supplied via the settings object and constrained to be a
    positive integer by ``CtsSettings``.
    """
    caplog.set_level(logging.DEBUG)
    file_path = _write_xml(tmp_path, "five.xml", PEOPLE_XML_5)
    list(
        process_xml_file(
            fake_settings(log_interval=log_interval),  # pyright: ignore[reportArgumentType]
            "person",
            parse_person_entry,
            file_path=file_path,
        )
    )

    interim_counts = [r.args[0] for r in caplog.records if r.msg == "Processed %d entries"]
    final_records = [r for r in caplog.records if r.msg == "Processed %d entries from %s"]

    assert interim_counts == expected_interim_counts
    if expected_final_count is None:
        assert final_records == []
    else:
        assert len(final_records) == 1
        assert final_records[0].args == (expected_final_count, file_path.name)


def test_process_xml_file_fail_missing_file(tmp_path: Path) -> None:
    """Verify a nonexistent file path raises FileNotFoundError when the generator is consumed."""
    missing_path = tmp_path / "does_not_exist.xml"
    with pytest.raises(FileNotFoundError, match="No such file or directory"):
        list(process_xml_file(fake_settings(), "person", parse_person_entry, file_path=missing_path))  # pyright: ignore[reportArgumentType]


def test_process_xml_file_fail_malformed_xml(tmp_path: Path) -> None:
    """Verify malformed (not well-formed) XML raises lxml's XMLSyntaxError during streaming."""
    file_path = _write_xml(tmp_path, "malformed.xml", PEOPLE_XML_MALFORMED)
    with pytest.raises(XMLSyntaxError, match="Opening and ending tag mismatch"):
        list(process_xml_file(fake_settings(), "person", parse_person_entry, file_path=file_path))  # pyright: ignore[reportArgumentType]


def test_process_xml_file_fail_parse_fn_raises_mid_stream(tmp_path: Path) -> None:
    """Verify an exception raised before the final buffer flush propagates without yielding buffered rows."""
    file_path = _write_xml(tmp_path, "five.xml", PEOPLE_XML_5)
    failing_parse_fn = make_failing_parse_fn(fail_on_call=3)
    generator: Generator[DataItemWithMeta, Any] = process_xml_file(
        fake_settings(),  # pyright: ignore[reportArgumentType]
        "person",
        failing_parse_fn,
        file_path=file_path,  # pyright: ignore[reportArgumentType]
    )

    collected: list[DataItemWithMeta] = []
    with pytest.raises(ValueError, match="boom on call 3"):  # noqa: PT012
        for item in generator:
            collected.append(item)  # noqa: PERF402

    assert _table_and_data(collected) == []


def test_process_xml_file_fail_parse_fn_returns_non_dict(tmp_path: Path) -> None:
    """Verify a parse_fn violating its dict[str, rows] contract raises AttributeError on .items()."""
    file_path = _write_xml(tmp_path, "two.xml", PEOPLE_XML_2)
    with pytest.raises(AttributeError, match="'list' object has no attribute 'items'"):
        list(process_xml_file(fake_settings(), "person", parse_returns_non_dict, file_path=file_path))  # pyright: ignore[reportArgumentType]


"""process_xml_file_to_dict"""


@pytest.mark.parametrize("gzip_compress", [False, True], ids=["plain", "gzip"])
def test_process_xml_file_to_dict_pass_matches_xmltodict_reference(tmp_path: Path, gzip_compress: bool) -> None:
    """Verify process_xml_file_to_dict matches direct xmltodict parsing of each streamed entry."""
    filename = "uniprot.xml.gz" if gzip_compress else "uniprot.xml"
    file_path = _write_xml(tmp_path, filename, UNIPROT_XML_2, gzip_compress=gzip_compress)
    xml_tag = f"{{{UNIPROT_NS}}}entry"

    raw_entries = [
        xmltodict.parse(tostring(entry), **DEFAULT_XMLTODICT_ARGS)["entry"]  # pyright: ignore[reportArgumentType]
        for entry in xml_module.stream_xml_file(file_path, xml_tag)
    ]
    # sanity check that the fixture actually exercises xmlns stripping, so a no-op filter
    # wouldn't slip through unnoticed.
    assert any(key.startswith("_xmlns") for entry in raw_entries for key in entry)
    expected_rows = [{"entry": {k: v for k, v in entry.items() if not k.startswith("_xmlns")}} for entry in raw_entries]

    settings = fake_xmltodict_settings(table_name="some_table", xml_tag=xml_tag)
    tagged = _table_and_data(process_xml_file_to_dict(settings, file_path=file_path))

    assert tagged == [("some_table", expected_rows)]


def test_process_xml_file_to_dict_pass_strips_only_xmlns_prefixed_keys(tmp_path: Path) -> None:
    """Verify the exact, hand-specified result: xmlns/xmlns:xsi attrs removed, other underscore-prefixed keys kept."""
    namespace = "http://example.com/ns"
    xml_with_xmlns_and_lookalike = (
        '<?xml version="1.0"?>\n'
        f'<person xmlns="{namespace}" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
        'id="1" xsinil="true"><name>Anne Example</name></person>\n'
    )
    file_path = _write_xml(tmp_path, "namespaced.xml", xml_with_xmlns_and_lookalike)
    settings = fake_xmltodict_settings(xml_tag=f"{{{namespace}}}person")

    items = _table_and_data(process_xml_file_to_dict(settings, file_path=file_path))

    # "_xmlns" and "_xmlns:xsi" (the real namespace declarations) are stripped; "_id" and the
    # underscore-prefixed-but-unrelated "_xsinil" (from a plain, non-namespace "xsinil" attribute)
    # are both left alone.
    assert items == [
        (
            "people",
            [{"person": {"_id": "1", "_xsinil": "true", "name": "Anne Example"}}],
        )
    ]


def test_process_xml_file_to_dict_pass_custom_table_name_and_xml_tag_are_respected(tmp_path: Path) -> None:
    """Verify rows land under settings.table_name, matched by settings.xml_tag, not any hardcoded name."""
    file_path = _write_xml(tmp_path, "two.xml", PEOPLE_XML_2)
    settings = fake_xmltodict_settings(table_name="library_entries", xml_tag="person")

    tagged = _table_and_data(process_xml_file_to_dict(settings, file_path=file_path))

    assert [table_name for table_name, _ in tagged] == ["library_entries"]
    assert tagged[0][1] == PEOPLE_XML_TO_DICT_PARSED[:2]


@pytest.mark.parametrize(
    ("buffer_size", "expected_batch_sizes"),
    [
        pytest.param(1, [1, 1, 1, 1, 1], id="one-row-batches"),
        pytest.param(2, [2, 2, 1], id="two-row-batches-with-remainder"),
        pytest.param(3, [3, 2], id="three-row-batches-with-remainder"),
        pytest.param(10, [5], id="larger-than-input"),
    ],
)
def test_process_xml_file_to_dict_pass_buffer_size_controls_batching(
    tmp_path: Path, buffer_size: int, expected_batch_sizes: list[int]
) -> None:
    """Verify rows are yielded in batches of at most `buffer_size`, with a final partial batch."""
    file_path = _write_xml(tmp_path, "five.xml", PEOPLE_XML_5)
    settings = fake_xmltodict_settings(buffer_size=buffer_size)

    items = _table_and_data(process_xml_file_to_dict(settings, file_path=file_path))

    assert [len(rows) for _, rows in items] == expected_batch_sizes
    assert [row for _, rows in items for row in rows] == PEOPLE_XML_TO_DICT_PARSED
    assert {table_name for table_name, _ in items} == {"people"}


def test_process_xml_file_to_dict_pass_no_matching_elements_yields_nothing_and_logs_no_progress(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """Verify a well-formed file with no elements matching xml_tag yields nothing and logs no progress."""
    caplog.set_level(logging.DEBUG)
    file_path = _write_xml(tmp_path, "empty.xml", PEOPLE_XML_EMPTY)
    settings = fake_xmltodict_settings()

    items = list(process_xml_file_to_dict(settings, file_path=file_path))

    assert items == []
    assert [r for r in caplog.records if r.msg.startswith("Processed")] == []


def test_process_xml_file_to_dict_pass_logs_reading_info_message_once(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """Verify the INFO-level 'Reading from <path>' message is logged exactly once with the correct path."""
    caplog.set_level(logging.INFO)
    file_path = _write_xml(tmp_path, "five.xml", PEOPLE_XML_5)
    settings = fake_xmltodict_settings()

    list(process_xml_file_to_dict(settings, file_path=file_path))

    info_records = [r for r in caplog.records if r.levelno == logging.INFO]
    assert len(info_records) == 1
    assert info_records[0].msg == "Reading from %s"
    assert info_records[0].args == (str(file_path),)


@pytest.mark.parametrize(
    ("log_interval", "expected_interim_counts", "expected_final_count"),
    [
        pytest.param(1, [1, 2, 3, 4, 5], None, id="divides_evenly_every_entry"),
        pytest.param(5, [5], None, id="divides_evenly_equal_to_count"),
        pytest.param(2, [2, 4], 5, id="remainder_leaves_trailing_entries"),
        pytest.param(10, [], 5, id="interval_larger_than_entry_count"),
    ],
)
def test_process_xml_file_to_dict_pass_log_interval_boundary(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
    log_interval: int,
    expected_interim_counts: list[int],
    expected_final_count: int | None,
) -> None:
    """Verify the exact entry counts logged at each interim checkpoint and in the trailing summary."""
    caplog.set_level(logging.DEBUG)
    file_path = _write_xml(tmp_path, "five.xml", PEOPLE_XML_5)
    settings = fake_xmltodict_settings(log_interval=log_interval)

    list(process_xml_file_to_dict(settings, file_path=file_path))

    interim_counts = [r.args[0] for r in caplog.records if r.msg == "Processed %d entries"]
    final_records = [r for r in caplog.records if r.msg == "Processed %d entries from %s"]

    assert interim_counts == expected_interim_counts
    if expected_final_count is None:
        assert final_records == []
    else:
        assert len(final_records) == 1
        assert final_records[0].args == (expected_final_count, file_path.name)


def test_process_xml_file_to_dict_pass_gzip_file(tmp_path: Path) -> None:
    """Verify gzip-compressed XML input (.gz suffix) is transparently decompressed and parsed."""
    plain_path = _write_xml(tmp_path, "plain.xml", PEOPLE_XML_5)
    gz_path = _write_xml(tmp_path, "compressed.xml.gz", PEOPLE_XML_5, gzip_compress=True)
    settings = fake_xmltodict_settings()

    plain_items = _table_and_data(process_xml_file_to_dict(settings, file_path=plain_path))
    gz_items = _table_and_data(process_xml_file_to_dict(fake_xmltodict_settings(), file_path=gz_path))

    plain_rows = [row for _, rows in plain_items for row in rows]
    gz_rows = [row for _, rows in gz_items for row in rows]
    assert gz_rows == plain_rows == PEOPLE_XML_TO_DICT_PARSED


def test_process_xml_file_to_dict_fail_missing_file(tmp_path: Path) -> None:
    """Verify a nonexistent file path raises FileNotFoundError when the generator is consumed."""
    missing_path = tmp_path / "does_not_exist.xml"
    settings = fake_xmltodict_settings()
    with pytest.raises(FileNotFoundError, match="No such file or directory"):
        list(process_xml_file_to_dict(settings, file_path=missing_path))


def test_process_xml_file_to_dict_fail_malformed_xml(tmp_path: Path) -> None:
    """Verify malformed (not well-formed) XML raises lxml's XMLSyntaxError during streaming."""
    file_path = _write_xml(tmp_path, "malformed.xml", PEOPLE_XML_MALFORMED)
    settings = fake_xmltodict_settings()
    with pytest.raises(XMLSyntaxError, match="Opening and ending tag mismatch"):
        list(process_xml_file_to_dict(settings, file_path=file_path))


# The <entry> subtree below matches tests/data/xsd/uniref_like.xsd: "property" may repeat under
# both "entry" and "dbReference" (list_paths), while "representativeMember" and "dbReference" may
# each occur at most once (single_paths) -- see test_xsd.py's
# test_find_list_and_single_child_paths_pass_uniref_like_schema for the fully-verified path sets.
# Every list-typed child below occurs exactly once in the data, so only schema knowledge (not
# occurrence count) can tell process_xml_file_to_dict it must still become a list.
UNIREF_LIKE_ENTRY_XML: Final[str] = """<?xml version="1.0"?>
<UniRef50>
    <entry>
        <property type="member count" value="1"/>
        <representativeMember>
            <dbReference>
                <property type="isSeed" value="true"/>
            </dbReference>
        </representativeMember>
    </entry>
</UniRef50>
"""


def _build_real_xmltodict_settings(tmp_path: Path, test_data_dir: Path, xsd_filename: str | None) -> XmlToDictSettings:
    """Build a real (non-mock) XmlToDictSettings, optionally with an xsd_file copied from tests/data/xsd.

    A real settings object is required here, rather than fake_xmltodict_settings, because
    xmltodict_args is a computed property with real schema-derived logic to exercise: a mock
    would either need that same logic re-implemented by hand (proving nothing about the actual
    implementation) or would unpack to `{}` by default, silently skipping the feature entirely.
    """
    input_dir = tmp_path / "input"
    input_dir.mkdir(exist_ok=True)
    output_dir = tmp_path / "output"
    output_dir.mkdir(exist_ok=True)
    log_config_file = tmp_path / "logging.json"
    log_config_file.write_text('{"version": 1}', encoding="utf-8")

    xsd_file = None
    if xsd_filename is not None:
        xsd_content = (test_data_dir / "xsd" / xsd_filename).read_text(encoding="utf-8")
        xsd_path = input_dir / xsd_filename
        xsd_path.write_text(xsd_content, encoding="utf-8")
        xsd_file = str(xsd_path)

    return XmlToDictSettings(
        input_dir=str(input_dir),
        output_dir=str(output_dir),
        use_destination="local_fs",
        log_config_file=str(log_config_file),
        dataset_name="xsd_test_dataset",
        table_name="entry",
        xml_tag="entry",
        xsd_file=xsd_file,
        buffer_size=100,
        log_interval=1000,
    )  # pyright: ignore[reportCallIssue]


@pytest.mark.parametrize(
    ("xsd_filename", "expected_entry"),
    [
        pytest.param(
            "uniref_like.xsd",
            {
                "property": [{"_type": "member count", "_value": "1"}],
                "representativeMember": {
                    "dbReference": {"property": [{"_type": "isSeed", "_value": "true"}]},
                },
            },
            id="with-xsd-file-forces-repeatable-schema-children-into-lists",
        ),
        pytest.param(
            None,
            {
                "property": {"_type": "member count", "_value": "1"},
                "representativeMember": {
                    "dbReference": {"property": {"_type": "isSeed", "_value": "true"}},
                },
            },
            id="without-xsd-file-keeps-xmltodicts-default-single-occurrence-dicts",
        ),
    ],
)
def test_process_xml_file_to_dict_pass_force_list_shapes_children_by_schema_occurs(
    tmp_path: Path,
    test_data_dir: Path,
    xsd_filename: str | None,
    expected_entry: dict[str, Any],
) -> None:
    """Verify children the schema allows to repeat become lists even with a single occurrence.

    With xsd_file set, every child the schema allows to repeat must be a list, regardless of how
    many times it actually occurs in this particular document, while children the schema allows
    only once stay plain dicts -- at every nesting depth. Without xsd_file, xmltodict's ordinary
    behaviour applies instead: a lone occurrence is never turned into a list, whether or not the
    schema would allow more.
    """
    file_path = _write_xml(tmp_path, "entry.xml", UNIREF_LIKE_ENTRY_XML)
    settings = _build_real_xmltodict_settings(tmp_path, test_data_dir, xsd_filename)

    items = _table_and_data(process_xml_file_to_dict(settings, file_path=file_path))

    assert items == [("entry", [{"entry": expected_entry}])]


"""process_xml_file_items"""


def _file_items(directory: Path) -> list[FileItemDict]:
    """Collect the file items emitted by the dlt filesystem source for XML files in a directory, in name order."""
    items = list(filesystem(bucket_url=str(directory), file_glob=DEFAULT_XML_FILE_GLOB))
    assert all(isinstance(item, FileItemDict) for item in items)
    return sorted(items, key=lambda item: item["file_name"])


def test_process_xml_file_items_pass_multiple_files_preserve_order(tmp_path: Path) -> None:
    """Verify files are processed in the order they are supplied, including gzipped files."""
    file_a = _write_xml(tmp_path, "a_00001.xml", PEOPLE_XML_2)
    file_b = _write_xml(tmp_path, "b_00002.xml.gz", PEOPLE_XML_5, gzip_compress=True)

    items = _table_and_data(
        process_xml_file_items(
            _file_items(tmp_path),
            fake_settings(input_dir=tmp_path),  # pyright: ignore[reportArgumentType]
            "person",
            parse_person_entry,
        )
    )

    sources = [row["source_file"] for _, rows in items for row in rows]
    assert sources == [str(file_a)] * 2 + [str(file_b)] * 5


def test_process_xml_file_items_pass_file_path_is_input_dir_joined_to_relative_path(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Verify the path given to the parser is input_dir plus the file's path relative to the bucket."""
    sub_dir = tmp_path / "sub"
    sub_dir.mkdir()
    _write_xml(sub_dir, "nested_00001.xml", PEOPLE_XML_2)
    items = list(filesystem(bucket_url=str(tmp_path), file_glob="**/*.xml"))
    process_file_mock = MagicMock(side_effect=lambda *_, **__: iter([]))
    monkeypatch.setattr(xml_module, "process_xml_file", process_file_mock)

    list(process_xml_file_items(items, fake_settings(input_dir="relative_dir"), "person", parse_person_entry))  # pyright: ignore[reportArgumentType]

    process_file_mock.assert_called_once()
    assert process_file_mock.call_args.kwargs["file_path"] == Path("relative_dir") / "sub" / "nested_00001.xml"


def test_process_xml_file_items_pass_no_items_yields_nothing() -> None:
    """Verify an empty iterable of file items produces no output."""
    assert list(process_xml_file_items([], fake_settings(), "person", parse_person_entry)) == []  # pyright: ignore[reportArgumentType]


@pytest.mark.parametrize("buffer_size", [1, 2, 50])
def test_process_xml_file_items_pass_settings_forwarded(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, buffer_size: int
) -> None:
    """Verify the settings object (carrying buffer_size/log_interval) and file path reach process_xml_file."""
    file_path = _write_xml(tmp_path, "two_00001.xml", PEOPLE_XML_2)
    settings = fake_settings(buffer_size=buffer_size, log_interval=CUSTOM_LOG_INTERVAL, input_dir=tmp_path)
    spy = MagicMock(wraps=xml_module.process_xml_file)
    monkeypatch.setattr(xml_module, "process_xml_file", spy)

    list(process_xml_file_items(_file_items(tmp_path), settings, "person", parse_person_entry))  # pyright: ignore[reportArgumentType]

    spy.assert_called_once_with(
        settings=settings, xml_tag="person", parse_fn=parse_person_entry, file_path=Path(str(file_path))
    )
    assert spy.call_args.kwargs["settings"].buffer_size == buffer_size
    assert spy.call_args.kwargs["settings"].log_interval == CUSTOM_LOG_INTERVAL


def test_process_xml_file_items_fail_missing_file_mid_iteration(tmp_path: Path) -> None:
    """Verify a missing file raises FileNotFoundError after earlier files were fully processed."""
    good_path = _write_xml(tmp_path, "a_00001.xml", PEOPLE_XML_2)
    items = _file_items(tmp_path)
    missing_item = FileItemDict({**items[0], "relative_path": "missing_00002.xml"})

    collected: list[DataItemWithMeta] = []
    with pytest.raises(FileNotFoundError, match="No such file or directory"):  # noqa: PT012
        for item in process_xml_file_items(
            [*items, missing_item],
            fake_settings(input_dir=tmp_path),  # pyright: ignore[reportArgumentType]
            "person",
            parse_person_entry,
        ):
            collected.append(item)  # noqa: PERF402

    tagged = _table_and_data(collected)
    assert [table for table, _ in tagged] == ["people"]
    assert [row["id"] for row in tagged[0][1]] == ["1", "2"]
    assert {row["source_file"] for row in tagged[0][1]} == {str(good_path)}


def test_process_xml_file_pass_calls_stream_xml_file_with_correct_args(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Verify process_xml_file invokes stream_xml_file exactly once with (file_path, xml_tag)."""
    stream_mock = MagicMock(return_value=iter([]))
    monkeypatch.setattr(xml_module, "stream_xml_file", stream_mock)
    file_path = Path("dummy.xml")

    list(process_xml_file(fake_settings(), "tag", MagicMock(), file_path=file_path))  # pyright: ignore[reportArgumentType]

    stream_mock.assert_called_once_with(file_path, "tag")


def test_process_xml_file_pass_calls_parse_fn_with_expected_keyword_args(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Verify parse_fn is invoked once per streamed entry with entry/file_path as kwargs."""
    fake_entry = MagicMock(spec=Element)
    monkeypatch.setattr(xml_module, "stream_xml_file", lambda _, __: iter([fake_entry]))
    parse_fn = MagicMock(return_value={})
    file_path = Path("dummy.xml")

    list(process_xml_file(fake_settings(), "tag", parse_fn, file_path=file_path))  # pyright: ignore[reportArgumentType]

    parse_fn.assert_called_once_with(entry=fake_entry, file_path=file_path)


def test_process_xml_file_pass_wraps_each_table_with_dlt_mark(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Verify each (table, rows) pair from parse_fn's return dict is forwarded to dlt.mark.with_table_name."""
    fake_entry = MagicMock(spec=Element)
    monkeypatch.setattr(xml_module, "stream_xml_file", lambda *_: iter([fake_entry]))
    rows_people = [{"id": 1}]
    rows_emails = [{"email": "a@example.com"}]
    parse_fn = MagicMock(return_value={"people": rows_people, "emails": rows_emails})
    mark_mock = MagicMock(side_effect=lambda rows, table: DataItemWithMeta(TableNameMeta(table), rows))
    monkeypatch.setattr(buffer_module.dlt.mark, "with_table_name", mark_mock)

    items = list(process_xml_file(fake_settings(), "tag", parse_fn, file_path=Path("f.xml")))  # pyright: ignore[reportArgumentType]

    mark_mock.assert_any_call(rows_people, "people")
    mark_mock.assert_any_call(rows_emails, "emails")
    assert len(items) == 2
    assert items[0].meta.table_name == "people"
    assert items[0].data == rows_people
    assert items[1].meta.table_name == "emails"
    assert items[1].data == rows_emails


def test_process_xml_file_pass_skips_parse_fn_when_no_entries(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Verify parse_fn is never invoked, and no items are yielded, when stream_xml_file yields no elements."""
    monkeypatch.setattr(xml_module, "stream_xml_file", lambda _fp, _tag: iter([]))
    parse_fn = MagicMock()

    items = list(process_xml_file(fake_settings(), "tag", parse_fn, file_path=Path("f.xml")))  # pyright: ignore[reportArgumentType]

    parse_fn.assert_not_called()
    assert items == []


def test_process_xml_file_fail_propagates_stream_xml_file_error(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Verify the first entry is parsed correctly before an exception raised mid-iteration by stream_xml_file propagates."""
    fake_entry = MagicMock(spec=Element)

    def _raising_stream(_fp: Path, _tag: str) -> Iterator[Element]:
        yield fake_entry
        msg = "disk read error"
        raise OSError(msg)

    monkeypatch.setattr(xml_module, "stream_xml_file", _raising_stream)
    parse_fn = MagicMock(return_value={})
    file_path = Path("f.xml")

    collected: list[DataItemWithMeta] = []
    with pytest.raises(OSError, match="disk read error"):  # noqa: PT012
        for item in process_xml_file(fake_settings(), "tag", parse_fn, file_path=file_path):  # pyright: ignore[reportArgumentType]
            collected.append(item)  # noqa: PERF402

    parse_fn.assert_called_once_with(entry=fake_entry, file_path=file_path)
    assert collected == []


# process_xml_file_batches - isolated unit tests (get_file_batches / process_xml_file mocked)
def test_process_xml_file_batches_pass_calls_get_file_batches_once_with_settings(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Verify get_file_batches is invoked exactly once, receiving the settings object unchanged."""
    settings = fake_settings()
    get_batches_mock = MagicMock(return_value=iter([]))
    monkeypatch.setattr(xml_module, "get_file_batches", get_batches_mock)

    items = list(process_xml_file_batches(settings, "tag", MagicMock()))  # pyright: ignore[reportArgumentType]

    get_batches_mock.assert_called_once_with(settings)
    assert items == []


def test_process_xml_file_batches_pass_calls_process_xml_file_per_file_in_order(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Verify process_xml_file is called once per file across all batches, in order, with forwarded args."""
    file_a, file_b, file_c = Path("a.xml"), Path("b.xml"), Path("c.xml")
    settings = fake_settings(log_interval=CUSTOM_BATCH_LOG_INTERVAL)
    monkeypatch.setattr(xml_module, "get_file_batches", lambda _: iter([[file_a, file_b], [file_c]]))
    process_file_mock = MagicMock(side_effect=lambda *_, **__: iter([]))
    monkeypatch.setattr(xml_module, "process_xml_file", process_file_mock)
    parse_fn = MagicMock()

    list(process_xml_file_batches(settings, "tag", parse_fn))  # pyright: ignore[reportArgumentType]

    assert [call.kwargs.get("file_path") for call in process_file_mock.call_args_list] == [
        file_a,
        file_b,
        file_c,
    ]
    for call in process_file_mock.call_args_list:
        assert call.kwargs.get("xml_tag") == "tag"
        assert call.kwargs.get("parse_fn") is parse_fn
        assert call.kwargs.get("settings") is settings
        assert call.kwargs.get("settings").log_interval == CUSTOM_BATCH_LOG_INTERVAL  # pyright: ignore[reportOptionalMemberAccess]


def test_process_xml_file_batches_pass_yields_items_from_process_xml_file_unchanged(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Verify items yielded by process_xml_file are passed through by process_xml_file_batches as-is."""
    sentinel_item = DataItemWithMeta(TableNameMeta("people"), [{"id": 1}])
    monkeypatch.setattr(xml_module, "get_file_batches", lambda _: iter([[Path("a.xml")]]))
    monkeypatch.setattr(xml_module, "process_xml_file", lambda *_, **__: iter([sentinel_item]))

    items = list(process_xml_file_batches(fake_settings(), "tag", MagicMock()))  # pyright: ignore[reportArgumentType]

    assert items == [sentinel_item]


def test_process_xml_file_batches_fail_propagates_process_xml_file_error(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Verify an exception raised inside process_xml_file for one file propagates without being swallowed."""

    def _raising_process(*_: Any, **__: Any) -> Iterator[DataItemWithMeta]:  # noqa: ANN401
        msg = "parse failure"
        raise RuntimeError(msg)
        yield  # pragma: no cover - unreachable, keeps this a generator function

    monkeypatch.setattr(xml_module, "get_file_batches", lambda _: iter([[Path("a.xml")]]))
    monkeypatch.setattr(xml_module, "process_xml_file", _raising_process)

    with pytest.raises(RuntimeError, match="parse failure"):
        list(process_xml_file_batches(fake_settings(), "tag", MagicMock()))  # pyright: ignore[reportArgumentType]
