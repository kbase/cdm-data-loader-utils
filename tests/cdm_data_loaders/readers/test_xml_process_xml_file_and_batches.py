"""Tests for `process_xml_file`, `process_xml_file_batches`, and `xml_to_dict_parse_fn` from cdm_data_loaders.readers.xml."""

import gzip
import logging
from collections.abc import Callable, Generator, Iterable, Iterator
from pathlib import Path
from typing import Any, Final
from unittest.mock import MagicMock

import pytest
import xmltodict
from dlt.extract.items import DataItemWithMeta, TableNameMeta
from lxml.etree import Element, XMLSyntaxError, tostring

import cdm_data_loaders.readers.xml as xml_module
from cdm_data_loaders.core.fields import DEFAULTS, BUFFER_SIZE, LOG_INTERVAL
from cdm_data_loaders.core.settings import BatchedFileInputSettings
from cdm_data_loaders.readers.xml import (
    process_xml_file,
    process_xml_file_batches,
    xml_to_dict_parse_fn,
)

ParseFn = Callable[..., dict[str, list[dict[str, Any]]]]

UNIPROT_NS: Final[str] = "http://uniprot.org/uniprot"

# Arbitrary non-default values used to prove settings are forwarded to process_xml_file.
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
) -> BatchedFileInputSettings:
    """Build a MagicMock stand-in for BatchedFileInputSettings.

    ``process_xml_file`` reads ``buffer_size`` and ``log_interval`` off the settings
    object, so these are exposed as real integers (rather than nested mocks) to make
    the arithmetic in the function under test behave correctly.
    """
    settings = MagicMock(spec=BatchedFileInputSettings)
    settings.buffer_size = buffer_size
    settings.log_interval = log_interval
    return settings


def test_process_xml_file_pass_single_table_single_entry(tmp_path: Path) -> None:
    """Verify a single matching element is parsed into one correctly table-tagged row."""
    file_path = _write_xml(tmp_path, "one.xml", PEOPLE_XML_2)
    items = list(process_xml_file(fake_settings(), "person", parse_person_entry, file_path=file_path))
    assert len(items) == 1
    parsed = _table_and_data(items)
    assert parsed[0][0] == "people"
    assert parsed[0][1] == [{"source_file": str(file_path), **p_data} for p_data in PEOPLE_PARSED[:2]]


def test_process_xml_file_pass_multiple_entries(tmp_path: Path) -> None:
    """Verify all entries in a multi-entry XML file are streamed and parsed in document order."""
    file_path = _write_xml(tmp_path, "five.xml", PEOPLE_XML_5)
    items = list(process_xml_file(fake_settings(), "person", parse_person_entry, file_path=file_path))
    parsed = _table_and_data(items)
    assert parsed[0][0] == "people"
    assert parsed[0][1] == [{"source_file": str(file_path), **p_data} for p_data in PEOPLE_PARSED]


@pytest.mark.parametrize("gzip_compress", [False, True], ids=["plain", "gzip"])
def test_process_xml_file_pass_xml_to_dict_parse_fn_implements_xmltodict(tmp_path: Path, gzip_compress: bool) -> None:
    """Verify xml_to_dict_parse_fn matches direct xmltodict parsing of each streamed entry."""
    filename = "people.xml.gz" if gzip_compress else "people.xml"
    file_path = _write_xml(tmp_path, filename, PEOPLE_XML_2, gzip_compress=gzip_compress)

    expected_rows = [xmltodict.parse(tostring(entry)) for entry in xml_module.stream_xml_file(file_path, "person")]
    tagged = _table_and_data(
        process_xml_file(fake_settings(), "person", xml_to_dict_parse_fn("xmltodict"), file_path=file_path)
    )

    assert tagged == [("xmltodict", expected_rows)]


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
    positive integer by ``BatchedFileInputSettings``.
    """
    file_path = _write_xml(tmp_path, "five.xml", PEOPLE_XML_5)

    items = _table_and_data(
        process_xml_file(
            fake_settings(buffer_size=buffer_size),
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
    items = list(process_xml_file(fake_settings(), "person", parse_person_multi_table, file_path=file_path))
    tagged = _table_and_data(items)
    assert [table for table, _ in tagged] == ["people", "emails"]
    assert [row["id"] for row in tagged[0][1]] == ["1", "2"]
    assert [row["email"] for row in tagged[1][1]] == ["a@example.com", "b@example.com"]


def test_process_xml_file_pass_gzip_file(tmp_path: Path) -> None:
    """Verify gzip-compressed XML input (.gz suffix) is transparently decompressed and parsed."""
    plain_path = _write_xml(tmp_path, "plain.xml", PEOPLE_XML_2)
    gz_path = _write_xml(tmp_path, "compressed.xml.gz", PEOPLE_XML_2, gzip_compress=True)

    plain_items = _table_and_data(process_xml_file(fake_settings(), "person", parse_person_entry, file_path=plain_path))
    gz_items = _table_and_data(process_xml_file(fake_settings(), "person", parse_person_entry, file_path=gz_path))

    plain_ids = [row["id"] for _, rows in plain_items for row in rows]
    gz_ids = [row["id"] for _, rows in gz_items for row in rows]
    assert gz_ids == plain_ids == ["1", "2"]


def test_process_xml_file_pass_namespaced_tag(tmp_path: Path) -> None:
    """Verify namespace-qualified tags (e.g. UniProt-style {ns}entry) are matched and parsed."""
    file_path = _write_xml(tmp_path, "uniprot.xml", UNIPROT_XML_2)
    xml_tag = f"{{{UNIPROT_NS}}}entry"
    items = _table_and_data(process_xml_file(fake_settings(), xml_tag, parse_uniprot_entry, file_path=file_path))
    accessions = [row["accession"] for _, rows in items for row in rows]
    assert accessions == ["P12345", "Q67890"]


def test_process_xml_file_pass_partial_empty_parse_results(tmp_path: Path) -> None:
    """Verify entries for which parse_fn returns an empty dict contribute nothing, and surviving rows are exact."""
    file_path = _write_xml(tmp_path, "missing_email.xml", PEOPLE_XML_2_ONE_MISSING_EMAIL)
    items = _table_and_data(
        process_xml_file(fake_settings(), "person", parse_person_skip_missing_email, file_path=file_path)
    )
    assert items == [("people", [{"id": "1", "email": "a@example.com"}])]


def test_process_xml_file_pass_no_matching_elements(tmp_path: Path) -> None:
    """Verify a well-formed XML file with no elements matching xml_tag yields nothing, without error."""
    file_path = _write_xml(tmp_path, "empty.xml", PEOPLE_XML_EMPTY)
    items = list(process_xml_file(fake_settings(), "person", parse_person_entry, file_path=file_path))
    assert items == []


def test_process_xml_file_pass_file_path_passed_through(tmp_path: Path) -> None:
    """Verify the correct file_path is forwarded to parse_fn for every entry."""
    file_path = _write_xml(tmp_path, "two.xml", PEOPLE_XML_2)
    items = _table_and_data(process_xml_file(fake_settings(), "person", parse_person_entry, file_path=file_path))
    for _, rows in items:
        assert rows[0]["source_file"] == str(file_path)


def test_process_xml_file_pass_logs_reading_info_message_once(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    """Verify the INFO-level 'Reading from <path>' message is logged exactly once with the correct path."""
    caplog.set_level(logging.INFO)
    file_path = _write_xml(tmp_path, "two.xml", PEOPLE_XML_2)
    list(process_xml_file(fake_settings(), "person", parse_person_entry, file_path=file_path))

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
    items = list(process_xml_file(fake_settings(), "person", parse_person_entry, file_path=file_path))

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
    positive integer by ``BatchedFileInputSettings``.
    """
    caplog.set_level(logging.DEBUG)
    file_path = _write_xml(tmp_path, "five.xml", PEOPLE_XML_5)
    list(
        process_xml_file(
            fake_settings(log_interval=log_interval),
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
        list(process_xml_file(fake_settings(), "person", parse_person_entry, file_path=missing_path))


def test_process_xml_file_fail_malformed_xml(tmp_path: Path) -> None:
    """Verify malformed (not well-formed) XML raises lxml's XMLSyntaxError during streaming."""
    file_path = _write_xml(tmp_path, "malformed.xml", PEOPLE_XML_MALFORMED)
    with pytest.raises(XMLSyntaxError, match="Opening and ending tag mismatch"):
        list(process_xml_file(fake_settings(), "person", parse_person_entry, file_path=file_path))


def test_process_xml_file_fail_parse_fn_raises_mid_stream(tmp_path: Path) -> None:
    """Verify an exception raised before the final buffer flush propagates without yielding buffered rows."""
    file_path = _write_xml(tmp_path, "five.xml", PEOPLE_XML_5)
    failing_parse_fn = make_failing_parse_fn(fail_on_call=3)
    generator: Generator[DataItemWithMeta, Any] = process_xml_file(
        fake_settings(), "person", failing_parse_fn, file_path=file_path
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
        list(process_xml_file(fake_settings(), "person", parse_returns_non_dict, file_path=file_path))


"""process_xml_file_batches"""


def test_process_xml_file_batches_pass_single_batch_single_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Verify a single batch containing one file is processed identically to calling process_xml_file directly."""
    file_path = _write_xml(tmp_path, "two.xml", PEOPLE_XML_2)
    monkeypatch.setattr(xml_module, "get_file_batches", lambda _: iter([[file_path]]))

    items = _table_and_data(process_xml_file_batches(fake_settings(), "person", parse_person_entry))
    ids = [row["id"] for _, rows in items for row in rows]
    assert ids == ["1", "2"]


def test_process_xml_file_batches_pass_multiple_batches_preserve_order(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Verify files across multiple batches are processed in batch order, then file order, within each."""
    file_a = _write_xml(tmp_path, "a.xml", PEOPLE_XML_2)
    file_b = _write_xml(tmp_path, "b.xml", PEOPLE_XML_5)
    monkeypatch.setattr(xml_module, "get_file_batches", lambda _: iter([[file_a], [file_b]]))

    items = _table_and_data(process_xml_file_batches(fake_settings(), "person", parse_person_entry))
    sources = [row["source_file"] for _, rows in items for row in rows]
    assert sources == [str(file_a)] * 2 + [str(file_b)] * 5


def test_process_xml_file_batches_pass_empty_batches(monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify a get_file_batches iterable with no batches produces no output items."""
    monkeypatch.setattr(xml_module, "get_file_batches", lambda _: iter([]))
    items = list(process_xml_file_batches(fake_settings(), "person", parse_person_entry))
    assert items == []


def test_process_xml_file_batches_pass_batch_with_empty_file_list(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Verify an empty file list within a batch is skipped gracefully and later batches still process."""
    file_path = _write_xml(tmp_path, "two.xml", PEOPLE_XML_2)
    monkeypatch.setattr(xml_module, "get_file_batches", lambda _: iter([[], [file_path]]))

    items = _table_and_data(process_xml_file_batches(fake_settings(), "person", parse_person_entry))
    ids = [row["id"] for _, rows in items for row in rows]
    assert ids == ["1", "2"]


@pytest.mark.parametrize("buffer_size", [1, 2, 50])
def test_process_xml_file_batches_pass_settings_forwarded(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, buffer_size: int
) -> None:
    """Verify the settings object (carrying buffer_size/log_interval) is forwarded to process_xml_file."""
    file_path = _write_xml(tmp_path, "two.xml", PEOPLE_XML_2)
    settings = fake_settings(buffer_size=buffer_size, log_interval=CUSTOM_LOG_INTERVAL)
    monkeypatch.setattr(xml_module, "get_file_batches", lambda _: iter([[file_path]]))
    spy = MagicMock(wraps=xml_module.process_xml_file)
    monkeypatch.setattr(xml_module, "process_xml_file", spy)

    list(process_xml_file_batches(settings, "person", parse_person_entry))

    spy.assert_called_once()
    forwarded_settings = spy.call_args.kwargs.get("settings")
    assert forwarded_settings is not None
    assert forwarded_settings is settings
    assert forwarded_settings.log_interval == CUSTOM_LOG_INTERVAL
    assert forwarded_settings.buffer_size == buffer_size


def test_process_xml_file_batches_pass_default_settings_forwarded(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Verify the default buffer_size/log_interval on the settings object reach process_xml_file."""
    file_path = _write_xml(tmp_path, "two.xml", PEOPLE_XML_2)
    settings = fake_settings()
    monkeypatch.setattr(xml_module, "get_file_batches", lambda _: iter([[file_path]]))
    spy = MagicMock(wraps=xml_module.process_xml_file)
    monkeypatch.setattr(xml_module, "process_xml_file", spy)

    list(process_xml_file_batches(settings, "person", parse_person_entry))

    forwarded_settings = spy.call_args.kwargs.get("settings")
    assert forwarded_settings is not None
    assert forwarded_settings is settings
    assert forwarded_settings.log_interval == DEFAULTS[LOG_INTERVAL]
    assert forwarded_settings.buffer_size == DEFAULTS[BUFFER_SIZE]


def test_process_xml_file_batches_pass_process_xml_file_not_called_for_empty_batches(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Verify process_xml_file is never invoked for empty file-list batches, and invoked exactly once for the real file."""
    file_path = Path("only_file.xml")
    monkeypatch.setattr(xml_module, "get_file_batches", lambda _: iter([[], [file_path], []]))
    process_file_mock = MagicMock(side_effect=lambda *_, **__: iter([]))
    monkeypatch.setattr(xml_module, "process_xml_file", process_file_mock)

    list(process_xml_file_batches(fake_settings(), "tag", MagicMock()))

    process_file_mock.assert_called_once()
    assert process_file_mock.call_args.kwargs.get("file_path") == file_path


def test_process_xml_file_batches_fail_missing_file_mid_batch(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify a missing file later in a batch raises FileNotFoundError after earlier files were fully processed."""
    good_path = _write_xml(tmp_path, "good.xml", PEOPLE_XML_2)
    missing_path = tmp_path / "missing.xml"
    monkeypatch.setattr(xml_module, "get_file_batches", lambda _: iter([[good_path, missing_path]]))

    collected: list[DataItemWithMeta] = []
    with pytest.raises(FileNotFoundError, match="No such file or directory"):  # noqa: PT012
        for item in process_xml_file_batches(fake_settings(), "person", parse_person_entry):
            collected.append(item)  # noqa: PERF402

    tagged = _table_and_data(collected)
    assert [table for table, _ in tagged] == ["people"]
    assert tagged[0][1] == [
        {
            "source_file": str(good_path),
            "id": "1",
            "name": "Anne Example",
            "email": "a@example.com",
        },
        {
            "source_file": str(good_path),
            "id": "2",
            "name": "Belinda Carlisle",
            "email": "b@example.com",
        },
    ]


def test_process_xml_file_batches_fail_get_file_batches_raises(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Verify an exception raised by get_file_batches itself propagates out of the generator."""

    def _raising_get_file_batches(_settings: BatchedFileInputSettings) -> Iterator[list[Path]]:
        msg = "invalid settings"
        raise RuntimeError(msg)

    monkeypatch.setattr(xml_module, "get_file_batches", _raising_get_file_batches)
    with pytest.raises(RuntimeError, match="invalid settings"):
        list(process_xml_file_batches(fake_settings(), "person", parse_person_entry))


def test_process_xml_file_pass_calls_stream_xml_file_with_correct_args(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Verify process_xml_file invokes stream_xml_file exactly once with (file_path, xml_tag)."""
    stream_mock = MagicMock(return_value=iter([]))
    monkeypatch.setattr(xml_module, "stream_xml_file", stream_mock)
    file_path = Path("dummy.xml")

    list(process_xml_file(fake_settings(), "tag", MagicMock(), file_path=file_path))

    stream_mock.assert_called_once_with(file_path, "tag")


def test_process_xml_file_pass_calls_parse_fn_with_expected_keyword_args(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Verify parse_fn is invoked once per streamed entry with entry/file_path as kwargs."""
    fake_entry = MagicMock(spec=Element)
    monkeypatch.setattr(xml_module, "stream_xml_file", lambda _, __: iter([fake_entry]))
    parse_fn = MagicMock(return_value={})
    file_path = Path("dummy.xml")

    list(process_xml_file(fake_settings(), "tag", parse_fn, file_path=file_path))

    parse_fn.assert_called_once_with(entry=fake_entry, file_path=file_path)


def test_process_xml_file_pass_wraps_each_table_with_dlt_mark(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Verify each (table, rows) pair from parse_fn's return dict is forwarded to dlt.mark.with_table_name."""
    fake_entry = MagicMock(spec=Element)
    monkeypatch.setattr(xml_module, "stream_xml_file", lambda fp, tag: iter([fake_entry]))
    rows_people = [{"id": 1}]
    rows_emails = [{"email": "a@example.com"}]
    parse_fn = MagicMock(return_value={"people": rows_people, "emails": rows_emails})
    mark_mock = MagicMock(side_effect=lambda rows, table: DataItemWithMeta(TableNameMeta(table), rows))
    monkeypatch.setattr(xml_module.dlt.mark, "with_table_name", mark_mock)

    items = list(process_xml_file(fake_settings(), "tag", parse_fn, file_path=Path("f.xml")))

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

    items = list(process_xml_file(fake_settings(), "tag", parse_fn, file_path=Path("f.xml")))

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
        for item in process_xml_file(fake_settings(), "tag", parse_fn, file_path=file_path):
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

    items = list(process_xml_file_batches(settings, "tag", MagicMock()))

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

    list(process_xml_file_batches(settings, "tag", parse_fn))

    assert [call.kwargs.get("file_path") for call in process_file_mock.call_args_list] == [
        file_a,
        file_b,
        file_c,
    ]
    for call in process_file_mock.call_args_list:
        assert call.kwargs.get("xml_tag") == "tag"
        assert call.kwargs.get("parse_fn") is parse_fn
        assert call.kwargs.get("settings") is settings
        assert call.kwargs.get("settings").log_interval == CUSTOM_BATCH_LOG_INTERVAL


def test_process_xml_file_batches_pass_yields_items_from_process_xml_file_unchanged(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Verify items yielded by process_xml_file are passed through by process_xml_file_batches as-is."""
    sentinel_item = DataItemWithMeta(TableNameMeta("people"), [{"id": 1}])
    monkeypatch.setattr(xml_module, "get_file_batches", lambda _: iter([[Path("a.xml")]]))
    monkeypatch.setattr(xml_module, "process_xml_file", lambda *_, **__: iter([sentinel_item]))

    items = list(process_xml_file_batches(fake_settings(), "tag", MagicMock()))

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
        list(process_xml_file_batches(fake_settings(), "tag", MagicMock()))
