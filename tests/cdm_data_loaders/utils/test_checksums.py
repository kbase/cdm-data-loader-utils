"""Tests for utils.checksums module — hashlib and CRC64/NVME checksum utilities."""

import base64
import copy
import dataclasses
import hashlib
from collections.abc import Callable, Generator
from io import BytesIO
from pathlib import Path
from typing import Any, Final

import pytest
import requests
from _pytest.mark.structures import ParameterSet
from awscrt.checksums import crc64nvme as _crc64nvme

from cdm_data_loaders.utils.checksums import (
    EXTRA_ALGORITHMS,
    ChecksumEntry,
    Crc64NvmeHasher,
    Hasher,
    HashingReader,
    compute_file_checksum,
    new_hasher,
    validate_checksum_fn,
)
from tests.conftest import DEFAULT_VCR_CONFIG

TEST_BUCKET: Final[str] = "test-bucket"
TEST_KEY: Final[str] = "path/to/file.pdf"
TEST_ECHO_URL: Final[str] = "https://httpbin.org/bytes/1024"  # deterministic-size real
EXPECTED_CRC64_BYTE_LEN = 8


TEST_DATASETS: list[ParameterSet] = [
    pytest.param(b"", id="empty"),
    pytest.param(b"Hello, world!", id="small"),
    pytest.param(b"x" * (2 * (1 << 20) + 123), id="multi-chunk"),
]
TEXT_TRANSFORMATIONS: list[ParameterSet] = [
    pytest.param(lambda x: x, id="identity"),
    pytest.param(lambda x: x.lower(), id="lower"),
    pytest.param(lambda x: x.upper(), id="upper"),
    pytest.param(lambda x: x.capitalize(), id="capitalize"),
    pytest.param(lambda x: x.title().swapcase(), id="togglecase"),
]

HASHLIB_ALGORITHMS: list[str] = [
    "sha1",
    "sha224",
    "sha256",
    "sha384",
    "sha512",
    "sha3_224",
    "sha3_256",
    "sha3_384",
    "sha3_512",
    "blake2b",
    "blake2s",
    "md5",
]
CRC64NVME_ALGORITHMS: list[str] = ["crc64nvme", "CRC64NVME", "Crc64Nvme"]
UNSUPPORTED_ALGORITHMS: list[str] = [
    "",
    "banana",
    "shake_128",
    "shake_256",
    "shake_n_vac",
]


@pytest.fixture
def requests_session() -> Generator[requests.Session, Any]:
    """A requests session."""
    with requests.Session() as s:
        yield s


@pytest.fixture(scope="module")
def vcr_config() -> dict[str, Any]:
    """VCR config for tests that make HTTP requests."""
    return {**DEFAULT_VCR_CONFIG}


"""Completely crap test for the EXTRA_ALGORITHMS constant"""


def test_extra_algorithms_pass_contains_only_crc64nvme() -> None:
    """EXTRA_ALGORITHMS is exactly the set of algorithms handled outside hashlib."""
    assert frozenset({"crc64nvme"}) == EXTRA_ALGORITHMS


"""Tests for the ChecksumEntry dataclass"""


def test_checksumentry_pass_construction_and_fields() -> None:
    """ChecksumEntry stores algorithm and value exactly as provided."""
    entry = ChecksumEntry(algorithm="sha256", value="flash_gordon")
    assert entry.algorithm == "sha256"
    assert entry.value == "flash_gordon"


def test_checksumentry_pass_equality() -> None:
    """Two ChecksumEntry instances with the same fields compare equal; differing fields compare unequal."""
    a = ChecksumEntry(algorithm="sha256", value="flash_gordon")
    b = ChecksumEntry(algorithm="sha256", value="flash_gordon")
    c = ChecksumEntry(algorithm="md5", value="flash_gordon")
    assert a == b
    assert a != c


def test_checksumentry_pass_hashable() -> None:
    """ChecksumEntry instances are hashable and usable in sets/dict keys, as implied by frozen=True."""
    a = ChecksumEntry(algorithm="sha256", value="flash_gordon")
    b = ChecksumEntry(algorithm="sha256", value="flash_gordon")
    assert hash(a) == hash(b)
    assert len({a, b}) == 1


def test_checksumentry_fail_frozen_immutable() -> None:
    """Mutating an existing field on a ChecksumEntry raises FrozenInstanceError."""
    entry = ChecksumEntry(algorithm="sha256", value="flash_gordon")
    with pytest.raises(dataclasses.FrozenInstanceError, match="cannot assign to field 'algorithm'"):
        entry.algorithm = "md5"  # type: ignore[misc]


def test_checksumentry_fail_frozen_prevents_new_attributes() -> None:
    """The frozen dataclass rejects attribute assignment even for names that aren't declared fields."""
    entry = ChecksumEntry(algorithm="sha256", value="flash_gordon")
    # can be one of two errors, depending on which version of CPython is in use.
    # should be a FrozenInstanceError but some versions produce a TypeError instead
    with pytest.raises(
        (TypeError, dataclasses.FrozenInstanceError),
        match=r"(is not an instance or subtype of type \(ChecksumEntry\)|cannot assign to field 'extra')",
    ):
        entry.extra = "not allowed"  # type: ignore[attr-defined]


"""Tests of new_hasher"""


@pytest.mark.parametrize("transformation", TEXT_TRANSFORMATIONS)
def test_new_hasher_pass_returns_hashlib_object_for_standard_algorithm(transformation: Callable) -> None:
    """new_hasher returns a hashlib-backed object whose digest matches hashlib directly."""
    h: Hasher = new_hasher(transformation("sha256"))
    assert isinstance(h, Hasher)
    assert not isinstance(h, Crc64NvmeHasher)
    assert h.hexdigest() == hashlib.sha256(b"").hexdigest()


@pytest.mark.parametrize("transformation", TEXT_TRANSFORMATIONS)
def test_new_hasher_pass_returns_crc64nvme_hasher(transformation: Callable) -> None:
    """new_hasher returns a Crc64NvmeHasher instance for the 'crc64nvme' algorithm."""
    h: Hasher = new_hasher(transformation("crc64nvme"))
    assert isinstance(h, Crc64NvmeHasher)
    assert isinstance(h, Hasher)


@pytest.mark.parametrize("transformation", TEXT_TRANSFORMATIONS)
@pytest.mark.parametrize("algorithm", UNSUPPORTED_ALGORITHMS)
def test_new_hasher_fail_unsupported_algorithm(algorithm: str, transformation: Callable) -> None:
    """new_hasher raises ValueError for an unsupported algorithm name."""
    with pytest.raises(ValueError, match="not supported"):
        new_hasher(transformation(algorithm))


@pytest.mark.parametrize("transformation", TEXT_TRANSFORMATIONS)
def test_new_hasher_fail_crc64nvme_include_extras_false(transformation: Callable) -> None:
    """new_hasher raises ValueError for the 'crc64nvme' algorithm when include_extras is False."""
    with pytest.raises(ValueError, match="not supported"):
        new_hasher(transformation("crc64nvme"), include_extras=False)


"""Tests for the Hasher Protocol"""


def test_hasher_protocol_fail_object_without_required_methods() -> None:
    """An object lacking update()/hexdigest() does not satisfy the Hasher protocol."""

    class NotAHasher:
        pass

    assert not isinstance(NotAHasher(), Hasher)


def test_hasher_protocol_pass_custom_duck_typed_object() -> None:
    """Any object exposing update() and hexdigest(), regardless of implementation, satisfies Hasher."""

    class CustomHasher:
        def update(self, data: bytes) -> None:
            pass

        def hexdigest(self) -> str:
            return "custom"

    assert isinstance(CustomHasher(), Hasher)


"""Tests of the Crc64NvmeHasher class"""


def test_crc64nvmehasher_pass_incremental_matches_single_shot() -> None:
    """Feeding data to Crc64NvmeHasher in small chunks matches a single update() call."""
    data: bytes = b"abcdefgh" * 1000
    incremental = Crc64NvmeHasher()
    for i in range(0, len(data), 7):
        incremental.update(data[i : i + 7])
    single = Crc64NvmeHasher()
    single.update(data)
    assert incremental.hexdigest() == single.hexdigest()
    assert incremental.b64digest() == single.b64digest()


def test_crc64nvmehasher_pass_digest_hexdigest_b64digest_consistent() -> None:
    """digest(), hexdigest(), and b64digest() are all consistent representations of the same value."""
    h = Crc64NvmeHasher()
    h.update(b"hello world")
    assert h.hexdigest() == h.digest().hex()
    assert h.b64digest() == base64.b64encode(h.digest()).decode()


def test_crc64nvmehasher_pass_empty_input_is_zero() -> None:
    """An unmodified Crc64NvmeHasher digests to eight zero bytes."""
    h = Crc64NvmeHasher()
    assert h.digest() == (0).to_bytes(8, "big")


def test_crc64nvmehasher_pass_update_returns_none() -> None:
    """update() returns None; the Hasher Protocol's `-> object` signature only requires this to be compatible."""
    h = Crc64NvmeHasher()
    assert h.update(b"data") is None


"""Tests for compute_file_checksum"""


@pytest.mark.parametrize("algorithm", HASHLIB_ALGORITHMS)
@pytest.mark.parametrize("transformation", TEXT_TRANSFORMATIONS)
@pytest.mark.parametrize("content", TEST_DATASETS)
def test_compute_file_checksum_pass_hashlib_algorithm_matches_hashlib(
    tmp_path: Path, algorithm: str, transformation: Callable, content: bytes
) -> None:
    """compute_file_checksum matches hashlib's own digest for supported hashlib algorithms."""
    file_path = tmp_path / "data.bin"
    file_path.write_bytes(content)
    expected = hashlib.new(algorithm, content).hexdigest()
    assert compute_file_checksum(file_path, transformation(algorithm)) == expected


@pytest.mark.parametrize("transformation", TEXT_TRANSFORMATIONS)
@pytest.mark.parametrize("content", TEST_DATASETS)
def test_compute_file_checksum_pass_crc64nvme_matches_direct_computation(
    tmp_path: Path, content: bytes, transformation: Callable
) -> None:
    """compute_file_checksum's crc64nvme result matches folding crc64nvme over the content directly.

    Relies on crc64nvme being a proper incremental CRC, i.e. chunk-boundary-independent.
    """
    file_path = tmp_path / "data.bin"
    file_path.write_bytes(content)
    expected_crc = _crc64nvme(content, 0)
    expected = base64.b64encode(expected_crc.to_bytes(8, byteorder="big")).decode()
    assert compute_file_checksum(file_path, transformation("crc64nvme")) == expected


@pytest.mark.parametrize("transformation", TEXT_TRANSFORMATIONS)
@pytest.mark.parametrize("algorithm", UNSUPPORTED_ALGORITHMS)
def test_compute_file_checksum_fail_unsupported_algorithm(
    tmp_path: Path, algorithm: str, transformation: Callable
) -> None:
    """An unsupported algorithm name, including shake* variants, raises ValueError before any file I/O happens."""
    file_path = tmp_path / "data.bin"
    file_path.write_bytes(b"content")
    with pytest.raises(ValueError, match="not supported"):
        compute_file_checksum(file_path, transformation(algorithm))


@pytest.mark.parametrize("path_type", [str, Path])
def test_compute_file_checksum_pass_accepts_str_or_path(tmp_path: Path, path_type: type) -> None:
    """compute_file_checksum accepts both str and Path arguments for the file location."""
    file_path = tmp_path / "data.bin"
    file_path.write_bytes(b"content")
    expected = hashlib.sha256(b"content").hexdigest()
    assert compute_file_checksum(path_type(file_path), "sha256") == expected


def test_compute_file_checksum_fail_missing_file(tmp_path: Path) -> None:
    """A nonexistent file path raises FileNotFoundError rather than being silently ignored."""
    missing_path = tmp_path / "does_not_exist.bin"
    with pytest.raises(FileNotFoundError, match="No such file or directory"):
        compute_file_checksum(missing_path, "sha256")


def test_compute_file_checksum_fail_directory_path(tmp_path: Path) -> None:
    """Passing a directory path raises IsADirectoryError instead of being silently mishandled."""
    with pytest.raises(IsADirectoryError, match="Is a directory"):
        compute_file_checksum(tmp_path, "sha256")


@pytest.mark.parametrize("size", [(1 << 20) - 1, 1 << 20, (1 << 20) + 1])
def test_compute_file_checksum_pass_chunk_boundary_sizes(tmp_path: Path, size: int) -> None:
    """Files whose size lands exactly on the internal 1 MiB read-chunk boundary hash correctly."""
    file_path: Path = tmp_path / "data.bin"
    content: bytes = b"y" * size
    file_path.write_bytes(content)
    assert compute_file_checksum(file_path, "sha256") == hashlib.sha256(content).hexdigest()


@pytest.mark.parametrize("size", [(1 << 20) - 1, 1 << 20, (1 << 20) + 1])
def test_compute_file_checksum_pass_crc64nvme_chunk_boundary_sizes(tmp_path: Path, size: int) -> None:
    """crc64nvme checksums remain correct for files whose size lands exactly on the internal 1 MiB read-chunk boundary."""
    file_path = tmp_path / "data.bin"
    content = b"z" * size
    file_path.write_bytes(content)
    expected_crc = _crc64nvme(content, 0)
    expected = base64.b64encode(expected_crc.to_bytes(8, byteorder="big")).decode()
    assert compute_file_checksum(file_path, "crc64nvme") == expected


"""Tests for validate_checksum_fn"""


@pytest.mark.parametrize("algorithm", HASHLIB_ALGORITHMS)
@pytest.mark.parametrize("transformation", TEXT_TRANSFORMATIONS)
def test_validate_checksum_fn_pass_supported_algorithm(algorithm: str, transformation: Callable) -> None:
    """Supported hashlib algorithms are accepted and returned lowercased."""
    assert validate_checksum_fn(transformation(algorithm)) == algorithm.lower()


@pytest.mark.parametrize("transformation", TEXT_TRANSFORMATIONS)
@pytest.mark.parametrize("algorithm", UNSUPPORTED_ALGORITHMS)
def test_validate_checksum_fn_fail_unsupported_algorithm(algorithm: str, transformation: Callable) -> None:
    """Unsupported algorithm names raise ValueError."""
    with pytest.raises(ValueError, match="not supported"):
        validate_checksum_fn(transformation(algorithm))


@pytest.mark.parametrize("transformation", TEXT_TRANSFORMATIONS)
def test_validate_checksum_fn_pass_crc64nvme_case_insensitive(transformation: Callable) -> None:
    """crc64nvme is recognised regardless of case, same as hashlib algorithm names."""
    assert validate_checksum_fn(transformation("crc64nvme")) == "crc64nvme"


@pytest.mark.parametrize("algorithm", HASHLIB_ALGORITHMS)
@pytest.mark.parametrize("transformation", TEXT_TRANSFORMATIONS)
def test_validate_checksum_fn_pass_supported_algorithm_include_extras_false(
    algorithm: str, transformation: Callable
) -> None:
    """Supported hashlib algorithms are accepted and returned lowercased with include_extras=False."""
    assert validate_checksum_fn(transformation(algorithm), include_extras=False) == algorithm.lower()


@pytest.mark.parametrize("transformation", TEXT_TRANSFORMATIONS)
def test_validate_checksum_fn_fail_crc64nvme_include_extras_false(transformation: Callable) -> None:
    """Unsupported algorithm names raise ValueError."""
    with pytest.raises(ValueError, match="not supported"):
        validate_checksum_fn(transformation("crc64nvme"), include_extras=False)


"""Tests for the HashingReader class"""


@pytest.mark.parametrize("data", [b"", b"hello world", b"x" * 5000])
def test_hashingreader_pass_hexdigest_matches_hashlib(data: bytes) -> None:
    """Reading all bytes through HashingReader produces the same digest as hashing directly."""
    reader = HashingReader(BytesIO(data), "sha256")
    read_data = reader.read(-1)
    assert read_data == data
    assert reader.hexdigest() == hashlib.sha256(data).hexdigest()


def test_hashingreader_pass_incremental_reads_accumulate_hash() -> None:
    """Multiple small reads update the running hash the same way as one big read."""
    data = b"a" * 10 + b"b" * 10 + b"c" * 10
    reader = HashingReader(BytesIO(data), "md5")
    chunks = [reader.read(10) for _ in range(3)]
    assert b"".join(chunks) == data
    assert reader.hexdigest() == hashlib.md5(data).hexdigest()  # noqa: S324


def test_hashingreader_fail_unsupported_algorithm() -> None:
    """Constructing a HashingReader with an unsupported algorithm raises ValueError."""
    with pytest.raises(ValueError, match="not supported"):
        HashingReader(BytesIO(b"data"), "not_a_real_algorithm")


def test_hashingreader_pass_proxies_unknown_attributes() -> None:
    """Attributes not defined on HashingReader are proxied through to the wrapped object."""
    raw = BytesIO(b"data")
    reader = HashingReader(raw, "md5")
    assert reader.tell() == raw.tell()


def test_hashingreader_fail_proxies_missing_attribute() -> None:
    """Accessing an attribute that exists on neither wrapper nor wrapped object raises AttributeError."""
    reader = HashingReader(BytesIO(b"data"), "md5")
    with pytest.raises(AttributeError, match="object has no attribute 'not_a_real_attribute'"):
        _ = reader.not_a_real_attribute


def test_hashingreader_pass_crc64nvme_b64digest_and_hexdigest_output() -> None:
    """HashingReader.b64digest() for crc64nvme returns S3's base64 ChecksumCRC64NVME format."""
    data = b"some bytes"
    reader = HashingReader(BytesIO(data), "crc64nvme")
    reader.read(-1)
    expected_crc: int = _crc64nvme(data, 0)
    # hexdigest returns the hex form
    hexpected = expected_crc.to_bytes(8, "big").hex()
    assert reader.hexdigest() == hexpected
    # b64digest returns S3's base64 form
    expected = base64.b64encode(expected_crc.to_bytes(8, byteorder="big")).decode()
    assert reader.b64digest() == expected


def test_hashingreader_fail_b64digest_unsupported_for_hashlib_algorithm() -> None:
    """b64digest() raises AttributeError for algorithms whose hasher doesn't implement it."""
    reader = HashingReader(BytesIO(b"data"), "sha256")
    with pytest.raises(AttributeError, match="object has no attribute 'b64digest'"):
        reader.b64digest()


def test_hashingreader_fail_getattr_before_init_does_not_recurse() -> None:
    """Accessing _raw/_hash before __init__ has run raises AttributeError, not RecursionError."""
    reader = HashingReader.__new__(HashingReader)
    with pytest.raises(AttributeError, match="object has no attribute '_raw'"):
        _ = reader.tell()  # any proxied access forces lookup of self._raw first


def test_hashingreader_pass_copy_does_not_recurse() -> None:
    """copy.copy on a HashingReader does not trigger infinite recursion via __getattr__."""
    reader = HashingReader(BytesIO(b"data"), "md5")
    duplicate = copy.copy(reader)
    assert duplicate.hexdigest() == reader.hexdigest()


def test_hashingreader_pass_proxies_method_with_arguments() -> None:
    """Proxied methods that take arguments (not just no-arg attributes) work through __getattr__."""
    raw = BytesIO(b"0123456789")
    reader = HashingReader(raw, "md5")
    reader.seek(5)
    assert raw.tell() == 5


def test_hashingreader_fail_context_manager_not_supported() -> None:
    """Dunder methods aren't proxied, so `with` doesn't work even though the wrapped object supports it."""
    reader = HashingReader(BytesIO(b"data"), "md5")
    with pytest.raises(TypeError, match="object does not support the context manager protocol"), reader:  # pyright: ignore[reportGeneralTypeIssues]
        pass


def test_hashingreader_pass_read_amt_zero_returns_empty_and_does_not_change_hash() -> None:
    """Reading zero bytes returns an empty chunk and leaves the running hash unchanged."""
    reader = HashingReader(BytesIO(b"data"), "md5")
    assert reader.read(0) == b""
    assert reader.hexdigest() == hashlib.md5(b"").hexdigest()  # noqa: S324


def test_hashingreader_fail_readline_bypasses_hashing() -> None:
    """readline() is proxied directly to the wrapped object and does NOT update the running hash.

    Documents a real gotcha: only read() is intercepted, so callers using
    readline()/readlines()/iteration will silently get an incomplete/incorrect digest.
    """
    data = b"line one\nline two\n"
    reader = HashingReader(BytesIO(data), "sha256")
    line = reader.readline()
    assert line == b"line one\n"
    # hash was NOT updated by readline(), so it still reflects zero bytes read
    assert reader.hexdigest() == hashlib.sha256(b"").hexdigest()


def test_hashingreader_pass_read_default_argument_reads_all() -> None:
    """Calling read() with no arguments reads all remaining data, same as read(-1)."""
    data = b"hello world"
    reader = HashingReader(BytesIO(data), "sha256")
    assert reader.read() == data
    assert reader.hexdigest() == hashlib.sha256(data).hexdigest()


def test_hashingreader_pass_progressive_digest_reflects_partial_reads() -> None:
    """hexdigest() can be called mid-stream, reflecting only bytes read so far, and updates with each subsequent read."""
    data = b"abcdefghij"
    reader = HashingReader(BytesIO(data), "sha256")
    first_half = reader.read(5)
    assert reader.hexdigest() == hashlib.sha256(first_half).hexdigest()
    second_half = reader.read(5)
    assert reader.hexdigest() == hashlib.sha256(first_half + second_half).hexdigest()


def test_hashingreader_pass_repeated_reads_past_eof_are_idempotent() -> None:
    """Once the wrapped stream is exhausted, further read() calls return empty bytes and leave the hash unchanged."""
    data = b"data"
    reader = HashingReader(BytesIO(data), "sha256")
    reader.read(-1)
    digest_after_exhaustion = reader.hexdigest()
    assert reader.read(-1) == b""
    assert reader.read(10) == b""
    assert reader.hexdigest() == digest_after_exhaustion == hashlib.sha256(data).hexdigest()


@pytest.mark.vcr
def test_hashing_reader_reads_real_http_response(requests_session: requests.Session) -> None:
    """HashingReader should stream and hash a real HTTP response body correctly."""
    resp = requests_session.get(TEST_ECHO_URL, stream=True)
    resp.raw.decode_content = True

    hasher = HashingReader(resp.raw, "sha256")
    body = hasher.read(-1)

    assert len(body) == 1024
    # sanity check: hash is deterministic and matches recomputation
    assert hasher.hexdigest() == hashlib.sha256(body).hexdigest()


@pytest.mark.vcr
def test_hashing_reader_handles_chunked_real_response(requests_session: requests.Session) -> None:
    """Hash computed via chunked reads must match hash from a single full read of the same cassette."""
    resp = requests_session.get(TEST_ECHO_URL, stream=True)
    resp.raw.decode_content = True

    hasher = HashingReader(resp.raw, "sha256")
    chunks = []
    while chunk := hasher.read(128):
        chunks.append(chunk)

    assert hasher.hexdigest() == hashlib.sha256(b"".join(chunks)).hexdigest()


@pytest.mark.vcr
def test_hashing_reader_handles_real_gzip_response(requests_session: requests.Session) -> None:
    """A real gzip-encoded response should be transparently decoded before hashing."""
    resp = requests_session.get("https://httpbin.org/gzip", stream=True)
    resp.raw.decode_content = True

    hasher = HashingReader(resp.raw, "sha256")
    body = hasher.read(-1)

    assert b'"gzipped": true' in body
