"""Path-related fields and utils for pipeline settings."""

from pathlib import Path
from typing import Annotated, Any, Final

from fsspec.core import split_protocol
from pydantic import AfterValidator, BeforeValidator, DirectoryPath, Field

VALID_PATH_ONE_CHAR: Final[str] = "A valid path must contain at least one character"
FILE_PROTOCOL: Final[str] = "file"
PATH_EMPTY: Final[str] = "Path must not be empty."
PATH_NOT_LOCAL: Final[str] = "Path must be on the local filesystem, not a URL"
DIR_NOT_CREATED: Final[str] = "Cannot create directory "


def url_protocol(url: str) -> str:
    """Return the canonical protocol of a path or URL. Plain paths count as 'file'.

    :param url: path or URL
    :type url: str
    :return: lower-case protocol name, with aliases (e.g. s3a -> s3) applied
    :rtype: str
    """
    protocol, _ = split_protocol(url)
    return (protocol or "file").lower()


def is_local(url: str) -> bool:
    """Whether a path or URL is on the local filesystem."""
    return url_protocol(url) == FILE_PROTOCOL


def check_local_path(value: Any) -> Any:  # noqa: ANN401
    """Strip whitespace and any local protocol from a string path.

    Values that are not strings are returned unchanged.

    :param value: path to check
    :type value: Any
    :raises ValueError: if the path is empty or not on the local filesystem
    :return: the path without a protocol
    :rtype: Any
    """
    if not isinstance(value, str):
        # pydantic's Path validation handles these
        return value
    value = value.strip()
    if not is_local(value):
        err_msg = f"{PATH_NOT_LOCAL}: {value!r}"
        raise ValueError(err_msg)
    _, path = split_protocol(value)
    # check after removing the protocol so that "file://" is rejected; Path("") would become Path(".")
    if not path:
        raise ValueError(PATH_EMPTY)
    return path


def to_absolute(path: Path) -> Path:
    """Return the absolute, resolved form of a path.

    :param path: the path to resolve
    :type  path: Path
    :return: absolute path
    :rtype : Path
    """
    return path.resolve()


def ensure_dir(path: Path) -> Path:
    """Create a directory, converting OS errors to ValueError so pydantic reports them.

    :param path: the directory to create
    :type  path: Path
    :return: the path object
    :rtype : Path
    """
    try:
        path.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        msg = f"{DIR_NOT_CREATED}{path!r}: {e}"
        raise ValueError(msg) from e
    return path


def normalise_dir(full_path: str | None) -> str:
    """Remove trailing slashes from a directory path or URL, keeping roots such as '/' and 'file:///'.

        Does not permit protocol-only paths, such as `https://`, `data:`, or `file://`

    :param full_path: directory path or URL; may include the protocol
    :type full_path: str | None
    :raises ValueError: if the path is empty or only contains a protocol
    :return: the normalised path
    :rtype : str
    """
    if full_path is None or not full_path.strip():
        raise ValueError(VALID_PATH_ONE_CHAR)
    full_path = full_path.strip()
    protocol, path = split_protocol(full_path)
    if not path:
        raise ValueError(PATH_EMPTY)
    stripped = path.rstrip("/")
    # path was '/' or multiple '/'s
    if path and not stripped:
        stripped = "/"
    sep = ""
    if protocol:
        sep = ":" if protocol == "data" else "://"
    return f"{protocol or ''}{sep}{stripped}"


LocalInputDir = Annotated[
    DirectoryPath,
    Field(
        description="Local directory containing file(s) to import",
    ),
    BeforeValidator(check_local_path),
    AfterValidator(to_absolute),
]

LocalOutputDir = Annotated[
    Path,
    Field(
        description="Local directory to save output and, if use_output_dir_for_pipeline_metadata is set, pipeline metadata to. Required. Use `use_destination` to set where dlt loads data.",
    ),
    BeforeValidator(check_local_path),
    AfterValidator(to_absolute),
    AfterValidator(ensure_dir),
]
