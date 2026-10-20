"""Draft-specific JSON Schema keyword selection for field metadata."""

from dataclasses import dataclass
from functools import cache
from typing import Final

from jsonschema import Draft7Validator

DEFAULT_METADATA_KEYWORDS: Final = frozenset({"title"})
REF_AND_IDENTITY_KEYWORDS: Final = frozenset(
    {"$anchor", "$defs", "$dynamicAnchor", "$dynamicRef", "$id", "$ref", "$schema", "$vocabulary", "definitions", "id"}
)
STRUCTURAL_OR_COMPOSITIONAL_KEYWORDS: Final = frozenset(
    {
        "additionalItems",
        "additionalProperties",
        "allOf",
        "anyOf",
        "contains",
        "contentSchema",
        "dependencies",
        "dependentRequired",
        "dependentSchemas",
        "else",
        "if",
        "items",
        "not",
        "oneOf",
        "patternProperties",
        "prefixItems",
        "properties",
        "propertyNames",
        "required",
        "then",
        "type",
        "unevaluatedItems",
        "unevaluatedProperties",
    }
)


def get_known_jsonschema_keywords(validator_cls: type) -> set[str]:
    """Combine draft-specific assertions with the Draft-07 annotation vocabulary."""
    keywords = set(Draft7Validator.META_SCHEMA["properties"]) | set(validator_cls.VALIDATORS)
    return keywords - REF_AND_IDENTITY_KEYWORDS


@cache
def metadata_keys_for(validator_cls: type) -> frozenset[str]:
    """Find nonstructural field metadata keywords for a validator class."""
    return frozenset(get_known_jsonschema_keywords(validator_cls) - STRUCTURAL_OR_COMPOSITIONAL_KEYWORDS)


@dataclass(frozen=True)
class ConversionContext:
    """Draft-specific metadata selection for a single emission."""

    validator_cls: type
    extra_metadata_keywords: frozenset[str] = frozenset()

    @property
    def metadata_keys(self) -> frozenset[str]:
        """Return the active draft's eligible metadata keywords."""
        return metadata_keys_for(self.validator_cls)

    @property
    def allowed_extra_metadata_keywords(self) -> frozenset[str]:
        """Select standard keywords and explicitly requested extension namespaces."""
        return frozenset(
            key for key in self.extra_metadata_keywords if key in self.metadata_keys or key.startswith("x-")
        )

    @property
    def invalid_extra_metadata_keywords(self) -> frozenset[str]:
        """Return ignored structural, identity and unknown keywords."""
        return self.extra_metadata_keywords - self.allowed_extra_metadata_keywords
