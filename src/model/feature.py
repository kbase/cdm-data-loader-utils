"""Class for representing feature data."""


class Feature:
    def __init__(
        self: "Feature",
        feature_id: str,
        sequence: str,
        description: str | None = None,
        aliases=None,
    ) -> None:
        self.id = feature_id
        self.seq = sequence
        self.description = description
        self.ontology_terms = {}
        self.aliases = aliases

    def add_ontology_term(self: "Feature", ontology_term: str, value: str) -> None:
        if ontology_term not in self.ontology_terms:
            self.ontology_terms[ontology_term] = []
        if value not in self.ontology_terms[ontology_term]:
            self.ontology_terms[ontology_term].append(value)
