"""Functions for reading XSD files."""

from pathlib import Path

import xmlschema


def load_schema(schema_path: str | Path) -> xmlschema.XMLSchema:
    """Load an XML schema document into an XMLSchema object.

    :param schema_path: path to the schema file
    :type schema_path: str | Path
    :return: loaded schema
    :rtype: xmlschema.XMLSchema
    """
    return xmlschema.XMLSchema(str(schema_path))


def find_list_and_single_child_paths(
    schema: xmlschema.XMLSchema,
) -> tuple[set[tuple[str, str]], set[tuple[str, str]]]:
    """Traverse an xmlschema XMLSchema and classify parent/child element paths.

    A child is put in list_paths if it can occur zero, one, or more times
    (maxOccurs is unbounded or greater than 1), meaning it should be
    represented as a list. Otherwise it goes in single_paths (occurs at most
    once, whether required or optional).

        :param schema: the schema to traverse
        :type schema: xmlschema.XMLSchema
        :return: tuple of (list_paths, single_paths)
        :rtype: tuple[set[tuple[str, str]], set[tuple[str, str]]]
    """
    list_paths: set[tuple[str, str]] = set()
    single_paths: set[tuple[str, str]] = set()
    visiting: set[int] = set()

    def walk(element: xmlschema.XsdElement) -> None:
        elem_type = element.type
        if elem_type is None:
            return

        model_group = elem_type.model_group
        if model_group is None:
            return

        # guards against infinite recursion for self-referential types
        type_key = id(elem_type)
        if type_key in visiting:
            return
        visiting.add(type_key)

        parent_name = element.local_name
        for child in model_group.iter_elements():
            if not isinstance(child, xmlschema.XsdElement):
                continue

            max_occurs = elem_type.overall_max_occurs(child)
            # never occurs => skip it
            if max_occurs == 0:
                continue

            path = (parent_name, child.local_name)
            if max_occurs is None or max_occurs > 1:
                list_paths.add(path)
            else:
                single_paths.add(path)

            walk(child)

        visiting.discard(type_key)

    for element in schema.elements.values():
        walk(element)

    return list_paths, single_paths
