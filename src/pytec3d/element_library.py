"""The element library: domains, fields and element types, read from the JSON resources of the solver.

``pytec3d/resources`` holds ``domains.json``, ``fields.json``, ``element_types.json`` and ``laws.json``. The C++ solver
compiles them (``src/tools/generate_registries.py``); this module reads the same files so that the Python side knows the
element types, their fields and the grid arrays a model needs. Set ``PYTEC3D_RESOURCES`` to read another directory.

Grid arrays of an element type (:func:`grid_arrays`):

* DOF fields: point array ``F`` (solution) and ``_F`` (Dirichlet values, NaN where free);
* nodal fields: point array ``F``; element fields: cell array ``F``;
* nodal loads (``mn.`` parameters): point array; element loads (``me.`` parameters): cell array.

A vector field has 3 components, a scalar one 1, whatever the active components of the element type.
"""

from __future__ import annotations

import json
import os
from collections.abc import Iterable
from dataclasses import dataclass
from functools import cache
from pathlib import Path
from typing import Any

from ._typing import Location

#: Element type names or indices (of ``element_types.json``); None for all of them.
ElementTypes = Iterable[str | int] | None


@dataclass(frozen=True)
class GridArray:
    """A grid array a model needs.

    Attributes:
        name: Array name.
        location: ``"node"`` (point data) or ``"elem"`` (cell data).
        components: 1 for a scalar, 3 for a vector.
        fill: Default value: 0, or NaN for the Dirichlet arrays ``_F``.
    """

    name: str
    location: Location
    components: int
    fill: float


def resources_dir() -> Path:
    """Directory of the JSON resources: ``$PYTEC3D_RESOURCES``, or ``resources`` inside the package."""
    override = os.environ.get("PYTEC3D_RESOURCES")
    return Path(override) if override else Path(__file__).resolve().parent / "resources"


def _load(name: str) -> dict[str, Any]:
    path = resources_dir() / name
    if not path.is_file():
        raise FileNotFoundError(f"Element library resource '{path}' not found (set PYTEC3D_RESOURCES).")
    with path.open(encoding="utf-8") as handle:
        data: dict[str, Any] = json.load(handle)
    return data


@cache
def domains() -> tuple[str, ...]:
    """Physical regions (``domains.json``), in index order."""
    return tuple(_load("domains.json")["domains"])


@cache
def fields() -> dict[str, str]:
    """Field name -> ``"scalar"`` or ``"vector"`` (``fields.json``)."""
    return dict(_load("fields.json"))


@cache
def element_types() -> dict[str, dict[str, Any]]:
    """Element type name -> definition, with ``inherits`` resolved, in index order (``element_types.json``).

    Raises:
        ValueError: If the ``inherits`` chain of an element type is circular.
    """
    definitions: dict[str, dict[str, Any]] = _load("element_types.json")["element_types"]
    resolved: dict[str, dict[str, Any]] = {}

    def resolve(name: str, chain: tuple[str, ...]) -> dict[str, Any]:
        if name in resolved:
            return resolved[name]
        if name in chain:
            raise ValueError(f"element type '{name}': circular inherits ({' -> '.join(chain + (name,))}).")
        definition = definitions[name]
        parent = definition.get("inherits")
        # A section given here replaces the parent's entries key by key.
        result: dict[str, Any] = (
            {} if parent is None
            else {k: dict(v) if isinstance(v, dict) else v for k, v in resolve(parent, chain + (name,)).items()}
        )
        for key, value in definition.items():
            if key == "inherits":
                continue
            if isinstance(value, dict) and isinstance(result.get(key), dict):
                result[key].update(value)
            else:
                result[key] = value
        resolved[name] = result
        return result

    return {name: resolve(name, ()) for name in definitions}


def element_type_names() -> list[str]:
    """Element type names; the index of a name is its ``element_type`` index in the model data."""
    return list(element_types())


def element_type_name(value: str | int) -> str:
    """The name of an element type given by name or by index.

    Raises:
        ValueError: If the name is unknown or the index out of range.
    """
    names = element_type_names()
    if isinstance(value, str):
        if value not in names:
            raise ValueError(f"Unknown element type {value!r}; known: {names}.")
        return value
    index = int(value)
    if not 0 <= index < len(names):
        raise ValueError(f"Element type index {index} out of range [0, {len(names)}).")
    return names[index]


def _names(types: ElementTypes) -> list[str]:
    return element_type_names() if types is None else [element_type_name(value) for value in types]


def _parameters(definition: dict[str, Any]) -> Iterable[tuple[str, str]]:
    """(prefix, name) of each declared parameter, e.g. ("mn", "AMPS") for ``"mn.AMPS"``."""
    for declaration in definition.get("parameters", {}):
        prefix, _, name = "".join(declaration.split()).partition(".")
        yield prefix, name


def _components(field: str) -> int:
    return 3 if fields().get(field) == "vector" else 1


def dof_fields(types: ElementTypes = None) -> list[str]:
    """DOF field names of the element types (all of the library when None), each once, in library order."""
    out: dict[str, None] = {}
    for name in _names(types):
        out.update(dict.fromkeys(element_types()[name].get("dof_fields", {})))
    return list(out)


def grid_arrays(types: ElementTypes = None) -> list[GridArray]:
    """Grid arrays the element types (all of the library when None) read and write, each once."""
    arrays: dict[tuple[str, str], GridArray] = {}

    def add(name: str, location: Location, fill: float, field: str | None = None) -> None:
        arrays.setdefault((location, name), GridArray(name, location, _components(field or name), fill))

    for type_name in _names(types):
        definition = element_types()[type_name]
        for field in definition.get("dof_fields", {}):
            add(field, "node", 0.0)
            add("_" + field, "node", float("nan"), field)
        for field in definition.get("nodal_fields", {}):
            add(field, "node", 0.0)
        for field in definition.get("element_fields", {}):
            add(field, "elem", 0.0)
        for prefix, name in _parameters(definition):
            if prefix == "mn":
                add(name, "node", 0.0)
            elif prefix == "me":
                add(name, "elem", 0.0)
    return list(arrays.values())


def load_arrays(types: ElementTypes = None) -> list[GridArray]:
    """Load arrays of the element types (all of the library when None): their ``me.`` and ``mn.`` parameters."""
    arrays: dict[tuple[str, str], GridArray] = {}
    for type_name in _names(types):
        for prefix, name in _parameters(element_types()[type_name]):
            if prefix in ("me", "mn"):
                location: Location = "elem" if prefix == "me" else "node"
                arrays.setdefault((location, name), GridArray(name, location, _components(name), 0.0))
    return list(arrays.values())
