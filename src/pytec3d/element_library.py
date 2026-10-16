"""The element library of the solver, read from the JSON resources it is compiled from.

``src/resources`` holds ``domains.json``, ``fields.json`` and ``element_types.json``. The C++ solver compiles
them (``src/tools/generate_registries.py``); this module reads the same files to know, on the Python side, the
element types, their fields and the grid arrays a model needs. Set ``PYTEC3D_RESOURCES`` to use another
directory.

Grid arrays of an element type (``grid_arrays``):

* DOF fields: point array ``F`` (solution) and ``_F`` (Dirichlet values, NaN where free);
* nodal fields: point array ``F``; element fields: cell array ``F``;
* nodal loads (``mn.`` parameters): point array; element loads (``me.`` parameters): cell array.

A vector field has 3 components, a scalar one 1, whatever the active components of the element type.
"""
from __future__ import annotations

import json
import os
from collections import OrderedDict
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Optional

#: Sections of an element type that "inherits" replaces key by key.
_FIELD_SECTIONS = ("dof_fields", "nodal_fields", "element_fields")


@dataclass(frozen=True)
class GridArray:
    """A grid array a model needs: where it lives, its components and its default value."""

    name: str
    location: str  # "node" or "elem"
    components: int
    fill: float


def resources_dir() -> Path:
    """Directory of the JSON resources: ``PYTEC3D_RESOURCES`` or ``src/resources`` next to the package."""
    env = os.environ.get("PYTEC3D_RESOURCES")
    if env:
        return Path(env)
    return Path(__file__).resolve().parents[1] / "resources"


def _load(name: str) -> dict:
    path = resources_dir() / name
    if not path.exists():
        raise FileNotFoundError(f"Element library resource '{path}' not found (set PYTEC3D_RESOURCES).")
    with open(path, encoding="utf-8") as handle:
        return json.load(handle, object_pairs_hook=OrderedDict)


@lru_cache(maxsize=None)
def domains() -> tuple[str, ...]:
    """Physical regions (domains.json), in index order."""
    return tuple(_load("domains.json")["domains"])


@lru_cache(maxsize=None)
def fields() -> "OrderedDict[str, str]":
    """Field name -> "scalar" or "vector" (fields.json)."""
    return _load("fields.json")


@lru_cache(maxsize=None)
def element_types() -> "OrderedDict[str, dict]":
    """Element type name -> definition with "inherits" resolved, in index order (element_types.json)."""
    definitions = _load("element_types.json")["element_types"]
    resolved: "OrderedDict[str, dict]" = OrderedDict()

    def resolve(name: str, chain: tuple[str, ...]) -> dict:
        if name in resolved:
            return resolved[name]
        if name in chain:
            raise ValueError(f"element type '{name}': circular inherits.")
        definition = definitions[name]
        parent = definition.get("inherits")
        result = OrderedDict() if parent is None else OrderedDict(
            (k, OrderedDict(v) if isinstance(v, dict) else v) for k, v in resolve(parent, chain + (name,)).items()
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

    for name in definitions:
        resolve(name, ())
    return OrderedDict((name, resolved[name]) for name in definitions)


def element_type_names() -> list[str]:
    """Element type names; the index of a name is its element_type index in the model data."""
    return list(element_types())


def element_type_name(value) -> str:
    """The name of an element type given by name or by index."""
    names = element_type_names()
    if isinstance(value, str):
        if value not in names:
            raise ValueError(f"Unknown element type {value!r}; known: {names}.")
        return value
    index = int(value)
    if not 0 <= index < len(names):
        raise ValueError(f"Element type index {index} out of range [0, {len(names)}).")
    return names[index]


def dof_fields(type_names: Optional[Iterable] = None) -> list[str]:
    """DOF field names of the element types (all of them when None), without repetition."""
    out: list[str] = []
    for name in _types(type_names):
        for field in element_types()[name].get("dof_fields", {}):
            if field not in out:
                out.append(field)
    return out


def grid_arrays(type_names: Optional[Iterable] = None) -> list[GridArray]:
    """Grid arrays the element types (all of them when None) read and write."""
    kinds = fields()
    arrays: "OrderedDict[tuple[str, str], GridArray]" = OrderedDict()

    def add(name: str, location: str, fill: float, field: Optional[str] = None) -> None:
        components = 3 if kinds.get(field or name) == "vector" else 1
        arrays.setdefault((location, name), GridArray(name, location, components, fill))

    for type_name in _types(type_names):
        definition = element_types()[type_name]
        for field in definition.get("dof_fields", {}):
            add(field, "node", 0.0)
            add("_" + field, "node", float("nan"), field)
        for field in definition.get("nodal_fields", {}):
            add(field, "node", 0.0)
        for field in definition.get("element_fields", {}):
            add(field, "elem", 0.0)
        for declaration in definition.get("parameters", {}):
            prefix, _, name = "".join(declaration.split()).partition(".")
            if prefix == "mn":
                add(name, "node", 0.0)
            elif prefix == "me":
                add(name, "elem", 0.0)
    return list(arrays.values())


def load_arrays(type_names: Optional[Iterable] = None) -> list[GridArray]:
    """Load arrays of the element types (all of them when None): their me. (cell) and mn. (point) parameters."""
    out: "OrderedDict[tuple[str, str], GridArray]" = OrderedDict()
    for type_name in _types(type_names):
        for declaration in element_types()[type_name].get("parameters", {}):
            prefix, _, name = "".join(declaration.split()).partition(".")
            if prefix in ("me", "mn"):
                location = "elem" if prefix == "me" else "node"
                components = 3 if fields().get(name) == "vector" else 1
                out.setdefault((location, name), GridArray(name, location, components, 0.0))
    return list(out.values())


def _types(type_names: Optional[Iterable]) -> list[str]:
    if type_names is None:
        return element_type_names()
    return [element_type_name(value) for value in type_names]
