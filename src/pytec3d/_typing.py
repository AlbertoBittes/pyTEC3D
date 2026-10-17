"""Type aliases shared by the modules of the package."""

from __future__ import annotations

from os import PathLike
from typing import Literal

#: Where an array lives: ``"node"`` for point data, ``"elem"`` for cell data.
Location = Literal["node", "elem"]

#: How a new selection combines with the current one (the APDL selection modes):
#: ``"s"`` replace it, ``"a"`` add to it, ``"r"`` reselect within it, ``"u"`` remove from it.
SelectionMode = Literal["s", "a", "r", "u"]

#: A filesystem path.
PathType = str | PathLike[str]

LOCATIONS: tuple[Location, ...] = ("node", "elem")
SELECTION_MODES: tuple[SelectionMode, ...] = ("s", "a", "r", "u")
