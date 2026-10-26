"""Exceptions raised by pytec3d.

Every exception of the package derives from :class:`PyTEC3DError`, so ``except PyTEC3DError`` catches them all. Each
also derives from the built-in exception that describes it (``ValueError``, ``RuntimeError``), so code written against
the built-in types keeps working.
"""


class PyTEC3DError(Exception):
    """Base class of the pytec3d exceptions."""


class ValidationError(PyTEC3DError, ValueError):
    """The model, its data, the solver settings or the time step are inconsistent.

    The message lists every problem found, one per line.
    """

    def __init__(self, title: str, problems: list[str]) -> None:
        self.problems = list(problems)
        super().__init__(f"{title} ({len(self.problems)} problem(s)):\n" + "\n".join(f"  - {p}" for p in self.problems))


class MissingArrayError(PyTEC3DError, ValueError):
    """A point or cell array the operation needs is not on the mesh."""

    def __init__(self, name: str, location: str, hint: str = "") -> None:
        self.name = name
        self.location = location
        where = "point" if location == "node" else "cell"
        super().__init__(f"The mesh has no {where} array '{name}'." + (f" {hint}" if hint else ""))


class ModelLockedError(PyTEC3DError, RuntimeError):
    """The model is loaded into a solver session: its mesh and model data cannot change."""


class SolverStateError(PyTEC3DError, RuntimeError):
    """An operation was called in the wrong state (no mesh, no solver session, closed solver)."""


class ConvergenceError(PyTEC3DError, RuntimeError):
    """An iterative procedure did not converge within its iteration limit."""
