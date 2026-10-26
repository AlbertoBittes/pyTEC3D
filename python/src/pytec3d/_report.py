"""Text tables of the Python log, in the style of the C++ solver tables, and the run provenance.

Tables are printed on MPI rank 0 (standard output) and appended to the log file of the run, if any.
"""

from __future__ import annotations

import getpass
import hashlib
import os
import platform
import re
import socket
import sys
from collections.abc import Sequence
from datetime import datetime

from ._logging import world_rank, write_log_file

_NUMBER = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$")


class Table:
    """A text table: a title, optional column headers, rows and section titles; the columns fit their content."""

    def __init__(self, title: str) -> None:
        self._title = title
        self._header: list[str] = []
        self._lines: list[tuple[bool, list[str]]] = []  # (is_section, cells)

    def header(self, cells: Sequence[object]) -> Table:
        self._header = [str(c) for c in cells]
        return self

    def row(self, cells: Sequence[object]) -> Table:
        self._lines.append((False, [str(c) for c in cells]))
        return self

    def section(self, title: str) -> Table:
        self._lines.append((True, [title]))
        return self

    def __str__(self) -> str:
        rows = [cells for is_section, cells in self._lines if not is_section]
        columns = max([len(self._header), 1, *(len(r) for r in rows)])
        widths = [0] * columns
        for cells in [self._header, *rows]:
            for c, cell in enumerate(cells):
                widths[c] = max(widths[c], len(cell))
        inner = sum(widths) + 3 * (columns - 1)
        spanning = max([len(self._title), *(len(cells[0]) for is_section, cells in self._lines if is_section)])
        if spanning > inner:
            widths[-1] += spanning - inner
            inner = spanning

        rule = "+" + "+".join("-" * (w + 2) for w in widths) + "+"
        full_rule = "+" + "-" * (inner + 2) + "+"

        def spanning_line(text: str) -> str:
            return f"| {text.ljust(inner)} |"

        def row_line(cells: list[str]) -> str:
            padded = [cells[c] if c < len(cells) else "" for c in range(columns)]
            return "|" + "|".join(
                f" {cell.rjust(w) if _NUMBER.match(cell) else cell.ljust(w)} " for cell, w in zip(padded, widths)
            ) + "|"

        out = [full_rule, spanning_line(self._title), rule]
        if self._header:
            out += [row_line(self._header), rule]
        for i, (is_section, cells) in enumerate(self._lines):
            if is_section:
                if i > 0:
                    out.append(rule)
                out += [spanning_line(cells[0]), rule]
            else:
                out.append(row_line(cells))
        if not self._lines or not self._lines[-1][0]:
            out.append(rule)
        return "\n".join(out)


def print_table(table: Table | str) -> None:
    """Print a table on rank 0 and append it to the log file."""
    if world_rank() != 0:
        return
    text = str(table)
    sys.stdout.write(text + "\n")
    sys.stdout.flush()
    write_log_file(text)


def sha256(path: str) -> str:
    """The SHA-256 of a file, hexadecimal."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _user() -> str:
    try:
        return getpass.getuser()
    except (KeyError, OSError):
        return "-"


def provenance_table(input_files: dict[str, str | None], package_version: str) -> Table:
    """Where and how the run was made, and the checksums of its input files.

    Args:
        input_files: Label -> path; None or missing files are shown as such.
        package_version: The pytec3d version.
    """
    import numpy

    table = Table("Run")
    table.row(["Started", datetime.now().astimezone().isoformat(timespec="seconds")])
    table.row(["Command", " ".join([os.path.basename(sys.executable), *sys.argv])])
    table.row(["Working directory", os.getcwd()])
    table.row(["Host / user", f"{socket.gethostname()} / {_user()}"])
    table.row(["Platform", platform.platform()])
    table.row(["Python / numpy", f"{platform.python_version()} / {numpy.__version__}"])
    table.row(["pytec3d", package_version])
    table.section("Input files (SHA-256)")
    for label, path in input_files.items():
        if path is None:
            table.row([label, "-"])
        elif not os.path.isfile(path):
            table.row([label, f"{path} (not found)"])
        else:
            table.row([label, path])
            table.row(["", sha256(path)])
    return table
