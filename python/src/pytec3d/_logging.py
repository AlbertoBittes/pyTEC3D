"""Logging of the package: one ``pytec3d`` logger, in the style of the C++ solver log, on MPI rank 0 only.

The level is ``INFO`` unless the environment variable ``PYTEC3D_LOG_LEVEL`` names another one (``DEBUG``,
``WARNING``...). Set ``NO_COLOR`` to disable the colored level labels.
"""

from __future__ import annotations

import logging
import os
from datetime import datetime
from functools import cache

#: Name of the package logger; every module logs to a child of it.
ROOT_LOGGER_NAME = "pytec3d"

#: Environment variables read for the log level, the first set wins.
_LEVEL_VARIABLES = ("PYTEC3D_LOG_LEVEL", "pyTEC3D_LOG_LEVEL")


@cache
def world_rank() -> int:
    """MPI rank of this process in ``COMM_WORLD``; 0 when mpi4py is not available.

    The rank cannot change during a run, so it is computed once.
    """
    try:
        from mpi4py import MPI
    except ImportError:
        return 0
    return int(MPI.COMM_WORLD.Get_rank())


class _SolverStyleFormatter(logging.Formatter):
    """``2026-01-01 12:00:00.000 [info] message``, the format of the C++ solver log."""

    _RESET = "\x1b[0m"
    _LEVEL_COLORS = {
        logging.DEBUG: "\x1b[36m",
        logging.INFO: "\x1b[32m",
        logging.WARNING: "\x1b[33m",
        logging.ERROR: "\x1b[31m",
        logging.CRITICAL: "\x1b[1;31m",
    }

    def __init__(self, use_color: bool) -> None:
        super().__init__()
        self._use_color = use_color

    def format(self, record: logging.LogRecord) -> str:
        timestamp = datetime.fromtimestamp(record.created).strftime("%Y-%m-%d %H:%M:%S")
        label = f"[{record.levelname.lower()}]"
        color = self._LEVEL_COLORS.get(record.levelno) if self._use_color else None
        if color is not None:
            label = f"{color}{label}{self._RESET}"
        rendered = f"{timestamp}.{int(record.msecs):03d} {label} {record.getMessage()}"
        if record.exc_info:
            rendered += "\n" + self.formatException(record.exc_info)
        if record.stack_info:
            rendered += "\n" + self.formatStack(record.stack_info)
        return rendered


class _RootRankFilter(logging.Filter):
    """Lets records through on MPI rank 0 only, as the C++ solver does."""

    def filter(self, record: logging.LogRecord) -> bool:
        return world_rank() == 0


class _CountingHandler(logging.Handler):
    """Counts the warnings and errors logged (on rank 0), for the end-of-run summary."""

    def __init__(self) -> None:
        super().__init__(logging.WARNING)
        self.warnings = 0
        self.errors = 0

    def emit(self, record: logging.LogRecord) -> None:
        if record.levelno >= logging.ERROR:
            self.errors += 1
        else:
            self.warnings += 1


_COUNTER = _CountingHandler()
_COUNTER.addFilter(_RootRankFilter())
_file_handler: logging.FileHandler | None = None


def _supports_color(stream: object) -> bool:
    if os.environ.get("NO_COLOR") is not None:
        return False
    isatty = getattr(stream, "isatty", None)
    try:
        return bool(isatty()) if callable(isatty) else False
    except (OSError, ValueError):
        return False


def _level_from_environment() -> int:
    for variable in _LEVEL_VARIABLES:
        name = os.environ.get(variable)
        if name:
            level = logging.getLevelName(name.strip().upper())
            if isinstance(level, int):
                return level
    return logging.INFO


@cache
def _configure() -> None:
    """Attach the handler to the package logger, once."""
    level = _level_from_environment()
    handler = logging.StreamHandler()
    handler.setLevel(level)
    handler.addFilter(_RootRankFilter())
    handler.setFormatter(_SolverStyleFormatter(_supports_color(handler.stream)))

    logger = logging.getLogger(ROOT_LOGGER_NAME)
    logger.setLevel(level)
    logger.handlers.clear()
    logger.addHandler(handler)
    logger.addHandler(_COUNTER)
    logger.propagate = False


def set_log_level(level: int | str) -> None:
    """Change the level of the package log.

    Args:
        level: A ``logging`` level (``logging.DEBUG``) or its name (``"DEBUG"``).

    Raises:
        ValueError: If ``level`` is not a known level name.
    """
    _configure()
    if isinstance(level, str):
        resolved = logging.getLevelName(level.strip().upper())
        if not isinstance(resolved, int):
            raise ValueError(f"Unknown log level {level!r}.")
        level = resolved
    logger = logging.getLogger(ROOT_LOGGER_NAME)
    logger.setLevel(level)
    for handler in logger.handlers:
        if handler is not _COUNTER:
            handler.setLevel(level)


def message_counts() -> tuple[int, int]:
    """``(warnings, errors)`` logged by the package so far (on rank 0; 0 on the other ranks)."""
    return _COUNTER.warnings, _COUNTER.errors


def open_log_file(path: str) -> None:
    """Also write the package log to ``path`` (appending), on rank 0; replaces the previous log file, if any."""
    global _file_handler
    _configure()
    close_log_file()
    if world_rank() != 0:
        return
    handler = logging.FileHandler(path, mode="a", encoding="utf-8")
    handler.setLevel(logging.getLogger(ROOT_LOGGER_NAME).level)
    handler.setFormatter(_SolverStyleFormatter(use_color=False))
    logging.getLogger(ROOT_LOGGER_NAME).addHandler(handler)
    _file_handler = handler


def close_log_file() -> None:
    """Stop writing the package log to a file."""
    global _file_handler
    if _file_handler is None:
        return
    logging.getLogger(ROOT_LOGGER_NAME).removeHandler(_file_handler)
    _file_handler.close()
    _file_handler = None


def write_log_file(text: str) -> None:
    """Append ``text`` (a table) to the log file as it is, if there is one."""
    handler = _file_handler
    if handler is None or handler.stream is None:
        return
    handler.acquire()
    try:
        handler.stream.write(text + "\n")
        handler.flush()
    finally:
        handler.release()


def get_logger(name: str) -> logging.Logger:
    """The logger ``name``, a child of the package logger (configured on first use).

    Args:
        name: Dotted logger name, e.g. ``"pytec3d.model"``.
    """
    _configure()
    return logging.getLogger(name)
