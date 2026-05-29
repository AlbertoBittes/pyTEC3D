import logging
import os
from datetime import datetime
from typing import Optional


_LOGGER_CONFIGURED = False


class _CppLikeFormatter(logging.Formatter):
    """Format Python logs to mirror the C++ logger style."""

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
        self._use_color = bool(use_color)

    def format(self, record: logging.LogRecord) -> str:
        level = record.levelname.lower()
        message = record.getMessage()
        ts = datetime.fromtimestamp(record.created).strftime("%Y-%m-%d %H:%M:%S")
        timestamp = f"{ts}.{int(record.msecs):03d}"
        prefix = f"[{level}] "

        if self._use_color:
            color = self._LEVEL_COLORS.get(record.levelno)
            if color is not None:
                prefix = f"{color}[{level}]{self._RESET} "

        rendered = f"{timestamp} {prefix}{message}"
        if record.exc_info:
            rendered = f"{rendered}\n{self.formatException(record.exc_info)}"
        if record.stack_info:
            rendered = f"{rendered}\n{self.formatStack(record.stack_info)}"
        return rendered


def _stream_supports_color(stream: object) -> bool:
    """Return True when ANSI color should be used for this output stream."""
    if os.environ.get("NO_COLOR") is not None:
        return False

    isatty = getattr(stream, "isatty", None)
    if callable(isatty):
        try:
            return bool(isatty())
        except Exception:
            return False
    return False


def _safe_world_rank() -> int:
    """Return MPI world rank when available.

    Returns:
        int: COMM_WORLD rank when mpi4py is available; otherwise 0.
    """
    try:
        from mpi4py import MPI
    except Exception:
        return 0

    try:
        return int(MPI.COMM_WORLD.Get_rank())
    except Exception:
        return 0


def _log_level_from_env(default: str = "INFO") -> int:
    """Resolve package log level from environment with a fallback default.

    Args:
        default: Fallback log-level name used when environment is not set.

    Returns:
        int: Python logging level numeric value.
    """
    if not isinstance(default, str) or not default.strip():
        raise ValueError("default must be a non-empty string log level name.")
    level_name = os.environ.get("APEC_M_LOG_LEVEL", default).upper()
    return getattr(logging, level_name, logging.INFO)


def configure_logging(level: Optional[int] = None, root_only: Optional[bool] = None) -> None:
    """Configure package logging with root-rank-only emission semantics.

    Args:
        level: Optional explicit numeric logging level.
        root_only: Retained for API compatibility; logs are emitted only on rank 0.

    Returns:
        None
    """
    global _LOGGER_CONFIGURED
    if _LOGGER_CONFIGURED:
        return

    if level is None:
        level = _log_level_from_env()
    elif not isinstance(level, int):
        raise TypeError("level must be an int logging level or None.")

    if root_only is not None and not isinstance(root_only, bool):
        raise TypeError("root_only must be a bool or None.")

    # Python logs should mirror the C++ behavior: only rank 0 emits output.
    class _RootOnlyFilter(logging.Filter):
        def filter(self, record: logging.LogRecord) -> bool:
            return _safe_world_rank() == 0

    handler = logging.StreamHandler()
    handler.setLevel(level)
    handler.addFilter(_RootOnlyFilter())
    handler.setFormatter(_CppLikeFormatter(_stream_supports_color(handler.stream)))

    package_logger = logging.getLogger("apec_m")
    package_logger.setLevel(level)
    package_logger.handlers.clear()
    package_logger.addHandler(handler)
    package_logger.propagate = False

    _LOGGER_CONFIGURED = True


def get_logger(name: str) -> logging.Logger:
    """Return a configured package logger by name.

    Args:
        name: Logger name.

    Returns:
        logging.Logger: Configured logger instance.
    """
    if not isinstance(name, str) or not name.strip():
        raise ValueError("name must be a non-empty logger name.")
    configure_logging()
    return logging.getLogger(name)
