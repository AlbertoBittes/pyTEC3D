import logging
import os
from datetime import datetime
from typing import Any, Optional


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


def log_python_stdio_redirect(logger: logging.Logger, rank: int) -> None:
    """Log stdout/stderr redirection for non-root MPI ranks."""
    logger.debug("Python stdout/stderr redirected to os.devnull for rank %d.", rank)


def log_solver_shutdown_import_error(logger: logging.Logger, exc: Exception) -> None:
    """Log failure to import the solver module during shutdown."""
    logger.error("Failed to import solver module for shutdown: %s", exc)


def log_solver_runtime_shutdown_failure(logger: logging.Logger, exc: Exception) -> None:
    """Log failure while shutting down the C++ solver runtime."""
    logger.exception("Solver runtime shutdown failed: %s", exc)


def log_mesh_summary_info(
    logger: logging.Logger,
    mesh: Any,
    mesh_name: str,
    mesh_path: Optional[str],
    cwd: str,
) -> None:
    """Emit a formatted mesh summary at INFO level."""
    logger.info(format_mesh_summary(mesh, mesh_name, mesh_path, cwd))


def log_paraview_console_package_root(logger: logging.Logger, package_parent: str) -> None:
    """Log the package root configured for ParaView console usage."""
    logger.debug("ParaView console package root set to %s", package_parent)


def log_paraview_init(logger: logging.Logger) -> None:
    """Log ParaView-side APEC_M initialization start."""
    logger.info("Initializing APEC_M in ParaView.")


def log_paraview_working_directory(logger: logging.Logger, workdir: str) -> None:
    """Log the working directory selected for ParaView sessions."""
    logger.info("APEC_M ParaView working directory: %s", workdir)


def format_mesh_summary(
    mesh: Any,
    mesh_name: str,
    mesh_path: Optional[str],
    cwd: str,
) -> str:
    """Render a deterministic mesh summary string in the package log style."""
    w = 70
    sep = " " + "-" * w

    n_pts = mesh.GetNumberOfPoints()
    n_cells = mesh.GetNumberOfCells()
    bounds = mesh.GetBounds()
    z_span = bounds[5] - bounds[4]
    spatial_dim = "2D" if z_span < 1e-10 else "3D"

    def _bstr(lo: float, hi: float) -> str:
        return f"[{lo:.6g}, {hi:.6g}]"

    def _array_lines(data_obj: Any, indent: str = "    ") -> list[str]:
        rows: list[str] = []
        for i in range(data_obj.GetNumberOfArrays()):
            arr = data_obj.GetArray(i)
            if arr is None:
                continue
            name = arr.GetName() or "(unnamed)"
            nc = arr.GetNumberOfComponents()
            dtype = arr.GetDataTypeAsString()
            comp_label = f"{nc} comp" + ("s" if nc != 1 else " ")
            rows.append(f"{indent}{name:<28}{comp_label:<10}{dtype}")
        return rows

    def _field_lines(fd: Any, indent: str = "    ") -> list[str]:
        rows: list[str] = []
        for i in range(fd.GetNumberOfArrays()):
            arr = fd.GetAbstractArray(i)
            if arr is None:
                continue
            name = arr.GetName() or "(unnamed)"
            nc = arr.GetNumberOfComponents()
            dtype = "unknown"
            if hasattr(arr, "GetDataTypeAsString"):
                try:
                    dtype = arr.GetDataTypeAsString()
                except Exception:
                    dtype = "unknown"
            comp_label = f"{nc} comp" + ("s" if nc != 1 else " ")
            rows.append(f"{indent}{name:<28}{comp_label:<10}{dtype}")
        return rows

    pd_obj = mesh.GetPointData()
    cd_obj = mesh.GetCellData()
    fd_obj = mesh.GetFieldData()

    lines: list[str] = []
    lines.append("")
    lines.append(sep)
    lines.append(f"  {'APEC-M  Mesh Summary':^{w}}")
    lines.append(sep)
    lines.append(f"  {'Job name':<18}: {mesh_name}")
    lines.append(f"  {'Working dir':<18}: {cwd}")
    lines.append(f"  {'Mesh file':<18}: {mesh_path or '(in-memory grid)'}")
    lines.append(sep)
    lines.append("  Geometry")
    lines.append(f"    {'Points':<20}: {n_pts:,}")
    lines.append(f"    {'Cells':<20}: {n_cells:,}")
    lines.append(f"    {'Spatial dimension':<20}: {spatial_dim}")
    lines.append(f"    {'Bounds X':<20}: {_bstr(bounds[0], bounds[1])}")
    lines.append(f"    {'Bounds Y':<20}: {_bstr(bounds[2], bounds[3])}")
    lines.append(f"    {'Bounds Z':<20}: {_bstr(bounds[4], bounds[5])}")
    lines.append(sep)
    lines.append(f"  Point arrays ({pd_obj.GetNumberOfArrays()})")
    lines.extend(_array_lines(pd_obj))
    lines.append(sep)
    lines.append(f"  Cell arrays ({cd_obj.GetNumberOfArrays()})")
    lines.extend(_array_lines(cd_obj))
    lines.append(sep)
    lines.append(f"  Field arrays ({fd_obj.GetNumberOfArrays()})")
    lines.extend(_field_lines(fd_obj))
    lines.append(sep)
    return "\n".join(lines)
