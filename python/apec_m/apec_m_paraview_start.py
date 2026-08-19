import sys
import os

from apec_m.logging_utils import get_logger

# Ensure the package root (<repo>/python) is on sys.path when run as a script.
PACKAGE_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PACKAGE_ROOT not in sys.path:
    sys.path.insert(0, PACKAGE_ROOT)

from apec_m import APEC_M


LOGGER = get_logger("apec_m.paraview_start")


def start_apec_for_paraview() -> APEC_M:
    """Create and return an APEC_M instance for ParaView script startup.

    Returns:
        APEC_M: Initialized package entry object.
    """
    LOGGER.info("Initializing APEC_M in ParaView.")
    return APEC_M()


apec = start_apec_for_paraview()
