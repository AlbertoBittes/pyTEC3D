import sys
import os

# Ensure the package root (<repo>/python) is on sys.path when run as a script.
PACKAGE_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PACKAGE_ROOT not in sys.path:
    sys.path.insert(0, PACKAGE_ROOT)

from apec_m import APEC_M

apec = APEC_M()
