"""Pytest configuration to avoid __pycache__ artifacts during test runs."""

from pathlib import Path
import shutil
import sys

# Prevent Python from writing .pyc files while tests run.
sys.dont_write_bytecode = True


def pytest_configure(_config):
    """Remove any existing __pycache__ directories before tests execute."""
    repo_root = Path(__file__).resolve().parent.parent
    for cache_dir in repo_root.rglob("__pycache__"):
        shutil.rmtree(cache_dir, ignore_errors=True)
