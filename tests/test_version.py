"""The version copies of the build files agree with VERSION (tools/set_version.py)."""
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_version_copies():
    env = dict(os.environ)
    env.pop("GITHUB_REF", None)
    result = subprocess.run([sys.executable, os.path.join(ROOT, "tools", "set_version.py"), "--check"],
                            capture_output=True, text=True, env=env)
    assert result.returncode == 0, result.stdout + result.stderr
