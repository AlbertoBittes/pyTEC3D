"""The version of pyTEC3D: VERSION (repository root) is the single source; this script copies it where the build tools
need it as a literal.

    python tools/set_version.py 0.4.0     writes 0.4.0 to VERSION and to every copy
    python tools/set_version.py           copies VERSION to every copy
    python tools/set_version.py --check   exit 1 if a copy differs from VERSION (tests/test_version.py, the CI);
                                          with GITHUB_REF=refs/tags/v<X> (a tag build) also when X differs

The copies: the version of the three distributions (python/, element_library/, solver/pyproject.toml: setuptools
cannot read a file outside the folder of its pyproject), the pin of pytec3d-element-library in python/pyproject.toml,
the context version of conda/recipe.yaml. solver/CMakeLists.txt reads solver/pyproject.toml, docs/conf.py reads VERSION.
"""
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VERSION_FILE = os.path.join(ROOT, "VERSION")
NUMBER = r"[0-9]+\.[0-9]+\.[0-9]+"
# file: (pattern whose group 2 is the version, how many matches)
COPIES = {
    "python/pyproject.toml": [(r'^(version = ")(' + NUMBER + ')(")', 1),
                              (r'^(\s*"pytec3d-element-library==)(' + NUMBER + ')(")', 1)],
    "element_library/pyproject.toml": [(r'^(version = ")(' + NUMBER + ')(")', 1)],
    "solver/pyproject.toml": [(r'^(version = ")(' + NUMBER + ')(")', 1)],
    "conda/recipe.yaml": [(r'^(  version: ")(' + NUMBER + ')(")', 1)],
}


def read_version():
    with open(VERSION_FILE) as f:
        version = f.read().strip()
    if not re.fullmatch(NUMBER, version):
        sys.exit(f"VERSION: '{version}' is not X.Y.Z")
    return version


def copies():
    """[(file, line number, version found)] of every copy."""
    found = []
    for rel, patterns in COPIES.items():
        with open(os.path.join(ROOT, rel)) as f:
            text = f.read()
        for pattern, count in patterns:
            matches = list(re.finditer(pattern, text, re.MULTILINE))
            if len(matches) != count:
                sys.exit(f"{rel}: {len(matches)} match(es) of {pattern}, expected {count}")
            found += [(rel, text.count("\n", 0, m.start()) + 1, m.group(2)) for m in matches]
    return found


def write(version):
    for rel, patterns in COPIES.items():
        path = os.path.join(ROOT, rel)
        with open(path) as f:
            text = f.read()
        for pattern, _ in patterns:
            text = re.sub(pattern, lambda m: m.group(1) + version + m.group(3), text, flags=re.MULTILINE)
        with open(path, "w") as f:
            f.write(text)


def check(version):
    problems = [f"{rel}:{line}: {v}" for rel, line, v in copies() if v != version]
    ref = os.environ.get("GITHUB_REF", "")
    if ref.startswith("refs/tags/v") and ref[len("refs/tags/v"):] != version:
        problems.append(f"tag {ref[len('refs/tags/'):]}")
    for p in problems:
        print(f"version {version} (VERSION) differs: {p}")
    return 1 if problems else 0


def main(argv):
    if argv == ["--check"]:
        return check(read_version())
    if len(argv) > 1:
        sys.exit(__doc__)
    if argv:
        if not re.fullmatch(NUMBER, argv[0]):
            sys.exit(f"'{argv[0]}' is not X.Y.Z")
        with open(VERSION_FILE, "w") as f:
            f.write(argv[0] + "\n")
    version = read_version()
    write(version)
    print(f"version {version}: " + ", ".join(sorted(COPIES)))
    return check(version)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
