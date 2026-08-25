"""
Package-hygiene tests: importing any shipped pimms module must be side-effect free,
and no file outside the test packages may look like a test to pytest.

These pin a bug class that has now happened twice. Developer scripts living inside
the package (``test_megagrank.py``, then ``randint_test.py`` / ``randneg_test.py``)
matched pytest's default collection patterns (``test_*.py`` / ``*_test.py``), so
``pytest pimms/`` imported them during collection and executed their top-level
bodies - one of which printed ten million random integers. The same bodies also ran
on a plain ``import``, and every one of these files ships in the wheel as an
importable module. The scripts are now ``dev_*``-named and ``__main__``-guarded;
these tests keep it that way.
"""

import pathlib
import subprocess
import sys
import fnmatch
import shlex

_PACKAGE_DIR = pathlib.Path(__file__).resolve().parents[1]
_REPO_ROOT = _PACKAGE_DIR.parent

# directories under pimms/ whose modules are NOT part of the shipped import surface
_EXCLUDED_PARTS = {"tests", "fast_kernels", "cython_backup", ".pytest_cache", "__pycache__"}


def _shipped_modules():
    """Dotted names of every pure-Python module shipped in the pimms package."""
    modules = []
    for path in sorted(_PACKAGE_DIR.rglob("*.py")):
        rel = path.relative_to(_PACKAGE_DIR)
        if any(part in _EXCLUDED_PARTS for part in rel.parts):
            continue
        parts = ["pimms"] + list(rel.parts)
        parts[-1] = parts[-1][:-3]                      # strip .py
        if parts[-1] == "__init__":
            parts = parts[:-1]
        modules.append(".".join(parts))
    return modules


def test_every_shipped_module_imports_without_side_effects():
    """Importing any pimms module must produce no stdout and must not crash.

    Run in ONE fresh interpreter so each module's import side effects (if any)
    actually fire - in the test process most modules are already imported, which
    would make an in-process check vacuous.
    """
    modules = _shipped_modules()
    assert "pimms.simulation" in modules                # sanity: the walk found the package

    script = (
        "import importlib, io, contextlib, sys\n"
        "failures = []\n"
        f"for name in {modules!r}:\n"
        "    buf = io.StringIO()\n"
        "    try:\n"
        "        with contextlib.redirect_stdout(buf):\n"
        "            importlib.import_module(name)\n"
        "    except Exception as e:\n"
        "        failures.append(f'{name}: raised {type(e).__name__}: {e}')\n"
        "        continue\n"
        "    if buf.getvalue():\n"
        "        failures.append(f'{name}: printed on import: {buf.getvalue()[:120]!r}')\n"
        "for f in failures:\n"
        "    print(f, file=sys.stderr)\n"
        "sys.exit(1 if failures else 0)\n"
    )
    result = subprocess.run([sys.executable, "-c", script],
                            capture_output=True, text=True, timeout=300)
    assert result.returncode == 0, f"import side effects detected:\n{result.stderr}"


def test_no_pytest_collectable_files_outside_the_test_packages():
    """Nothing outside pimms/tests and pimms/lemonade/tests may look like a test.

    pytest's default patterns are ``test_*.py`` and ``*_test.py``; a stray dev
    script matching either is executed at collection time by ``pytest pimms/``.
    (``testpaths`` protects a bare ``pytest`` from the repo root, but not an
    explicit path argument.)
    """
    offenders = []
    for path in _PACKAGE_DIR.rglob("*.py"):
        rel = path.relative_to(_PACKAGE_DIR)
        if "tests" in rel.parts or ".pytest_cache" in rel.parts:
            continue
        name = path.name
        if name.startswith("test_") or name.endswith("_test.py"):
            offenders.append(str(rel))
    assert offenders == [], (
        f"files pytest would collect outside the test packages: {offenders} - "
        "rename them (dev_*.py) so `pytest pimms/` cannot execute them"
    )


def test_manifest_excludes_generated_simulation_outputs():
    """Ignored run artifacts must never be swept into a distribution.

    ``graft pimms`` operates on the filesystem, so the regression suite leaves
    hundreds of otherwise ignored files eligible for packaging unless every
    output basename/suffix is explicitly excluded.
    """
    lines = (_REPO_ROOT / "MANIFEST.in").read_text().splitlines()
    global_excludes = []
    includes = []
    for raw_line in lines:
        tokens = shlex.split(raw_line, comments=True)
        if not tokens:
            continue
        if tokens[0] == "global-exclude":
            global_excludes.extend(tokens[1:])
        elif tokens[0] == "include":
            includes.extend(tokens[1:])

    generated = [
        "ENERGY.dat",
        "traj.xtc",
        "START.pdb",
        "log.txt",
        "restart.pimms",
        "parameters_used.prm",
        "absolute_energies_of_angles.txt",
        "pytest_test_12_log.txt",
    ]
    for filename in generated:
        assert any(fnmatch.fnmatch(filename, pattern) for pattern in global_excludes), (
            f"MANIFEST.in would package generated simulation output {filename}"
        )

    assert "pimms/data/look_and_say.dat" in includes
