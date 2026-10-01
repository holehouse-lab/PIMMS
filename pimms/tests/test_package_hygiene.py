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

import importlib.util
import pathlib
import re
import subprocess
import sys
import sysconfig
import fnmatch
import shlex
import shutil

import pytest

_PACKAGE_DIR = pathlib.Path(__file__).resolve().parents[1]
_REPO_ROOT = _PACKAGE_DIR.parent

# directories under pimms/ whose modules are NOT part of the shipped import surface
_EXCLUDED_PARTS = {"tests", "fast_kernels", "cython_backup", ".pytest_cache", "__pycache__"}

_REGRESSION_CONFTEST = _PACKAGE_DIR / "tests" / "simulation_tests" / "conftest.py"


def _require_source_checkout() -> pathlib.Path:
    """Return the repository root, or skip the calling test from an installed copy.

    The tests ship in the wheel, where there is no ``MANIFEST.in``, ``pyproject.toml``
    or ``docs/`` beside the package. A test that reads those files can only run from
    a source checkout (or an unpacked sdist), so from an install it is skipped with
    a reason rather than failed for a file that was never meant to be there.

    Returns
    -------
    pathlib.Path
        The directory holding ``pyproject.toml``, ``MANIFEST.in`` and ``pimms/``.
    """
    if not ((_REPO_ROOT / "MANIFEST.in").is_file() and (_REPO_ROOT / "pyproject.toml").is_file()):
        pytest.skip(
            f"needs the packaging files of a source checkout; {_REPO_ROOT} is an installed copy"
        )
    return _REPO_ROOT


def _manifest_commands(root: pathlib.Path) -> dict[str, list[str]]:
    """Read ``MANIFEST.in`` into a mapping of command to the patterns it is given.

    Parameters
    ----------
    root : pathlib.Path
        Directory holding ``MANIFEST.in``.

    Returns
    -------
    dict[str, list[str]]
        For each manifest command (``include``, ``graft``, ``global-exclude``, ...),
        every pattern passed to it anywhere in the file, in file order.
    """
    commands: dict[str, list[str]] = {}
    for raw_line in (root / "MANIFEST.in").read_text().splitlines():
        tokens = shlex.split(raw_line, comments=True)
        if tokens:
            commands.setdefault(tokens[0], []).extend(tokens[1:])
    return commands


def _pyproject(root: pathlib.Path) -> dict:
    """Parse ``pyproject.toml``, skipping the calling test if no TOML reader exists.

    ``tomllib`` is in the standard library from Python 3.11. On 3.10 we use ``tomli``
    when it happens to be installed and skip otherwise - it is not a PIMMS dependency.

    Parameters
    ----------
    root : pathlib.Path
        Directory holding ``pyproject.toml``.

    Returns
    -------
    dict
        The parsed file.
    """
    try:
        import tomllib
    except ImportError:                                  # Python 3.10
        tomllib = pytest.importorskip("tomli")
    return tomllib.loads((root / "pyproject.toml").read_text())


def _floor(requirements: list[str], name: str) -> tuple[int, ...]:
    """Return the ``>=`` floor a requirement list sets for one distribution.

    Parameters
    ----------
    requirements : list[str]
        Requirement strings as written in ``pyproject.toml`` (``"numpy>=1.23"``).
    name : str
        Distribution name, compared without regard to case.

    Returns
    -------
    tuple[int, ...]
        The floor as a tuple of integers, ``()`` if the requirement has no ``>=``.

    Raises
    ------
    AssertionError
        If ``name`` is not in the list exactly once.
    """
    matches = [r for r in requirements
               if re.match(rf"\s*{re.escape(name)}\s*($|[<>=!~;\[])", r, flags=re.IGNORECASE)]
    assert len(matches) == 1, f"expected exactly one requirement for {name}, got {matches}"
    found = re.search(r">=\s*([0-9]+(?:\.[0-9]+)*)", matches[0])
    return tuple(int(part) for part in found.group(1).split(".")) if found else ()


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
    manifest = _manifest_commands(_require_source_checkout())
    global_excludes = manifest.get("global-exclude", [])
    includes = manifest.get("include", [])

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


def test_manifest_excludes_generated_c_and_the_sdist_can_regenerate_it():
    """The Cython-generated C must not ship, and leaving it out must be safe.

    ``setup.py`` runs ``cythonize`` on import, which writes a ``.c`` next to every
    ``.pyx`` and makes it the extension's source, so ``graft pimms`` swept ten
    generated files (16 MB of a 23 MB wheel) into both the sdist and the wheel. In an
    sdist they also decide, by file timestamp, whether the build environment's Cython
    is used at all. Excluding them is only right if an install from the sdist can
    regenerate them, so this pins the three things that makes true: the exclusion,
    the ``.pyx`` sources still being packaged, and Cython 3 being a build requirement.
    """
    root = _require_source_checkout()
    manifest = _manifest_commands(root)

    # every extension source in the tree, found from the tree rather than from a list
    pyx_sources = sorted(p.relative_to(root).as_posix() for p in (root / "pimms").rglob("*.pyx"))
    assert "pimms/mega_crank_fast.pyx" in pyx_sources          # sanity: the walk found them
    assert "pimms/lemonade/kernels/_pbc.pyx" in pyx_sources

    for pyx in pyx_sources:
        generated = pyx[:-len(".pyx")] + ".c"
        assert any(fnmatch.fnmatch(pathlib.PurePosixPath(generated).name, pattern)
                   for pattern in manifest.get("global-exclude", [])), (
            f"MANIFEST.in would package the Cython-generated {generated}"
        )
        # the .pyx itself must survive: grafted, and matched by no exclusion
        assert any(pyx.startswith(tree + "/") for tree in manifest.get("graft", []))
        assert not any(fnmatch.fnmatch(pathlib.PurePosixPath(pyx).name, pattern)
                       for pattern in manifest.get("global-exclude", []))
        assert not any(pyx.startswith(tree + "/") for tree in manifest.get("prune", []))

    build_requires = _pyproject(root)["build-system"]["requires"]
    assert _floor(build_requires, "Cython") >= (3, 0), (
        "the sdist ships no generated C, so Cython 3 has to be a declared build "
        f"requirement; build-system.requires is {build_requires}"
    )


def test_dependency_floors_are_the_demonstrated_ones():
    """The runtime floors must be versions that can be installed and that work.

    Two floors were wrong. ``mdtraj>=1.10`` admitted 1.10.0, whose PDB reader raises
    ``KeyError`` on an atom serial that has wrapped past 99,999, so any system of
    100,000 beads or more could not be read back. ``numpy>=1.21`` could never be
    selected, because mdtraj 1.10.x itself requires numpy >= 1.23.
    """
    root = _require_source_checkout()
    dependencies = _pyproject(root)["project"]["dependencies"]

    assert _floor(dependencies, "mdtraj") >= (1, 10, 1)
    assert _floor(dependencies, "numpy") >= (1, 23)
    assert _floor(dependencies, "scipy") >= (1, 9)

    # the installation page must state the floors pyproject.toml declares
    page = (root / "docs" / "installation.rst").read_text()
    for name in ("numpy", "scipy", "mdtraj"):
        floor = ".".join(str(part) for part in _floor(dependencies, name))
        assert re.search(rf"``{name}``\s*\(≥ {re.escape(floor)}[;)]", page), (
            f"docs/installation.rst does not state {name} ≥ {floor}"
        )


def test_git_archive_build_reports_the_real_version(tmp_path):
    """A tree made by ``git archive`` must report its version, not ``0+unknown``.

    A GitHub tarball has no ``.git`` directory, so versioningit's default method
    found nothing to describe and every such build - and so ``PIMMS --version``,
    ``keyfile_used.kf`` and the restart file's ``PIMMS_VERSION`` - said ``0+unknown``.
    The fix has two halves that only work together: ``pyproject.toml`` carries a
    describe placeholder, and ``.gitattributes`` marks that file ``export-subst`` so
    ``git archive`` fills the placeholder in. Here we write the file ``git archive``
    would produce and ask versioningit what version it reads from it.
    """
    root = _require_source_checkout()
    text = (root / "pyproject.toml").read_text()

    # `tags` is required: the release tags are lightweight, and a plain describe
    # only considers annotated ones. git expands at most one describe placeholder
    # per archive, so there must be exactly one.
    placeholder = "$Format:%(describe:tags)$"
    assert text.count(placeholder) == 1
    vcs = _pyproject(root)["tool"]["versioningit"]["vcs"]
    assert vcs == {"method": "git-archive", "describe-subst": placeholder}

    attributes = [line.split() for line in (root / ".gitattributes").read_text().splitlines()
                  if line.strip() and not line.lstrip().startswith("#")]
    assert ["pyproject.toml", "export-subst"] in attributes

    versioningit = pytest.importorskip("versioningit")
    archive = tmp_path / "archive"
    archive.mkdir()
    # three commits past a tag, then exactly on a tag: the two strings `git describe
    # --tags` gives, and the versions the default format makes of them
    for described, expected in (("v9.8.7-3-gabc1234", "9.8.7.post3+gabc1234"), ("v9.8.7", "9.8.7")):
        (archive / "pyproject.toml").write_text(text.replace(placeholder, described))
        assert versioningit.get_version(archive) == expected


def _load_regression_conftest(conftest_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch):
    """Import a copy of the regression-suite conftest from an arbitrary location.

    Parameters
    ----------
    conftest_path : pathlib.Path
        Where the copy has been written. The helpers under test work out where they
        are from their own ``__file__``, so the location is the input of the test.
    monkeypatch : pytest.MonkeyPatch
        Used to register the module for the duration of the test (``dataclasses``
        needs the defining module in ``sys.modules``).

    Returns
    -------
    module
        The imported copy.
    """
    name = "_pimms_regression_conftest_copy"
    spec = importlib.util.spec_from_file_location(name, conftest_path)
    module = importlib.util.module_from_spec(spec)
    monkeypatch.setitem(sys.modules, name, module)
    spec.loader.exec_module(module)
    return module


def _fake_tree(root: pathlib.Path, *files: str) -> None:
    """Create empty files (and their parent directories) under ``root``.

    Parameters
    ----------
    root : pathlib.Path
        Directory to create the files in.
    *files : str
        Paths relative to ``root``.

    Returns
    -------
    None
    """
    for relative in files:
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("")


def test_regression_harness_does_not_adopt_an_enclosing_checkout(tmp_path, monkeypatch):
    """An installed copy of the regression suite must test the installed PIMMS.

    The harness used to walk up from its own location until it met a directory with
    ``setup.py`` and ``scripts/PIMMS``. For a wheel installed into an environment
    that lives inside a checkout, the walk left ``site-packages``, found the
    enclosing checkout, and the 15 regression scenarios then ran that checkout's
    script against that checkout's ``pimms`` - reporting on a build nobody had asked
    about. We rebuild that layout and check which executable and which import root
    the harness chooses.
    """
    checkout = tmp_path / "checkout"
    _fake_tree(checkout, "setup.py", "scripts/PIMMS", "env/bin/PIMMS")
    site_packages = checkout / "env" / "lib" / "site-packages"
    installed = site_packages / "pimms" / "tests" / "simulation_tests" / "conftest.py"
    installed.parent.mkdir(parents=True)
    shutil.copyfile(_REGRESSION_CONFTEST, installed)

    real_get_path = sysconfig.get_path
    monkeypatch.setattr(
        sysconfig, "get_path",
        lambda name, *a, **k: str(checkout / "env" / "bin") if name == "scripts" else real_get_path(name, *a, **k),
    )
    harness = _load_regression_conftest(installed, monkeypatch)

    command = harness._resolve_pimms_command()
    assert command == [sys.executable, str((checkout / "env" / "bin" / "PIMMS").resolve())], (
        f"an installed regression suite would run {command[-1]}"
    )
    assert harness._import_root() == site_packages.resolve()
    assert harness._source_tree_root() is None


def test_regression_harness_runs_the_working_tree_from_a_checkout(tmp_path, monkeypatch):
    """From a source checkout the harness must keep running ``scripts/PIMMS``.

    This is the behaviour the walk-up was written for, and it must survive the fix:
    a ``PIMMS`` installed in the environment is a different build from the tree
    under test, so it must not be picked while a checkout's own script is there.
    """
    checkout = tmp_path / "checkout"
    _fake_tree(checkout, "setup.py", "scripts/PIMMS", "elsewhere/bin/PIMMS")
    in_tree = checkout / "pimms" / "tests" / "simulation_tests" / "conftest.py"
    in_tree.parent.mkdir(parents=True)
    shutil.copyfile(_REGRESSION_CONFTEST, in_tree)

    real_get_path = sysconfig.get_path
    monkeypatch.setattr(
        sysconfig, "get_path",
        lambda name, *a, **k: str(checkout / "elsewhere" / "bin") if name == "scripts" else real_get_path(name, *a, **k),
    )
    harness = _load_regression_conftest(in_tree, monkeypatch)

    assert harness._resolve_pimms_command() == [sys.executable, str((checkout / "scripts" / "PIMMS").resolve())]
    assert harness._import_root() == checkout.resolve()


def test_regression_harness_says_so_when_an_install_has_no_executable(tmp_path, monkeypatch):
    """An installed suite that cannot find ``PIMMS`` must fail with the places it looked.

    A regression suite that cannot run must not pass, and must not fall back to some
    other tree; the message has to tell the user what to fix.
    """
    site_packages = tmp_path / "env" / "lib" / "site-packages"
    installed = site_packages / "pimms" / "tests" / "simulation_tests" / "conftest.py"
    installed.parent.mkdir(parents=True)
    shutil.copyfile(_REGRESSION_CONFTEST, installed)
    empty_bin = tmp_path / "env" / "bin"
    empty_bin.mkdir()

    real_get_path = sysconfig.get_path
    monkeypatch.setattr(
        sysconfig, "get_path",
        lambda name, *a, **k: str(empty_bin) if name == "scripts" else real_get_path(name, *a, **k),
    )
    monkeypatch.setattr(sys, "executable", str(empty_bin / "python"))
    monkeypatch.setattr(shutil, "which", lambda *a, **k: None)
    harness = _load_regression_conftest(installed, monkeypatch)

    with pytest.raises(RuntimeError) as excinfo:
        harness._resolve_pimms_command()
    message = str(excinfo.value)
    assert str(empty_bin / "PIMMS") in message and "PATH" in message
    assert str(site_packages.resolve()) in message


def _fake_regression_suite(tmp_path: pathlib.Path) -> tuple[pathlib.Path, pathlib.Path, dict[str, bytes]]:
    """Build a throwaway checkout holding the regression harness and one fixture.

    The checkout has a stand-in ``pimms`` package (so the harness's "which pimms
    would the run import" probe has something to find), a stand-in ``scripts/PIMMS``
    that writes the files a real run writes, and a ``test_1`` fixture holding both
    input files and outputs left behind by an earlier run.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Directory to build the checkout in.

    Returns
    -------
    tuple[pathlib.Path, pathlib.Path, dict[str, bytes]]
        The path of the harness copy, the fixture directory, and the name and
        content of every file placed in the fixture directory.
    """
    checkout = tmp_path / "checkout"
    _fake_tree(checkout, "setup.py", "pimms/__init__.py", "pimms/tests/__init__.py")
    script = checkout / "scripts" / "PIMMS"
    script.parent.mkdir()
    script.write_text(
        "import pathlib\n"
        "assert pathlib.Path('KEYFILE.kf').is_file() and pathlib.Path('params.prm').is_file()\n"
        "pathlib.Path('ENERGY.dat').write_text('10\\t-1.0\\n20\\t-2.0\\n')\n"
        "pathlib.Path('absolute_energies_of_angles.txt').write_text('fresh\\n')\n"
        "pathlib.Path('restart.pimms').write_text('fresh\\n')\n"
        "print('stand-in run complete')\n"
    )
    suite = checkout / "pimms" / "tests" / "simulation_tests"
    fixture = suite / "test_1"
    fixture.mkdir(parents=True)
    harness_path = suite / "conftest.py"
    shutil.copyfile(_REGRESSION_CONFTEST, harness_path)
    (suite / "expected_output").mkdir()
    (suite / "expected_output" / "ENERGY.dat.final_lines.txt").write_text("test_1\t20\t-2.0\n")

    contents = {
        # inputs
        "KEYFILE.kf": b"input\n", "params.prm": b"input\n", "in.pimms": b"input\n", "frz.in": b"input\n",
        # outputs of an earlier run; the angle table stands for the one that is tracked
        "absolute_energies_of_angles.txt": b"tracked\n", "ENERGY.dat": b"5\t-9.0\n", "START.pdb": b"old\n",
        "traj.xtc": b"old\n", "restart.pimms": b"old\n", "log.txt": b"old\n",
    }
    for name, data in contents.items():
        (fixture / name).write_bytes(data)
    return harness_path, fixture, contents


def test_regression_harness_runs_outside_the_fixture_directories(tmp_path, monkeypatch):
    """A regression run must read its fixture directory and write somewhere else.

    The harness used to delete the previous outputs from the fixture directory -
    which took the git-tracked ``test_11/absolute_energies_of_angles.txt`` with them
    - and then run PIMMS there. Every test run dirtied the tree, and two runs on
    one checkout removed each other's files part-way through. We run the real
    harness on a stand-in suite and compare the fixture directory before and after.
    """
    harness_path, fixture, contents = _fake_regression_suite(tmp_path)
    harness = _load_regression_conftest(harness_path, monkeypatch)
    run_dir = tmp_path / "run"
    run_dir.mkdir()

    returned_dir, final_lines, line_counts = harness._run_single_testset(1, run_dir)

    # the fixture directory is byte-for-byte what it was, and nothing was added
    assert {p.name: p.read_bytes() for p in fixture.iterdir()} == contents

    # the run happened in run_dir: inputs copied in, stale outputs not, new outputs there
    assert returned_dir == run_dir
    written = {p.name: p.read_bytes() for p in run_dir.iterdir()}
    for name in ("KEYFILE.kf", "params.prm", "in.pimms", "frz.in"):
        assert written[name] == contents[name]
    assert "START.pdb" not in written and "traj.xtc" not in written and "log.txt" not in written
    assert written["restart.pimms"] == b"fresh\n"
    assert written["absolute_energies_of_angles.txt"] == b"fresh\n"
    assert b"stand-in run complete" in written["pytest_test_1_log.txt"]

    # and the values compared against the baseline are the new run's, not the stale file's
    assert final_lines == {"ENERGY.dat": "20\t-2.0"}
    assert line_counts == {"ENERGY.dat": 2}


def test_regression_harness_reports_the_log_of_a_failed_run(tmp_path, monkeypatch):
    """When a scenario fails, the message must say where its log is, and the log must be there."""
    harness_path, fixture, contents = _fake_regression_suite(tmp_path)
    script = harness_path.parents[3] / "scripts" / "PIMMS"
    script.write_text("import sys\nsys.stderr.write('stand-in failure\\n')\nsys.exit(3)\n")
    harness = _load_regression_conftest(harness_path, monkeypatch)
    run_dir = tmp_path / "run"
    run_dir.mkdir()

    with pytest.raises(AssertionError) as excinfo:
        harness._run_single_testset(1, run_dir)

    log_path = run_dir / "pytest_test_1_log.txt"
    assert "return code 3" in str(excinfo.value) and str(log_path) in str(excinfo.value)
    assert "stand-in failure" in log_path.read_text()
    assert {p.name: p.read_bytes() for p in fixture.iterdir()} == contents
