"""
Tests for the custom-analysis loader (``ANALYSIS_MODULE`` / ``ANA_CUSTOM``).

These exercise :func:`pimms.file_utilities.custom_analysis_module_import`, which
loads and validates a user-supplied ``analysis_function(step, lattice)`` and is the
entry point for runtime custom analysis code.
"""

import textwrap
import sys

import pytest

from pimms import file_utilities
from pimms.latticeExceptions import KeyFileException


def _write_module(tmp_path, name, body):
    """Write a python source file into tmp_path and return its path (as str)."""
    path = tmp_path / name
    path.write_text(textwrap.dedent(body))
    return str(path)


# --------------------------------------------------------------------------- #
# happy path
# --------------------------------------------------------------------------- #
def test_valid_module_returns_callable(tmp_path):
    path = _write_module(
        tmp_path,
        "good_analysis.py",
        """
        def analysis_function(step, lattice):
            return (step, lattice)
        """,
    )
    fn = file_utilities.custom_analysis_module_import(path)
    assert callable(fn)
    assert fn(7, "LATTICE") == (7, "LATTICE")


def test_valid_module_with_extra_defaulted_arg(tmp_path):
    # extra parameters are allowed provided they have defaults (PIMMS passes two)
    path = _write_module(
        tmp_path,
        "extra_arg.py",
        """
        def analysis_function(step, lattice, verbose=False):
            return step
        """,
    )
    fn = file_utilities.custom_analysis_module_import(path)
    assert fn(3, None) == 3


def test_valid_module_with_var_args(tmp_path):
    path = _write_module(
        tmp_path,
        "var_args.py",
        """
        def analysis_function(*args, **kwargs):
            return len(args)
        """,
    )
    fn = file_utilities.custom_analysis_module_import(path)
    assert fn(1, None) == 2


def test_module_can_import_sibling_helper(tmp_path):
    _write_module(tmp_path, "helper_mod.py", "VALUE = 42\n")
    path = _write_module(
        tmp_path,
        "uses_helper.py",
        """
        import helper_mod
        def analysis_function(step, lattice):
            return helper_mod.VALUE
        """,
    )
    fn = file_utilities.custom_analysis_module_import(path)
    assert fn(0, None) == 42


def test_module_can_lazily_import_sibling_helper_without_polluting_path(tmp_path):
    _write_module(tmp_path, "lazy_helper_mod.py", "VALUE = 43\n")
    path = _write_module(
        tmp_path,
        "uses_lazy_helper.py",
        """
        def analysis_function(step, lattice):
            import lazy_helper_mod
            return lazy_helper_mod.VALUE
        """,
    )

    fn = file_utilities.custom_analysis_module_import(path)

    assert str(tmp_path) not in sys.path
    assert fn(0, None) == 43
    assert str(tmp_path) not in sys.path


def test_sibling_helper_wins_over_same_named_module_earlier_on_path(tmp_path, monkeypatch):
    """A local helper must not be silently replaced by an unrelated module."""
    analysis_dir = tmp_path / "analysis"
    other_dir = tmp_path / "other"
    analysis_dir.mkdir()
    other_dir.mkdir()
    helper_name = "pimms_rereview_shadow_helper"
    _write_module(other_dir, f"{helper_name}.py", "VALUE = 'wrong'\n")
    _write_module(analysis_dir, f"{helper_name}.py", "VALUE = 'local'\n")
    path = _write_module(
        analysis_dir,
        "uses_shadowed_helper.py",
        f"""
        import {helper_name}
        def analysis_function(step, lattice):
            return {helper_name}.VALUE
        """,
    )
    monkeypatch.syspath_prepend(str(other_dir))
    sys.modules.pop(helper_name, None)

    fn = file_utilities.custom_analysis_module_import(path)

    assert fn(0, None) == "local"
    assert str(analysis_dir) not in sys.path


def test_sanitised_path_collisions_still_get_distinct_module_names(tmp_path):
    dash_dir = tmp_path / "a-b"
    underscore_dir = tmp_path / "a_b"
    dash_dir.mkdir()
    underscore_dir.mkdir()
    p1 = _write_module(
        dash_dir, "analysis.py",
        "def analysis_function(step, lattice):\n    return __name__\n",
    )
    p2 = _write_module(
        underscore_dir, "analysis.py",
        "def analysis_function(step, lattice):\n    return __name__\n",
    )

    assert file_utilities.custom_analysis_module_import(p1)(0, None) != (
        file_utilities.custom_analysis_module_import(p2)(0, None)
    )


def test_same_basename_different_dirs_do_not_collide(tmp_path):
    # two different files both called analysis.py must load independently
    d1 = tmp_path / "a"
    d2 = tmp_path / "b"
    d1.mkdir()
    d2.mkdir()
    p1 = _write_module(d1, "analysis.py",
                       "def analysis_function(step, lattice):\n    return 'A'\n")
    p2 = _write_module(d2, "analysis.py",
                       "def analysis_function(step, lattice):\n    return 'B'\n")
    fn1 = file_utilities.custom_analysis_module_import(p1)
    fn2 = file_utilities.custom_analysis_module_import(p2)
    assert fn1(0, None) == "A"
    assert fn2(0, None) == "B"


# --------------------------------------------------------------------------- #
# failure modes - each must raise a clear KeyFileException
# --------------------------------------------------------------------------- #
def test_missing_file_raises(tmp_path):
    with pytest.raises(KeyFileException, match="Could not find"):
        file_utilities.custom_analysis_module_import(str(tmp_path / "nope.py"))


@pytest.mark.parametrize("bad", [None, "", "   ", 123])
def test_non_string_or_empty_path_raises(bad):
    with pytest.raises(KeyFileException, match="ANALYSIS_MODULE must be a path"):
        file_utilities.custom_analysis_module_import(bad)


def test_syntax_error_raises_clear_error(tmp_path):
    path = _write_module(
        tmp_path,
        "broken_syntax.py",
        """
        def analysis_function(step, lattice)   # missing colon
            return step
        """,
    )
    with pytest.raises(KeyFileException, match="raised .* while being loaded"):
        file_utilities.custom_analysis_module_import(path)


def test_error_at_import_time_raises_clear_error(tmp_path):
    path = _write_module(
        tmp_path,
        "raises_on_import.py",
        """
        raise ValueError("boom at import")
        def analysis_function(step, lattice):
            return step
        """,
    )
    with pytest.raises(KeyFileException, match="Failed to import"):
        file_utilities.custom_analysis_module_import(path)


def test_system_exit_at_import_is_reported_as_a_keyfile_error(tmp_path):
    path = _write_module(
        tmp_path,
        "exits_on_import.py",
        """
        raise SystemExit(7)
        def analysis_function(step, lattice):
            return step
        """,
    )
    with pytest.raises(KeyFileException, match="raised SystemExit"):
        file_utilities.custom_analysis_module_import(path)


def test_missing_analysis_function_raises(tmp_path):
    path = _write_module(
        tmp_path,
        "no_entry.py",
        """
        def some_other_function(step, lattice):
            return step
        """,
    )
    with pytest.raises(KeyFileException, match="does not define an"):
        file_utilities.custom_analysis_module_import(path)


def test_analysis_function_not_callable_raises(tmp_path):
    path = _write_module(
        tmp_path,
        "not_callable.py",
        """
        analysis_function = 42
        """,
    )
    with pytest.raises(KeyFileException, match="is not callable"):
        file_utilities.custom_analysis_module_import(path)


@pytest.mark.parametrize(
    "sig",
    [
        "def analysis_function(lattice):",            # too few
        "def analysis_function():",                   # none
        "def analysis_function(a, b, c):",            # too many required
    ],
)
def test_wrong_signature_raises(tmp_path, sig):
    path = _write_module(
        tmp_path,
        "bad_sig.py",
        f"""
        {sig}
            return None
        """,
    )
    with pytest.raises(KeyFileException, match="analysis_function.* has signature"):
        file_utilities.custom_analysis_module_import(path)
