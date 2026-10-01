"""
Shared fixtures for the PIMMS test suite.
"""

import os

import pytest


@pytest.fixture(autouse=True)
def _restore_working_directory():
    """Put the working directory back after every test.

    A PIMMS run writes its outputs into the working directory, so many tests
    ``os.chdir`` into a temporary directory first. A test that does not change
    back leaves every later test running inside that directory, where any
    files it left behind (``restart.pimms``, ``*.dat``) are visible to code
    that resolves a relative path against the working directory. That made a
    lemonade test that expects a missing restart file to raise pass or fail
    depending on which tests happened to run before it.

    Yields
    ------
    None
        Control passes to the test; the original directory is restored after.
    """
    cwd = os.getcwd()
    try:
        yield
    finally:
        os.chdir(cwd)
