"""
Unit and regression test for the pimms package.
"""

# Import package, test suite, and other packages as needed
import pimms
import pytest
import sys

def test_pimms_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "pimms" in sys.modules
