# devtools

These files came with the MolSSI computational-molecular-science cookiecutter template this repository was created from in 2019, and supported a Travis-CI setup that has since been removed. Nothing in the repository uses them any more - no CI configuration, build file or script refers to this directory - and they are kept for reference only.

## Contents

* `travis-ci/before_install.sh` - installs Miniconda on a Travis-CI worker (the old `before_install` step).
* `scripts/create_conda_env.py` - creates a conda environment from a YAML environment file with a given name and Python version (`python create_conda_env.py -n NAME -p PYTHON_VERSION ENV_FILE`).
* `conda-envs/test_env.yaml` - a minimal conda test environment (`python`, `pip`, `pytest`, `pytest-cov`, `codecov`). It does not list PIMMS' own dependencies, so on its own it cannot run the test suite.

## Building and testing PIMMS

See the installation page of the documentation (`docs/installation.rst`) for how to build PIMMS, compile the Cython kernels and run the test suite with `pytest`.
