#!/bin/zsh
#
# Clean rebuild + editable install of PIMMS.
#
#     ./build.sh uv      reinstall with uv:   uv pip install -e . --no-deps --reinstall
#     ./build.sh pip     reinstall with pip:  python -m pip install -e . --force-reinstall --no-deps
#
# Cython's cythonize() SKIPS regenerating a .c file when that .c is newer than its .pyx,
# and build_ext can reuse cached .o files under build/. So a plain reinstall might NOT
# pick up .pyx changes. Removing the generated C (pimms/*.c and the lemonade kernel),
# the compiled extensions (*.so), and the build/ object cache guarantees that EVERY
# .pyx is recompiled from scratch.

# stop on the first real error (e.g. a failed compile) so it is visible
set -e

usage() {
    echo "usage: ./build.sh {uv|pip}" >&2
    echo "  uv    clean rebuild, then 'uv pip install -e . --no-deps --reinstall'" >&2
    echo "  pip   clean rebuild, then 'python -m pip install -e . --force-reinstall --no-deps'" >&2
    exit 2
}

# exactly one argument naming the installer; check the tool works BEFORE deleting
# anything, so a typo or a missing tool never leaves the tree without its extensions.
if [[ $# -ne 1 ]]; then
    usage
fi
installer=$1
case $installer in
    uv)
        if ! command -v uv > /dev/null 2>&1; then
            echo "build.sh: 'uv' is not on PATH - install it, or run './build.sh pip'" >&2
            exit 1
        fi
        ;;
    pip)
        # a uv-managed venv ships without pip, so 'python -m pip' can be missing
        # even though python itself works
        if ! python -m pip --version > /dev/null 2>&1; then
            echo "build.sh: 'python -m pip' does not work in this environment - run './build.sh uv', or install pip into it" >&2
            exit 1
        fi
        ;;
    *)
        usage
        ;;
esac

# the sweep and the install both assume the repo root, whatever directory the
# script was launched from
cd "${0:A:h}"

# remove generated C source, compiled extensions, and cached object files.
# (find tolerates the "no matches" case, unlike a bare shell glob under set -e.)
find pimms pimms/lemonade/kernels -maxdepth 1 \( -name '*.so' -o -name '*.c' \) -delete
rm -rf build/

# editable install, rebuilding ONLY PIMMS (--no-deps leaves other env packages
# untouched; the reinstall flag forces a fresh build - and, for uv, refreshes its
# build cache - rather than reusing what the installer thinks is already there).
case $installer in
    uv)
        uv pip install -e . --no-deps --reinstall
        ;;
    pip)
        python -m pip install -e . --force-reinstall --no-deps
        ;;
esac
