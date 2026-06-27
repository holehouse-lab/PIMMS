#!/bin/zsh
#
# Clean rebuild + editable install of PIMMS.
#
# Why the deletes matter: Cython's cythonize() SKIPS regenerating a .c file when
# that .c is newer than its .pyx, and build_ext can reuse cached .o files under
# build/. So a plain reinstall may NOT pick up .pyx changes. Removing the
# generated C (pimms/*.c), the compiled extensions (pimms/*.so), and the build/
# object cache guarantees that EVERY .pyx is recompiled from scratch.

# stop on the first real error (e.g. a failed compile) so it is visible
set -e

# remove generated C source, compiled extensions, and cached object files.
# (find tolerates the "no matches" case, unlike a bare shell glob under set -e.)
find pimms -maxdepth 1 \( -name '*.so' -o -name '*.c' \) -delete
rm -rf build/

# editable install, rebuilding ONLY PIMMS (--no-deps leaves other env packages
# untouched; --reinstall forces a fresh build and refreshes uv's build cache).
# If you use plain pip instead of uv, the equivalent is:
#     python -m pip install -e . --force-reinstall --no-deps
uv pip install -e . --no-deps --reinstall
