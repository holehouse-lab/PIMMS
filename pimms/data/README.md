# pimms/data

Two data files that ship inside the installed package (they are swept in by `graft pimms` in `MANIFEST.in`). Nothing in PIMMS reads either of them: no module, test, documentation page or demo loads anything from this directory. The only reference to it in code is a packaging test (`pimms/tests/test_package_hygiene.py`) that checks `MANIFEST.in` keeps shipping `look_and_say.dat`.

## Contents

* `gcf_rje23_v14.prm` - a copy of the residue-level amino-acid parameter file used by the single-chain protein demo; it is identical to `demo_keyfiles/single_chain_protein_demo/gcf_rje23_v14.prm`, which is the copy the demo actually uses. To use this one, point `PARAMETER_FILE` at it.
* `look_and_say.dat` - the first entries of the "Look and Say" integer series (OEIS [A005150](https://oeis.org/A005150)). It is a placeholder from the cookiecutter template the repository was created from and has nothing to do with the PIMMS model.
