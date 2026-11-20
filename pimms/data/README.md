# pimms/data

Two data files that ship inside the installed package (they are swept in by `graft pimms` in `MANIFEST.in`). The only code that reads anything here is `pdb2pimms`, which uses `gcf_rje23_v14.prm`; a packaging test (`pimms/tests/test_package_hygiene.py`) checks that `MANIFEST.in` keeps shipping `look_and_say.dat`.

## Contents

* `gcf_rje23_v14.prm` - the residue-level amino-acid parameter file used by the single-chain protein demo; it is identical to `demo_keyfiles/single_chain_protein_demo/gcf_rje23_v14.prm`, which is the copy the demo actually uses. It was parameterised at `TEMPERATURE : 50`. `pdb2pimms --protein-mode` copies it next to the keyfile it writes and runs it at that temperature. To use it yourself, point `PARAMETER_FILE` at it.
* `look_and_say.dat` - the first entries of the "Look and Say" integer series (OEIS [A005150](https://oeis.org/A005150)). It is a placeholder from the cookiecutter template the repository was created from and has nothing to do with the PIMMS model.
