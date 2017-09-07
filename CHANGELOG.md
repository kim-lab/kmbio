# Changelog

## 2.0.2

- Entities now keep children inside a single `OrderedDict` named `_children`. Previously there was a list `_child_list` and a dictionary `_child_dict`.
- Iteration now returns objects instead of ids. For example, `for chain in model` iterates over chain objects instead of chain ids. This is consistent with `Bio.PDB`.
- Refactor Cython code in `_mmcif_to_dict.pyx` for > 2x speed gain.
- Add `test_mmcif_vs_mmcif_ref` tests, since now those tests are fast enough (sort of).

## 2.0.1

- Fix lingering bugs with bioassembly generators.

## 2.0.0

- Restructured the repository. Move files into different submodules according to their function.
- Add jack-of-all-trades functions `load` and `save`.
- Add bioassembly generators to `PDBParser` and `MMCIFParser` classes.
- Many other changes (see commits since the fork from biopython).
