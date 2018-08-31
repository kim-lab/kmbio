# History

## 2.0.8 (2018-08-31)

* Allow chain ids of length 2 when saving PDB files.

## 2.0.7 (2018-08-17)

* Allow loading structures from `ffindex` archives using the `ff://` scheme. For example:

    ```python
    structure = kmbio.PDB.load("ff:///path/to/structures?4dkl.cif.gz")
    ```

* So that prettier leaves *this* alone...

## 2.0.6 (2018-08-06)

* Added `history.md` file.
