# History

## 2.0.7 (2018-08-17)

* Allow loading structures from `ffindex` archives using the `ff://` scheme. For example:

    ```python
    structure = kmbio.PDB.load("ff:///path/to/structures?4dkl.cif.gz")
    ```

* So that prettier leaves *this* alone...

## 2.0.6 (2018-08-06)

* Added `history.md` file.
