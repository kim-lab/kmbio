from collections import OrderedDict
from pathlib import Path
from urllib.parse import urlparse

import pytest

from kmbio.PDB.io.loaders import get_parser
from kmbio.PDB.utils import open_url, sort_ordered_dict

TESTS_DIR = Path(__file__).absolute().parent


@pytest.mark.parametrize(
    "unsorted_dict, sorted_dict",
    [
        (OrderedDict([(1, 10), (3, 30), (2, 20)]), OrderedDict([(1, 10), (2, 20), (3, 30)])),
        (
            OrderedDict([("b", "bye bye"), ("a", "hello")]),
            OrderedDict([("a", "hello"), ("b", "bye bye")]),
        ),
    ],
)
def test_sort_ordered_dict(unsorted_dict, sorted_dict):
    assert unsorted_dict != sorted_dict
    print(unsorted_dict, len(sorted_dict))
    sort_ordered_dict(unsorted_dict)
    assert unsorted_dict == sorted_dict


@pytest.mark.parametrize(
    "pdb_url, pdb_type",
    sum(
        [
            [
                (f"http://files.rcsb.org/download/{pdb_id}.pdb.gz", "pdb"),
                (f"http://files.rcsb.org/download/{pdb_id}.cif.gz", "cif"),
                # (f"http://mmtf.rcsb.org/v1.0/full/{pdb_id}", "mmtf"),
                (f"http://www.ebi.ac.uk/pdbe/entry-files/download/pdb{pdb_id}.ent", "pdb"),
                (f"http://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_id}.cif", "cif"),
                (
                    "ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/"
                    f"{pdb_id[1:3]}/pdb{pdb_id}.ent.gz",
                    "pdb",
                ),
                (
                    "ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/"
                    f"{pdb_id[1:3]}/{pdb_id}.cif.gz",
                    "cif",
                ),
                (f"ff://{TESTS_DIR}/test_utils/structures-subset?{pdb_id}.cif.gz", "cif"),
            ]
            for pdb_id in ["1glv", "1nip", "1ud3"]
        ],
        [],
    ),
)
def test_open_url(pdb_url, pdb_type):
    with open_url(pdb_url) as fh:
        parser = get_parser(pdb_type)
        parser.get_structure(fh)
