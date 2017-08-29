import logging
import os

import pytest

from common import (ATOM_DEFINED_TWICE_PDBS, LOCAL_REMOTE_MISMATCH, MISSING,
                    PDB_IDS, random_subset)
from kmbio.PDB import allequal, DEFAULT_ROUTES, load
from kmbio.PDB.exceptions import BioassemblyError
from kmbio.PDB.io.loaders import guess_pdb_type

logger = logging.getLogger(__name__)


def _get_local_url(pdb_id, pdb_type):
    if pdb_type in ['mmtf']:
        raise NotImplementedError("This route does not support '{}' file format!".format(pdb_type))
    URL = '{pdb_data_dir}/structures/divided/{pdb_format}/{pdb_filename}'
    # export PDB_DATA_DIR='/home/kimlab1/database_data/pdb/data/data/'
    pdb_data_dir = os.getenv("PDB_DATA_DIR")
    pdb_format = {'pdb': 'pdb', 'cif': 'mmCIF'}[pdb_type]
    pdb_filename = {
        'pdb': '{pdb_id_middle}/pdb{pdb_id}.ent.gz',
        'cif': '{pdb_id_middle}/{pdb_id}.cif.gz'
    }[pdb_type].format(
        pdb_id_middle=pdb_id[1:3], pdb_id=pdb_id)
    url = URL.format(pdb_data_dir=pdb_data_dir, pdb_format=pdb_format, pdb_filename=pdb_filename)
    return url


@pytest.mark.parametrize(
    "url, pdb_type", [('ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/dk/4dkl.cif.gz',
                       'cif')])
def test_guess_pdb_type(url, pdb_type):
    assert guess_pdb_type(url) == pdb_type


@pytest.mark.parametrize("pdb_id, pdb_type, bioassembly_id",
                         [(pdb_id, pdb_type, bioassembly_id)
                          for pdb_id in PDB_IDS for pdb_type in ['pdb', 'cif']
                          for bioassembly_id in ([0] if pdb_type == 'pdb' else [0, 1])
                          if (pdb_id, pdb_type) not in LOCAL_REMOTE_MISMATCH
                          if (pdb_id, pdb_type, bioassembly_id) not in MISSING])
def test_equal(pdb_id, pdb_type, bioassembly_id):
    """Make sure that loading local and remote files produces the same result."""
    filename = '{}.{}'.format(pdb_id, pdb_type)
    logger.debug(filename)
    structures = []
    exceptions = []
    for route in DEFAULT_ROUTES:
        try:
            structure = load(route + filename, bioassembly_id=bioassembly_id)
            structures.append(structure)
            exceptions.append(None)
        except BioassemblyError as exception:
            structures.append(None)
            exceptions.append(str(type(exception)))
    if any(s is not None for s in structures):
        assert all(e is None for e in exceptions), (list(DEFAULT_ROUTES), structures, exceptions)
        for s in structures[1:]:
            assert allequal(structure, structures[0])
    else:
        for e in exceptions:
            assert e == exceptions[0]


@pytest.mark.parametrize("pdb_id_1, pdb_id_2, pdb_type",
                         [(pdb_id_1, pdb_id_2, pdb_type)
                          for pdb_id_1 in PDB_IDS for pdb_id_2 in PDB_IDS
                          for pdb_type in ['pdb', 'cif']
                          if pdb_id_1 != pdb_id_2 if (pdb_id_1, pdb_type, 0) not in MISSING
                          if (pdb_id_2, pdb_type, 0) not in MISSING])
def test_notequal(pdb_id_1, pdb_id_2, pdb_type):
    """Make sure that structures that should be different are different."""
    s1 = load('rcsb://{}.{}'.format(pdb_id_1, pdb_type))
    s2 = load('rcsb://{}.{}'.format(pdb_id_2, pdb_type))
    assert not allequal(s1, s2)


@pytest.mark.parametrize("pdb_id", random_subset(ATOM_DEFINED_TWICE_PDBS))
def test_atom_defined_twice(pdb_id):
    """Tests for the ``Atom defined twice`` error."""
    s = load('rcsb://{}.{}'.format(pdb_id, 'cif'))
    assert s
