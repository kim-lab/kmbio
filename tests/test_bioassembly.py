import contextlib
import gzip
import logging
import os
import tempfile
import urllib.request

import pytest

from kmbio.PDB import allequal, MMCIFParser, PDBParser

logger = logging.getLogger(__name__)

URL = "ftp://ftp.wwpdb.org/pub/pdb/data/"

# (pdb_id, bioassembly_id)
TEST_DATA = [
    ('1y0y', 1),
]


@contextlib.contextmanager
def gzip_open_url(url, mode='w+t'):
    with tempfile.TemporaryFile(mode) as ofh:
        with urllib.request.urlopen(url) as ifh:
            ofh.write(gzip.decompress(ifh.read()).decode('utf-8'))
        ofh.seek(0)
        yield ofh
    return None


@pytest.mark.parametrize("pdb_id, bioassembly_id", TEST_DATA)
def test_mmcif_to_pdb(pdb_id, bioassembly_id):
    mmcif_url = (URL + "structures/divided/mmCIF/{}/{}.cif.gz".format(pdb_id[1:3], pdb_id))
    logger.info(mmcif_url)

    pdb_bioassembly_url = (
        URL + "biounit/PDB/divided/{}/{}.pdb{}.gz".format(pdb_id[1:3], pdb_id, bioassembly_id))
    logger.info(pdb_bioassembly_url)

    with gzip_open_url(mmcif_url) as fh:
        mmcif_structure = MMCIFParser(ignore_auth_id=False).get_structure(
            fh, bioassembly_id=bioassembly_id)

    with gzip_open_url(pdb_bioassembly_url) as fh:
        pdb_bioassembly_structure = PDBParser().get_structure(fh)

    assert allequal(mmcif_structure, pdb_bioassembly_structure)


# PDB_DATA_DIR = '/home/kimlab1/database_data/pdb/data/data/'
@pytest.mark.skipif(
    'PDB_DATA_DIR' not in os.environ,
    reason="set PDB_DATA_DIR environment variable to run this test!")
@pytest.mark.parametrize("pdb_id, bioassembly_id, ignore_auth_id",
                         [(*t, tf) for t in TEST_DATA for tf in [True, False]])
def test_mmcif_to_mmcif(pdb_id, bioassembly_id, ignore_auth_id):
    mmcif_file = (os.environ['PDB_DATA_DIR'] + "structures/divided/mmCIF/{}/{}.cif.gz".format(
        pdb_id[1:3], pdb_id))
    logger.info(mmcif_file)

    mmcif_bioassembly_file = (
        os.environ['PDB_DATA_DIR'] + "structures/divided/mmCIF/{}/{}-{}.cif".format(
            pdb_id[1:3], pdb_id, bioassembly_id))
    logger.info(mmcif_bioassembly_file)

    with gzip.open(mmcif_file, 'rt') as ifh:
        mmcif_structure = MMCIFParser(ignore_auth_id=False).get_structure(
            ifh, bioassembly_id=bioassembly_id)

    with open(mmcif_bioassembly_file, 'rt') as ifh:
        mmcif_bioassembly_structure = MMCIFParser(ignore_auth_id=False).get_structure(ifh)

    assert allequal(mmcif_structure, mmcif_bioassembly_structure)
