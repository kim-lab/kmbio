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
    ('2co0', 1),
]


@pytest.mark.parametrize("pdb_id, bioassembly_id", TEST_DATA)
def test_mmcif_to_pdb(pdb_id, bioassembly_id):
    mmcif_url = (URL + "structures/divided/mmCIF/{}/{}.cif.gz".format(pdb_id[1:3], pdb_id))
    logger.info(mmcif_url)

    pdb_bioassembly_url = (
        URL + "biounit/PDB/divided/{}/{}.pdb{}.gz".format(pdb_id[1:3], pdb_id, bioassembly_id))
    logger.info(pdb_bioassembly_url)

    with urllib.request.urlopen(mmcif_url) as ifh, \
            tempfile.TemporaryFile('w+t') as ofh:
        ofh.write(gzip.decompress(ifh.read()).decode('utf-8'))
        ofh.seek(0)
        mmcif_structure = MMCIFParser(ignore_auth_id=False).get_structure(
            ofh, bioassembly_id=bioassembly_id)

    with urllib.request.urlopen(pdb_bioassembly_url) as ifh, \
            tempfile.TemporaryFile('w+t') as ofh:
        ofh.write(gzip.decompress(ifh.read()).decode('utf-8'))
        ofh.seek(0)
        pdb_bioassembly_structure = PDBParser().get_structure(ofh)

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
