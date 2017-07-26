import datetime
import gzip
import logging

import pytest
import urllib.request

from kmbio.PDB.parsers._mmcif_to_dict import MMCIF2Dict as mmcif_to_dict
# from kmbio.PDB.parsers.mmcif_to_dict import MMCIF2Dict as mmcif_to_dict  # for reference

logger = logging.getLogger(__name__)

WWPDB_URL = "ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/mmCIF"

# Times are for the reference Python MMCIF2Dict
PDB_ID_REFS = [
    ('4dkl', 3),
    ('4ug0', 100),
    ('4v6x', 120),
    ('5lks', 110),
    ('1arr', 0.5),
    ('1dvf', 2),
    ('3mbp', 1.5),
    ('4p6f', 250),
]
MULTIPLIER = 1


@pytest.fixture(params=PDB_ID_REFS)
def cif_file_ref_time(request, tmpdir_factory):
    cif_id, ref_time = request.param
    url = WWPDB_URL + "/{}/{}.cif.gz".format(cif_id[1:3], cif_id)
    fn = tmpdir_factory.mktemp('cif_data').join(cif_id + '.cif')
    with urllib.request.urlopen(url) as ifh, open(str(fn), 'wb') as ofh:
        ofh.write(gzip.decompress(ifh.read()))
    return str(fn), ref_time


def test_speed(cif_file_ref_time):
    cif_file, ref_time = cif_file_ref_time
    start_time = datetime.datetime.now()
    mmcif_to_dict(cif_file)
    running_time = (datetime.datetime.now() - start_time).total_seconds()
    logger.info(
        "Processing file '%s' took %s seconds (%.2f faster than Python).",
        cif_file, running_time, ref_time / running_time)
    # assert running_time < (ref_time / MULTIPLIER)
