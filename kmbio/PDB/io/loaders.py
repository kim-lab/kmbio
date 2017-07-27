import functools
import io
import logging
import os.path as op

from kmbio.PDB import MMCIFParser, Model, PDBParser, Structure
from kmtools import system_tools

logger = logging.getLogger(__name__)


def load(pdb_file, pdb_id=None, pdb_type='cif', biounit=False, pdb_mirror=None):
    """Load local PDB file.

    Parameters
    ----------
    pdb_file : :obj:`str`
        File to load
    pdb_id : :obj:`str`, optional
        PDB_ID to assign to the loaded structure.
    pdb_type: :obj:`str`, one of: `{'pdb', 'cif'}`, optional
        Type of PDB file. If ``None``, try all available types until one works.

    Returns
    -------
    :class:`kmbio.PDB.Structure`
        Protein structure.

    Load Examples
    -------------
    >>> import urllib.request
    >>> pdb_file = op.join(tempfile.gettempdir(), '4dkl.pdb')
    >>> r = urllib.request.urlretrieve('http://files.rcsb.org/download/4dkl.pdb', pdb_file)
    >>> load(pdb_file)
    <Structure id=4dkl>

    Fetch Examples
    --------------
    >>> load('wwpdb://4dkl')
    <Structure id=4dkl>
    >>> load('wwpdb://4dkl.cif')
    <Structure id=4dkl>
    """
    if pdb_file.startswith(('rcsb://', 'wwpdb://')):
        pdb_data = _fetch_structure(pdb_id, pdb_type, biounit, pdb_mirror)
        parser = get_pdb_parser(pdb_type)
        structure = parser.get_structure(io.StringIO(pdb_data.decode('utf-8')), pdb_id)
        return structure

    if pdb_id is None:
        pdb_id = _guess_pdb_id(pdb_file)

    if pdb_type is not None:
        pdb_types = [pdb_type]
    else:
        pdb_types = ['pdb', 'cif']

    with system_tools.open_compressed(pdb_file, mode='rt') as ifh:
        for pdb_type in pdb_types:
            parser = get_pdb_parser(pdb_type)
            try:
                structure = parser.get_structure(ifh, pdb_id)
                break
            except KeyError:
                logger.info("Count not load structure using the %s parser.", parser)
                continue
    return structure


def _guess_pdb_id(pdb_file):
    """Extract the PDB id from a PDB file.

    Examples
    --------
    >>> _guess_pdb_id('4dkl.pdb')
    '4dkl'
    >>> _guess_pdb_id('/data/structures/divided/pdb/26/pdb126d.ent.gz')
    '126d'
    >>> _guess_pdb_id('/tmp/100d.cif.gz')
    '100d'
    """
    pdb_id = op.basename(pdb_file)
    for extension in ['.gz', '.pdb', '.ent', '.cif']:
        pdb_id = pdb_id.partition(extension)[0]
    if len(pdb_id) == 7 and (pdb_id.startswith('ent') or pdb_id.startswith('pdb')):
        pdb_id = pdb_id[3:]
        assert len(pdb_id) == 4
    pdb_id = pdb_id.lower()
    pdb_id = pdb_id.replace('.', '')
    return pdb_id


def _guess_pdb_type(pdb_file):
    """Guess PDB file type from file name.

    Examples
    >>> _guess_pdb_type('4dkl.pdb')
    'pdb'
    >>> _guess_pdb_type('/tmp/4dkl.cif.gz')
    'cif'
    """
    file_name, file_ext = op.splitext(pdb_file)
    while file_ext:
        file_ext = file_ext.lower()
        if file_ext in ['.pdb']:
            return 'pdb'
        if file_ext in ['.cif', '.mmcif']:
            return 'cif'
        file_name, file_ext = op.splitext(file_name)
    return None


def _get_wwpdb_url(pdb_id, pdb_type='cif', biounit_id=0, wwpdb_base_url=None):
    """Get PDB filename (inlcuding path) from a local mirror of the PDB repository.

    This assumes that you have mirrored the PDB ftp repository:
    ftp://ftp.wwpdb.org/pub/pdb/data/

    Examples
    --------
    >>> _get_wwpdb_url('4DKL', 'pdb', biounit=False, wwpdb_base_url='')
    'structures/divided/pdb/dk/pdb4dkl.ent.gz'
    >>> _get_wwpdb_url('4dkl', 'cif', biounit=False, wwpdb_base_url='')
    'structures/divided/mmCIF/dk/4dkl.cif.gz'
    >>> _get_wwpdb_url('4DKL', 'pdb', biounit=True, wwpdb_base_url='/tmp')
    '/tmp/biounit/PDB/divided/dk/4dkl.pdb1.gz'
    >>> _get_wwpdb_url('4NWR', 'cif', biounit=True)
    'ftp://ftp.wwpdb.org/pub/pdb/data/biounit/mmCIF/divided/nw/4nwr-assembly1.cif.gz'
    """
    if wwpdb_base_url is None:
        wwpdb_base_url = 'ftp://ftp.wwpdb.org/pub/pdb/data/'

    if pdb_type.startswith('pdb') and not biounit_id:
        subfolder = op.join('structures', 'divided', 'pdb')
        filename = 'pdb{}.ent.gz'.format(pdb_id.lower())
    elif pdb_type.startswith('pdb') and biounit_id:
        subfolder = op.join('biounit', 'PDB', 'divided')
        filename = '{}.pdb1.gz'.format(pdb_id.lower())
    elif pdb_type.startswith('cif') and not biounit_id:
        subfolder = op.join('structures', 'divided', 'mmCIF')
        filename = '{}.cif.gz'.format(pdb_id.lower())
    elif pdb_type.startswith('cif') and biounit_id:
        subfolder = op.join('biounit', 'mmCIF', 'divided')
        filename = '{}-assembly1.cif.gz'.format(pdb_id.lower())
    else:
        raise Exception

    pdb_filename = op.join(wwpdb_base_url, subfolder, pdb_id[1:3].lower(), filename)
    return pdb_filename


def _get_pdb_url(pdb_id, pdb_type='cif', mirror='rcsb'):
    """Download PDB from RCSB or EBI website.

    .. todo:: Add an option for biounit urls.

    .. note:: Not used at the moment.
    """
    if mirror == 'rcsb':
        return 'http://www.rcsb.org/pdb/files/{}.pdb'.format(pdb_id.lower())
    elif mirror == 'ebi':
        return 'http://www.ebi.ac.uk/pdbe/entry-files/download/pdb{}.ent'.format(pdb_id.lower())
    else:
        raise Exception("Unsupported mirror: {}.".format(mirror))


def get_pdb_parser(pdb_type):
    """Get kmbioPython PDB parser appropriate for `pdb_type`."""
    if pdb_type == 'pdb':
        return PDBParser()
    elif pdb_type == 'cif':
        return MMCIFParser()
    else:
        raise Exception


@functools.lru_cache(maxsize=512)
def _fetch_structure(pdb_id, pdb_type='cif', biounit_id=0, pdb_mirror=None):
    url = _get_wwpdb_url(
        pdb_id, pdb_type=pdb_type, biounit_id=biounit_id, wwpdb_base_url=pdb_mirror)
    logger.debug("url: %s", url)

    pdb_data = system_tools.read_url(url)

    return pdb_data


def structure_from_chain(structure_id, model_id, chain):
    structure = Structure(structure_id)
    model = Model(model_id)
    structure.add(model)
    model.add(chain.copy())
    return structure
