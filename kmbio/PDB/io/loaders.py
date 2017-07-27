import functools
import inspect
import logging
import os.path as op
import re
import string
import warnings

from kmbio.PDB import MMCIFParser, MMTFParser, Parser, PDBParser, Structure, open_url

from .routes import DEFAULT_ROUTES

logger = logging.getLogger(__name__)


def load(pdb_file, **kwargs) -> Structure:
    """Load local PDB file.

    Parameters
    ----------
    pdb_file : :obj:`str`
        File to load
    kwargs : :obj:`dict`
        Optional keyword arguments to be passed to the parser ``__init__`` and ``get_structure``
        methods.

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
    for default_route in DEFAULT_ROUTES:
        if pdb_file.startswith(default_route):
            pdb_filename = pdb_file.partition(default_route)[-1]
            pdb_id, _, pdb_type = pdb_filename.rpartition('.')
            pdb_file = DEFAULT_ROUTES[default_route](pdb_id, pdb_type)
            break

    if pdb_file.startswith('file://'):
        pdb_file = pdb_file.partition('file://')[-1]

    logger.debug('pdb_file: %s', pdb_file)
    pdb_type = guess_pdb_type(pdb_file)
    parser = _get_parser(pdb_type, **kwargs)
    with open_url(pdb_file) as fh:
        structure = parser.get_structure(fh)
    return structure


def guess_pdb_id(pdb_file):
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


def guess_pdb_type(pdb_file):
    """Guess PDB file type from file name.

    Examples
    >>> _guess_pdb_type('4dkl.pdb')
    'pdb'
    >>> _guess_pdb_type('/tmp/4dkl.cif.gz')
    'cif'
    """
    for chunk in re.split('/|\.|:', pdb_file):
        chunk = chunk.lower().strip(string.digits)
        if chunk in ['pdb', 'ent']:
            return 'pdb'
        elif chunk in ['cif', 'mmcif']:
            return 'cif'
        elif chunk in ['mmtf']:
            return 'mmtf'
    raise Exception("Count not guess pdb type for file '{}'!".format(pdb_file))


def _get_parser(pdb_type, **kwargs) -> Parser:
    """Get kmbioPython PDB parser appropriate for `pdb_type`."""
    if pdb_type == 'pdb':
        Parser = PDBParser
    elif pdb_type == 'cif':
        kwargs.setdefault('ignore_auth_id', True)
        Parser = MMCIFParser
    elif pdb_type == 'mmtf':
        Parser = MMTFParser
    else:
        raise Exception("Wrong pdb_type: '{}'".format(pdb_type))
    init_params = set(inspect.signature(Parser).parameters)
    parser = Parser(**{k: kwargs.pop(k) for k in list(kwargs) if k in init_params})
    func_params = set(inspect.signature(parser.get_structure).parameters)
    parser.get_structure = functools.partial(
        parser.get_structure, **{k: kwargs.pop(k)
                                 for k in list(kwargs) if k in func_params})
    if kwargs:
        warnings.warn("Not all arguments where used during the call to _get_parser! (kwargs = {})".
                      format(kwargs))
    return parser
