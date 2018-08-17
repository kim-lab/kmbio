import functools
import inspect
import logging
import os.path as op
import string
import warnings
from pathlib import Path
from typing import Type
from urllib.parse import urlparse

from kmbio.PDB import MMCIFParser, MMTFParser, Parser, PDBParser, Structure, open_url

from .routes import DEFAULT_ROUTES

logger = logging.getLogger(__name__)


def load(pdb_file: str, structure_id: str = None, **kwargs) -> Structure:
    """Load local PDB file.

    Args:
        pdb_file: File to load.
        kwargs: Optional keyword arguments to be passed to the parser
            ``__init__`` and ``get_structure`` methods.

    Load example:
        >>> import urllib.request
        >>> pdb_file = op.join(tempfile.gettempdir(), '4dkl.pdb')
        >>> r = urllib.request.urlretrieve('http://files.rcsb.org/download/4dkl.pdb', pdb_file)
        >>> load(pdb_file)
        <Structure id=4dkl>

    Fetch example:
        >>> load('wwpdb://4dkl')
        <Structure id=4dkl>
        >>> load('wwpdb://4dkl.cif')
        <Structure id=4dkl>
    """
    if isinstance(pdb_file, Path):
        pdb_file = pdb_file.as_posix()

    pdb_id = guess_pdb_id(pdb_file)
    pdb_type = guess_pdb_type(pdb_file)

    scheme = urlparse(pdb_file).scheme
    if scheme in DEFAULT_ROUTES:
        pdb_file = DEFAULT_ROUTES[scheme](pdb_id, pdb_type)

    parser = get_parser(pdb_type, **kwargs)

    with open_url(pdb_file) as fh:
        structure = parser.get_structure(fh)
        if not structure.id:
            structure.id = pdb_id

    return structure


def guess_pdb_id(pdb_file: str) -> str:
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
    for extension in [".gz", ".pdb", ".ent", ".cif"]:
        pdb_id = pdb_id.partition(extension)[0]
    if len(pdb_id) == 7 and (pdb_id.startswith("ent") or pdb_id.startswith("pdb")):
        pdb_id = pdb_id[3:]
        assert len(pdb_id) == 4
    pdb_id = pdb_id.lower()
    pdb_id = pdb_id.replace(".", "")
    return pdb_id


def guess_pdb_type(pdb_file: str) -> str:
    """Guess PDB file type from file name.

    Examples
    --------
    >>> _guess_pdb_type('4dkl.pdb')
    'pdb'
    >>> _guess_pdb_type('/tmp/4dkl.cif.gz')
    'cif'
    """
    for suffix in reversed(Path(pdb_file).suffixes):
        suffix = suffix.lower().strip(string.digits)
        if suffix in [".pdb", ".ent"]:
            return "pdb"
        elif suffix in [".cif", ".mmcif"]:
            return "cif"
        elif suffix in [".mmtf"]:
            return "mmtf"
    raise Exception(f"Could not guess pdb type for file '{pdb_file}'!")


def get_parser(pdb_type: str, **kwargs) -> Parser:
    """Get `kmbio.PDB` parser appropriate for `pdb_type`."""
    MyParser: Type[Parser]
    if pdb_type == "pdb":
        MyParser = PDBParser
    elif pdb_type == "cif":
        kwargs.setdefault("use_auth_id", False)
        MyParser = MMCIFParser
    elif pdb_type == "mmtf":
        MyParser = MMTFParser
    else:
        raise Exception("Wrong pdb_type: '{}'".format(pdb_type))
    init_params = set(inspect.signature(MyParser).parameters)
    parser = MyParser(  # type: ignore
        **{k: kwargs.pop(k) for k in list(kwargs) if k in init_params}
    )
    func_params = set(inspect.signature(parser.get_structure).parameters)
    parser.get_structure = functools.partial(  # type: ignore
        parser.get_structure, **{k: kwargs.pop(k) for k in list(kwargs) if k in func_params}
    )
    if kwargs:
        warnings.warn(
            f"Not all arguments where used during the call to _get_parser! (kwargs = {kwargs})"
        )
    return parser
