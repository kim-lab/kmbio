# flake8: noqa
import warnings

try:
    from ._mmcif_to_dict import mmcif2dict
except ImportError:
    warnings.warn("Cound not import cythonized `mmcif2dict` function. Performance will suffer!")
    from .mmcif_to_dict import MMCIF2Dict as mmcif2dict

from .bioassembly import ProcessRemark350, get_mmcif_bioassembly_data
from .parser import Parser
from .pdb_parser import PDBParser
from .mmcif_parser import MMCIFParser, FastMMCIFParser
from .mmtf_parser import MMTFParser
