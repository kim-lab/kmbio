from collections import OrderedDict
from typing import Callable, Dict


def get_rcsb_url(pdb_id: str, pdb_type: str) -> str:
    if pdb_type in ["pdb", "cif"]:
        url = f"http://files.rcsb.org/download/{pdb_id}.{pdb_type}.gz"
    elif pdb_type in ["mmtf"]:
        url = f"http://mmtf.rcsb.org/v1.0/full/{pdb_id}"
    else:
        raise TypeError(f"This route does not support '{pdb_type}' file format!")
    return url


def get_ebi_url(pdb_id: str, pdb_type: str) -> str:
    if pdb_type == "pdb":
        pdb_filename = f"pdb{pdb_id}.ent"
    elif pdb_type == "cif":
        pdb_filename = f"{pdb_id}.cif"
    else:
        raise TypeError(f"This route does not support '{pdb_type}' file format!")
    url = f"http://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_filename}"
    return url


def get_wwpdb_url(pdb_id: str, pdb_type: str) -> str:
    pdb_id_middle = pdb_id[1:3]
    if pdb_type == "pdb":
        pdb_filename = f"{pdb_id_middle}/pdb{pdb_id}.ent.gz"
        pdb_format = "pdb"
    elif pdb_type == "cif":
        pdb_filename = f"{pdb_id_middle}/{pdb_id}.cif.gz"
        pdb_format = "mmCIF"
    else:
        raise TypeError(f"This route does not support '{pdb_type}' file format!")
    url = f"ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/{pdb_format}/{pdb_filename}"
    return url


DEFAULT_ROUTES: Dict[str, Callable[[str, str], str]] = OrderedDict(
    [("rcsb", get_rcsb_url), ("ebi", get_ebi_url), ("wwpdb", get_wwpdb_url)]
)
