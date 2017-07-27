from collections import OrderedDict


def _get_rcsb_url(pdb_id, pdb_type):
    if pdb_type in ['pdb', 'cif']:
        URL = 'http://files.rcsb.org/download/{pdb_id}.{pdb_type}.gz'
    elif pdb_type in ['mmtf']:
        URL = 'http://mmtf.rcsb.org/v1.0/full/{pdb_id}'
    else:
        raise Exception("Wrong pdb_type: '{}'".format(pdb_type))
    url = URL.format(pdb_id=pdb_id, pdb_type=pdb_type)
    return url


def _get_ebi_url(pdb_id, pdb_type):
    if pdb_type in ['mmtf']:
        raise NotImplementedError("This route does not support '{}' file format!".format(pdb_type))
    URL = 'http://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_filename}'
    pdb_filename = {
        'pdb': 'pdb{pdb_id}.ent',
        'cif': '{pdb_id}.cif'
    }[pdb_type].format(pdb_id=pdb_id)
    url = URL.format(pdb_filename=pdb_filename)
    return url


def _get_wwpdb_url(pdb_id, pdb_type):
    if pdb_type in ['mmtf']:
        raise NotImplementedError("This route does not support '{}' file format!".format(pdb_type))
    URL = 'ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/{pdb_format}/{pdb_filename}'
    pdb_format = {'pdb': 'pdb', 'cif': 'mmCIF'}[pdb_type]
    pdb_filename = {
        'pdb': '{pdb_id_middle}/pdb{pdb_id}.ent.gz',
        'cif': '{pdb_id_middle}/{pdb_id}.cif.gz'
    }[pdb_type].format(
        pdb_id_middle=pdb_id[1:3], pdb_id=pdb_id)
    url = URL.format(pdb_format=pdb_format, pdb_filename=pdb_filename)
    return url


DEFAULT_ROUTES = OrderedDict([
    ('rcsb://', _get_rcsb_url),
    ('ebi://', _get_ebi_url),
    ('wwpdb://', _get_wwpdb_url),
])
