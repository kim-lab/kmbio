import re
import warnings

import pandas as pd

from kmbio.PDB import Model, Structure
from kmbio.PDB.exceptions import BioassemblyNotFoundError


class ProcessRemark350:

    RE_BIOMOLECULE = re.compile('BIOMOLECULE: ([0-9]+)')
    RE_CHAINS = re.compile('APPLY THE FOLLOWING TO CHAINS: ([a-zA-Z0-9, ]+)')
    RE_BIOMT = re.compile(
        'BIOMT([0-9]+)\s+([0-9]+)\s+([0-9\.]+)\s+([0-9\.]+)\s+([0-9\.]+)\s+([0-9\.]+)')

    def __init__(self):
        self._biomolecule = None
        self._chains = None
        self._biomt = None

        self.bioassembly_data = {}

    def process_lines(self, lines):
        for line in lines:
            assert line.startswith('REMARK 350 ')
            line = line[11:].strip()
            biomolecule = self.RE_BIOMOLECULE.findall(line)
            if biomolecule:
                self._process_biomolecule(biomolecule[0])
                continue
            chains = self.RE_CHAINS.findall(line)
            if chains:
                self._process_chains(chains[0])
                continue
            biomt = self.RE_BIOMT.findall(line)
            if biomt:
                self._process_biomt(biomt[0])
                continue
        self._flush()

    def _process_biomolecule(self, biomolecule):
        if self._biomolecule is None:
            assert not self._chains and not self._biomt and not self.bioassembly_data
            self._biomolecule = biomolecule
        else:
            assert self._chains and self._biomt
            self._biomolecule = biomolecule

    def _process_chains(self, chains):
        assert not self._chains
        self._chains = chains.replace(' ', '')

    def _process_biomt(self, biomt):
        matrix_id, transformation_id, matrix_0, matrix_1, matrix_2, vector = biomt
        assert matrix_id in ['1', '2', '3']
        if matrix_id == '1':
            row = {}
            self._biomt.append(row)
        else:
            row = biomt[-1]
            # Sanity check on the data that we already have
            assert all([
                '_pdbx_struct_oper_list.matrix[{}][{}]'.format(k, i) in row
                for k in range(1, int(matrix_id)) for i in range(1, 3)
            ])
            assert all(
                ['_pdbx_struct_oper_list.vector[{}]'.format(k) for k in range(1, int(matrix_id))])
        row['_pdbx_struct_oper_list.matrix[{}][1]'.format(matrix_id)] = matrix_0
        row['_pdbx_struct_oper_list.matrix[{}][2]'.format(matrix_id)] = matrix_1
        row['_pdbx_struct_oper_list.matrix[{}][3]'.format(matrix_id)] = matrix_2
        row['_pdbx_struct_oper_list.vector[{}]'.format(matrix_id)] = vector

    def _flush(self):
        assert str(self._biomolecule) not in self.bioassembly_data
        self.bioassembly_data[str(self._biomolecule)] = {
            'chains': self._chains,
            'transformations': [((get_rotation(row), get_translation(row)) for row in self._biomt)]
        }
        self._chains = None
        self._biomt = None


def mmcif_key_to_dataframe(sdict, key):
    """Convert all entries in ``sdict`` with keys starting with ``key`` into a `pandas.DataFrame`.

    Parameters
    ----------
    sdict : `dict`
        Dictionary of of MMCIF elements.
    key : `str`
        Element of the dictionary to convert to a dataframe.

    Examples
    --------
    >>> mmcif_key_to_dataframe({'a.1': [1, 2], 'a.2': [3, 4]}, 'a')
       a.1  a.2
    0    1    3
    1    2    4
    """
    data = {k: v for k, v in sdict.items() if k.startswith(key + '.')}
    if not data:
        return pd.DataFrame()
    data_element = next(iter(data.values()))
    if isinstance(data_element, (list, tuple)):
        assert all(len(data_element) == len(v) for v in data.values())
        _data = []
        for row in zip(* [v for v in data.values()]):
            _data.append(dict(zip(data.keys(), row)))
        data = _data
    else:
        data = [data]
    return pd.DataFrame(data)


def get_rotation(row):
    """Generate a rotation matrix from elements in dictionary ``row``.

    Examples
    ---------
    >>> sdict = { \
        '_pdbx_struct_oper_list.matrix[{}][{}]'.format(i // 3 + 1, i % 3 + 1): i \
        for i in range(9) \
    }
    >>> get_rotation(sdict)
    [[0.0, 1.0, 2.0], [3.0, 4.0, 5.0], [6.0, 7.0, 8.0]]
    """
    return [[
        float(row['_pdbx_struct_oper_list.matrix[1][1]']),
        float(row['_pdbx_struct_oper_list.matrix[2][1]']),
        float(row['_pdbx_struct_oper_list.matrix[3][1]'])
    ], [
        float(row['_pdbx_struct_oper_list.matrix[1][2]']),
        float(row['_pdbx_struct_oper_list.matrix[2][2]']),
        float(row['_pdbx_struct_oper_list.matrix[3][2]'])
    ], [
        float(row['_pdbx_struct_oper_list.matrix[1][3]']),
        float(row['_pdbx_struct_oper_list.matrix[2][3]']),
        float(row['_pdbx_struct_oper_list.matrix[3][3]'])
    ]]


def get_translation(row):
    """Generate a translation matrix from elements in dictionary ``row``.

    Examples
    ---------
    >>> sdict = {'_pdbx_struct_oper_list.vector[{}]'.format(i): i for i in range(1, 4)}
    >>> get_translation(sdict)
    [1.0, 2.0, 3.0]
    """
    return [
        float(row['_pdbx_struct_oper_list.vector[1]']),
        float(row['_pdbx_struct_oper_list.vector[2]']),
        float(row['_pdbx_struct_oper_list.vector[3]']),
    ]


def get_chain_and_transformation_ids(sdict, bioassembly_id, ignore_auth_id=True):
    try:
        _pdbx_struct_assembly_gen = \
            mmcif_key_to_dataframe(sdict, '_pdbx_struct_assembly_gen') \
            .set_index('_pdbx_struct_assembly_gen.assembly_id') \
            .loc[str(bioassembly_id)]
    except KeyError:
        raise BioassemblyNotFoundError
    chain_ids = \
        _pdbx_struct_assembly_gen['_pdbx_struct_assembly_gen.asym_id_list'] \
        .split(',')
    transformation_ids = \
        _pdbx_struct_assembly_gen['_pdbx_struct_assembly_gen.oper_expression'] \
        .split(',')

    if not ignore_auth_id:
        label_id_to_auth_id_df = \
            mmcif_key_to_dataframe(sdict, '_atom_site')[
                ['_atom_site.label_asym_id', '_atom_site.auth_asym_id']
            ].drop_duplicates()
        label_id_to_auth_id_map = dict(label_id_to_auth_id_df.values.tolist())
        if len(label_id_to_auth_id_map) != len(label_id_to_auth_id_df):
            raise Exception("Cound not reliably map 'label_asym_id' to 'auth_asym_id'!")
        chain_ids = sorted({label_id_to_auth_id_map[chain_id] for chain_id in chain_ids})
    return chain_ids, transformation_ids


def get_transformations(sdict, transformation_ids):
    _pdbx_struct_oper_list = \
        mmcif_key_to_dataframe(sdict, '_pdbx_struct_oper_list') \
        .set_index('_pdbx_struct_oper_list.id') \
        .loc[transformation_ids]
    _pdbx_struct_oper_list['rotation'] = \
        _pdbx_struct_oper_list.apply(get_rotation, axis=1)
    _pdbx_struct_oper_list['translation'] = \
        _pdbx_struct_oper_list.apply(get_translation, axis=1)
    transformations = \
        _pdbx_struct_oper_list[['rotation', 'translation']] \
        .values.tolist()
    return transformations


def transform_structure(structure, chain_ids, transformations):
    bioassembly = Structure(structure.id)
    for i, (rotation, translation) in enumerate(transformations):
        model = Model(i)
        bioassembly.add(model)
        for chain_id in sorted(chain_ids):
            chain = structure[0][chain_id].copy()
            model.add(chain)
        model.transform(rotation, translation)
    return bioassembly


def generate_bioassembly(sdict, structure, bioassembly_id, ignore_auth_id=True):
    if bioassembly_id == 0:
        warnings.warn("bioassembly with id {} is the original structure!".format(bioassembly_id))
        return structure
    chain_ids, transformation_ids = get_chain_and_transformation_ids(sdict, bioassembly_id,
                                                                     ignore_auth_id)
    transformations = get_transformations(sdict, transformation_ids)
    structure = transform_structure(structure, chain_ids, transformations)
    return structure
