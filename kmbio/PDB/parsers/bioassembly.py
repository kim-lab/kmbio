import logging
import re
from collections import OrderedDict, namedtuple

import numpy as np
import pandas as pd

from kmbio.PDB import BioassemblyError, Model, Structure, uniqueify

logger = logging.getLogger(__name__)

Transformation = namedtuple('Transformation', 'transformation_id, rotation, translation')


# === Common ===
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


def apply_bioassembly(structure, bioassembly_data):
    logger.debug("apply_bioassembly(%s, %s)", structure, bioassembly_data)
    bioassembly = Structure(structure.id)
    for chain_id, transformations in bioassembly_data.items():
        for transformation_id, rotation, translation in transformations:
            transformation_idx = transformation_id - 1
            rotation = np.array(rotation, dtype=np.float64)
            translation = np.array(translation, dtype=np.float64)
            try:
                model = bioassembly[transformation_idx]
            except KeyError:
                model = Model(transformation_idx)
                bioassembly.add(model)
            chain = structure[0][chain_id].copy()
            chain.transform(rotation, translation)
            model.add(chain)
    return bioassembly


# === PDB ===
class ProcessRemark350:

    RE_BIOMOLECULE = re.compile('^BIOMOLECULE: +([0-9]+)')
    RE_CHAINS = re.compile('^APPLY THE FOLLOWING TO CHAINS: +([a-zA-Z0-9, ]+)')
    RE_CHAINS_EXTRA = re.compile('^AND CHAINS: +([a-zA-Z0-9, ]+)')
    RE_BIOMT = re.compile(
        '^BIOMT([0-9]+)\s+([0-9]+)\s+([-0-9\.]+)\s+([-0-9\.]+)\s+([-0-9\.]+)\s+([-0-9\.]+)')

    def __init__(self):
        self._biomolecule = None
        self._chains = None
        self._biomt = None

        self.bioassembly_data = OrderedDict()

    def process_lines(self, lines):
        for line in lines:
            assert line.startswith('REMARK 350')
            line = line[11:].strip()
            logger.debug('[%s]', line)
            biomolecule = self.RE_BIOMOLECULE.findall(line)
            if biomolecule:
                self._process_biomolecule(biomolecule[0])
                continue
            chains = self.RE_CHAINS.findall(line)
            if chains:
                self._process_chains(chains[0])
                continue
            chains_extra = self.RE_CHAINS_EXTRA.findall(line)
            if chains_extra:
                self._process_chains_extra(chains_extra[0])
                continue
            biomt = self.RE_BIOMT.findall(line)
            if biomt:
                self._process_biomt(biomt[0])
                continue
        self._flush()
        return self.bioassembly_data

    def _process_biomolecule(self, biomolecule):
        logger.debug('_process_biomolecule(%s)', biomolecule)
        if self._biomolecule is None:
            assert not self._chains and not self._biomt and not self.bioassembly_data
            self._biomolecule = biomolecule
        else:
            assert self._chains and self._biomt
            self._flush()
            self._biomolecule = biomolecule

    def _process_chains(self, chains):
        logger.debug('_process_chains(%s)', chains)
        if self._chains is None:
            assert self._biomt is None
            self._chains = []
        else:
            assert len(self._chains) == len(self._biomt)
        self._chains.append(chains.strip(', ').replace(' ', '').split(','))

    def _process_chains_extra(self, chains_extra):
        logger.debug('_process_chains_extra(%s)', chains_extra)
        self._chains[-1].extend(chains_extra.strip(', ').replace(' ', '').split(','))

    def _process_biomt(self, biomt):
        logger.debug('_process_biomt(%s)', biomt)
        matrix_id, transformation_id, matrix_0, matrix_1, matrix_2, vector = biomt
        logger.debug('matrix_id: %s', matrix_id)
        assert matrix_id in ['1', '2', '3']
        if self._biomt is None:
            self._biomt = []
        # New row
        if matrix_id == '1':
            # New transformation
            if len(self._chains) > len(self._biomt):
                transformations = []
                self._biomt.append(transformations)
            else:
                transformations = self._biomt[-1]
            row = {'transformation_id': int(transformation_id)}
            transformations.append(row)
        else:
            transformations = self._biomt[-1]
            row = transformations[-1]
            self._validate_row(row, transformation_id, matrix_id)
        assert len(self._biomt) == len(self._chains)
        logger.debug('row: %s', row)
        row['_pdbx_struct_oper_list.matrix[{}][1]'.format(matrix_id)] = matrix_0
        row['_pdbx_struct_oper_list.matrix[{}][2]'.format(matrix_id)] = matrix_1
        row['_pdbx_struct_oper_list.matrix[{}][3]'.format(matrix_id)] = matrix_2
        row['_pdbx_struct_oper_list.vector[{}]'.format(matrix_id)] = vector

    def _validate_row(self, row, transformation_id, matrix_id):
        """Sanity check on the data that we already have in ``row``."""
        assert row['transformation_id'] == int(transformation_id)
        previous_matrix_data_is_present = [
            '_pdbx_struct_oper_list.matrix[{}][{}]'.format(k, i) in row
            for k in range(1, int(matrix_id)) for i in range(1, 3)
        ]
        previous_vector_data_is_present = [
            '_pdbx_struct_oper_list.vector[{}]'.format(k) in row for k in range(1, int(matrix_id))
        ]
        assert all(previous_matrix_data_is_present) and all(previous_vector_data_is_present)

    def _flush(self):
        logger.debug('_flush')
        assert str(self._biomolecule) not in self.bioassembly_data
        assert len(self._chains) == len(self._biomt)
        for chain_ids, transformations in zip(self._chains, self._biomt):
            transformations = [
                Transformation(t['transformation_id'], get_rotation(t), get_translation(t))
                for t in transformations
            ]
            for chain_id in chain_ids:
                self.bioassembly_data \
                    .setdefault(str(self._biomolecule), {}) \
                    .setdefault(chain_id, []) \
                    .extend(transformations)
        self._chains = None
        self._biomt = None


# === mmCIF ===
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


def get_label_id_to_auth_id_mapping(sdict):
    label_id_to_auth_id_df = \
        mmcif_key_to_dataframe(sdict, '_atom_site')[
            ['_atom_site.label_asym_id', '_atom_site.auth_asym_id']
        ].drop_duplicates()
    label_id_to_auth_id = dict(label_id_to_auth_id_df.values.tolist())
    if len(label_id_to_auth_id) != len(label_id_to_auth_id_df):
        raise Exception("Cound not reliably map 'label_asym_id' to 'auth_asym_id'!")
    return label_id_to_auth_id


def get_mmcif_bioassembly_data(sdict, use_auth_id=False):
    """Extract chain ids and transformations for each bioassembly from mmCIF data.

    See also
    --------
    :class:`ProcessRemark350`
    """
    # _pdbx_struct_assembly_gen
    df = mmcif_key_to_dataframe(sdict, '_pdbx_struct_assembly_gen') \
        .set_index('_pdbx_struct_assembly_gen.assembly_id')
    df = df[df['_pdbx_struct_assembly_gen.oper_expression'] != 'P']
    _pdbx_struct_assembly_gen = df

    # _pdbx_struct_oper_list
    df = mmcif_key_to_dataframe(sdict, '_pdbx_struct_oper_list') \
        .set_index('_pdbx_struct_oper_list.id')
    df = df[df.index.str.isdigit()]
    df['transformation_id'] = df.index.astype(int)
    df['rotation'] = df.apply(get_rotation, axis=1)
    df['translation'] = df.apply(get_translation, axis=1)
    _pdbx_struct_oper_list = df
    logger.debug("_pdbx_struct_oper_list: %s", _pdbx_struct_oper_list)

    if use_auth_id:
        label_id_to_auth_id = get_label_id_to_auth_id_mapping(sdict)

    bioassembly_data = OrderedDict()
    for bioassembly_id in _pdbx_struct_assembly_gen.index.drop_duplicates():
        logger.debug("bioassembly_id: %s", bioassembly_id)
        df = _pdbx_struct_assembly_gen.loc[bioassembly_id:bioassembly_id]
        for chain_ids, transformation_ids in zip(df['_pdbx_struct_assembly_gen.asym_id_list'],
                                                 df['_pdbx_struct_assembly_gen.oper_expression']):
            logger.debug("chain_ids: %s, transformation_ids: %s", chain_ids, transformation_ids)
            chain_ids = chain_ids.split(',')
            transformation_ids = _decode_transformation_ids(transformation_ids)
            logger.debug("(transformed) chain_ids: %s, transformation_ids: %s", chain_ids,
                         transformation_ids)
            transformations = \
                _pdbx_struct_oper_list \
                .loc[map(str, transformation_ids),
                     ['transformation_id', 'rotation', 'translation']] \
                .itertuples(index=False, name='Transformation')
            transformations = list(transformations)  # Required to materialize generator
            if use_auth_id:
                chain_ids = uniqueify([label_id_to_auth_id[c] for c in chain_ids])
            for chain_id in chain_ids:
                if use_auth_id:
                    existing_transformation_ids = set(
                        t.transformation_id
                        for t in bioassembly_data.get(str(bioassembly_id), {}).get(chain_id, []))
                else:
                    # There should be no duplicates in this case
                    existing_transformation_ids = set()
                new_transformations = [
                    t for t in transformations
                    if t.transformation_id not in existing_transformation_ids
                ]
                bioassembly_data \
                    .setdefault(str(bioassembly_id), OrderedDict()) \
                    .setdefault(chain_id, []) \
                    .extend(new_transformations)
    return bioassembly_data


def _decode_transformation_ids(transformation_ids):
    transformation_ids = transformation_ids.strip('()')
    try:
        return [int(transformation_ids)]
    except ValueError:
        pass
    try:
        return [int(v) for v in transformation_ids.split(',')]
    except ValueError:
        pass
    try:
        start, stop = transformation_ids.split('-')
        return list(range(int(start), int(stop) + 1))
    except ValueError:
        pass
    raise BioassemblyError("Could not parse transformation_ids: '{}'".format(transformation_ids))
