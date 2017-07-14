# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""The structure class, representing a macromolecular structure."""
import numpy as np
import pandas as pd

from kmbio.PDB.Entity import Entity


class Structure(Entity):
    """
    The Structure class contains a collection of Model instances.
    """
    level = "S"

    def __repr__(self):
        return "<Structure id=%s>" % self.id

    def __lt__(self, other):
        return self.id.lower() < other.id.lower()

    def __le__(self, other):
        return self.id.lower() <= other.id.lower()

    def __eq__(self, other):
        return self.id.lower() == other.id.lower()

    def __ne__(self, other):
        return self.id.lower() != other.id.lower()

    def __ge__(self, other):
        return self.id.lower() >= other.id.lower()

    def __gt__(self, other):
        return self.id.lower() > other.id.lower()

    def get_models(self):
        for m in self.values():
            yield m

    def get_chains(self):
        for m in self.get_models():
            for c in m.values():
                yield c

    def get_residues(self):
        for c in self.get_chains():
            for r in c.values():
                yield r

    def get_atoms(self):
        for r in self.get_residues():
            for a in r.values():
                yield a

    def extract_models(self, model_ids):
        structure = Structure(self.id)
        for model_id in model_ids:
            structure.add(self[model_id].copy())
        return structure

    def to_dataframe(self) -> (pd.DataFrame, np.ndarray):
        """Convert this structure into a pandas DataFrame

        The output of this method is intended to be compatible
        with the :mdtraj:`Topology.from_dataframe` method.

        Returns
        -------
        atoms : pandas.DataFrame
            The atoms in the structure, represented as a data frame.
        bonds : np.ndarray
            The bonds in this structure, represented as an n_bonds x 2 array
            of the indices of the atoms involved in each bond.
        """
        raise NotImplementedError