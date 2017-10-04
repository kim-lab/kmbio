# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""The structure class, representing a macromolecular structure."""
import numpy as np
import pandas as pd

from .entity import Entity


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

    def extract_models(self, model_ids):
        # TODO: Not sure if this is neccessary
        structure = Structure(self.id)
        for model_id in model_ids:
            structure.add(self[model_id].copy())
        return structure

    def select(self, models=None, chains=None, residues=None, hetatms=None):
        """This method allows you to select things from structures using a variety of queries.

        In particular, you should be able to select one or more chains,
        and all HETATMs that are within a certain distance of those chains.
        """
        raise NotImplementedError

    def to_dataframe(self) -> (pd.DataFrame, np.ndarray):
        """Convert this structure into a pandas DataFrame

        The output of this method is intended to be compatible
        with :meth:`mdtraj.Topology.from_dataframe`.

        Returns
        -------
        atoms : :class:`pandas.DataFrame`
            The atoms in the structure, represented as a data frame.
        bonds : :class:`numpy.ndarray`
            The bonds in this structure, represented as an n_bonds x 2 array
            of the indices of the atoms involved in each bond.
        """
        raise NotImplementedError

    @staticmethod
    def from_dataframe(atoms, bonds=None):
        """Generate a new structure from a dataframe of atoms and an array of bonds."""
        raise NotImplementedError

    @property
    def models(self):
        for m in self:
            yield m

    @property
    def chains(self):
        for m in self:
            for c in m:
                yield c

    @property
    def residues(self):
        for c in self.chains:
            for r in c:
                yield r

    @property
    def atoms(self):
        for r in self.residues:
            for a in r:
                yield a
