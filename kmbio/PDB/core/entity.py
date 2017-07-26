# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
Base class for Residue, Chain, Model and Structure classes.

It is a simple container class, with list and dictionary like properties.
"""
from abc import abstractmethod
import logging
from copy import copy

from kmbio.PDB.exceptions import PDBConstructionException

logger = logging.getLogger(__name__)


class Index:
    def __init__(self, entity):
        self._entity = entity

    def __getitem__(self, idx):
        """Return the child with given id."""
        if not isinstance(idx, (int, slice)):
            idx = slice(idx.start, idx.stop)
        return self._entity._child_list[idx]

    def __delitem__(self, idx):
        """Remove a child."""
        items = self[idx]
        for item in items if isinstance(items, list) else (items, ):
            del self._entity[item.id]


class Entity(object):
    """
    Basic container object. Structure, Model, Chain and Residue
    are subclasses of Entity. It deals with storage and lookup.
    """

    def __init__(self, id, children=None):
        self._id = id
        self._full_id = None
        self.parent = None
        self._child_list = []
        self._child_dict = {}
        # Dictionary that keeps additional properties
        self.xtra = {}
        #
        self.ix = Index(self)
        if children is not None:
            self.add(children)

    # Special methods

    def __getitem__(self, id):
        """Return the child with given id."""
        return self._child_dict[id]

    def __setitem__(self, id, item):
        """Add a child."""
        assert id == item.id
        self.add([item])

    def __delitem__(self, id):
        """Remove a child."""
        child = self._child_dict.pop(id)
        self._child_list.remove(child)
        child.parent = None

    def __contains__(self, id):
        """True if there is a child element with the given id."""
        return id in self._child_dict

    def __iter__(self):
        """Iterate over all children."""
        for child in self._child_list:
            yield child.id

    def __len__(self):
        """Return the number of children."""
        return len(self._child_list)

    def values(self):
        """Overrides the `values` method from `MutableMapping` for speed."""
        return self._child_list

    # Private methods

    def _reset_full_id(self):
        """Reset the full_id.

        Sets the full_id of this entity and
        recursively of all its children to None.
        This means that it will be newly generated
        at the next call to get_full_id.
        """
        for child in self._child_list:
            try:
                child._reset_full_id()
            except AttributeError:
                pass  # Atoms do not cache their full ids.
        self._full_id = None

    # Public methods

    @property
    def id(self):
        return self._id

    @id.setter
    def id(self, new_id):
        """Change the id of this entity.

        This will update the child_dict of this entity's parent
        and invalidate all cached full ids involving this entity.

        @raises: ValueError
        """
        if self.parent:
            if new_id in self.parent:
                raise ValueError(
                    "Cannot change id from `{0}` to `{1}`. "
                    "The id `{1}` is already used for a sibling of this entity."
                    .format(self._id, new_id))
            del self.parent._child_dict[self._id]
            self.parent._child_dict[new_id] = self
        self._id = new_id
        self._reset_full_id()

    @property
    @abstractmethod
    def level(self):
        """Return level in hierarchy.

        A - atom
        R - residue
        C - chain
        M - model
        S - structure
        """
        raise NotImplementedError

    def pop(self, id):
        """Remove and return a child."""
        child = self._child_dict.pop(id)
        self._child_list.remove(child)
        child.parent = None
        return child

    def clear(self):
        for child in self.values():
            child.parent = None
        self._child_list.clear()
        self._child_dict.clear()
        self.xtra.clear()

    def add(self, entities):
        """Add a child to the Entity."""
        # Single entity
        if entities is None:
            logger.info("Trying to add a 'None' child to {}".format(self))
            return
        if isinstance(entities, (Entity, DisorderedEntityWrapper)):
            entities = [entities]
        elif not isinstance(entities, list):
            # Like a generator...
            entities = list(entities)
        if any(c.id in self for c in entities):
            raise PDBConstructionException(
                "Some of the entities are defined twice")
        if len({c.id for c in entities}) < len(entities):
            raise PDBConstructionException(
                "Some of the entities are duplicates")
        for entity in entities:
            entity.parent = self
        self._child_list.extend(entities)
        self._child_dict.update({c.id: c for c in entities})

    def insert(self, pos, entities):
        """Add a child to the Entity at a specified position."""
        # Single entity
        if isinstance(entities, (Entity, DisorderedEntityWrapper)):
            entities = [entities]
        elif not isinstance(entities, list):
            # Like a generator...
            entities = list(entities)
        if any(c.id in self for c in entities):
            raise PDBConstructionException(
                "Some of the entities are defined twice")
        if len({c.id for c in entities}) < len(entities):
            raise PDBConstructionException(
                "Some of the entities are duplicates")
        for entity in entities:
            entity.parent = self
        self._child_list[pos:pos] = entities
        self._child_dict.update({c.id: c for c in entities})

    @property
    def full_id(self):
        """Return the full id.

        The full id is a tuple containing all id's starting from
        the top object (Structure) down to the current object. A full id for
        a Residue object e.g. is something like:

        ("1abc", 0, "A", (" ", 10, "A"))

        This corresponds to:

        Structure with id "1abc"
        Model with id 0
        Chain with id "A"
        Residue with id (" ", 10, "A")

        The Residue id indicates that the residue is not a hetero-residue
        (or a water) because it has a blank hetero field, that its sequence
        identifier is 10 and its insertion code "A".
        """
        if self._full_id is None:
            entity_id = self.id
            lst = [entity_id]
            parent = self.parent
            while parent is not None:
                entity_id = parent.id
                lst.append(entity_id)
                parent = parent.parent
            lst.reverse()
            self._full_id = tuple(lst)
        return self._full_id

    def transform(self, rot, tran):
        """
        Apply rotation and translation to the atomic coordinates.

        Example:
                >>> rotation=rotmat(pi, Vector(1, 0, 0))
                >>> translation=array((0, 0, 1), 'f')
                >>> entity.transform(rotation, translation)

        @param rot: A right multiplying rotation matrix
        @type rot: 3x3 Numeric array

        @param tran: the translation vector
        @type tran: size 3 Numeric array
        """
        for o in self.values():
            o.transform(rot, tran)

    def copy(self):
        shallow = copy(self)  # copy class type, etc.
        # Need a generator from self because lazy evaluation:
        Entity.__init__(shallow, shallow.id, (c.copy() for c in self.values()))
        shallow.xtra = self.xtra.copy()
        return shallow


class DisorderedEntityWrapper(object):
    """
    This class is a simple wrapper class that groups a number of equivalent
    Entities and forwards all method calls to one of them (the currently selected
    object). DisorderedResidue and DisorderedAtom are subclasses of this class.

    E.g.: A DisorderedAtom object contains a number of Atom objects,
    where each Atom object represents a specific position of a disordered
    atom in the structure.
    """

    def __init__(self, id):
        self.id = id
        self._child_dict = {}  # this is more of a sibling dict
        self.selected_child = None
        self._parent = None
        self.disordered = 2

    # Special methods

    def __getattr__(self, method):
        """Forward the method call to the selected child."""
        if not hasattr(self, 'selected_child'):
            # Avoid problems with pickling
            # Unpickling goes into infinite loop!
            raise AttributeError
        return getattr(self.selected_child, method)

    def __getitem__(self, id):
        """Return the child with the given id."""
        return self.selected_child[id]

    # XXX Why doesn't this forward to selected_child?
    # (NB: setitem was here before getitem, iter, len, sub)
    def __setitem__(self, id, child):
        """Add a child, associated with a certain id."""
        self._child_dict[id] = child

    def __contains__(self, id):
        """True if the child has the given id."""
        return (id in self.selected_child)

    def __iter__(self):
        """Iterate over the children."""
        for item in self.selected_child:
            yield item

    def __len__(self):
        """Return the number of children."""
        return len(self.selected_child)

    def __sub__(self, other):
        """Subtraction with another object."""
        return self.selected_child - other

    # Public methods

    @property
    def parent(self):
        """Return parent."""
        return self._parent

    @parent.setter
    def parent(self, parent):
        if parent is None:
            # Detach parent
            self.parent = None
            for child in self.disordered_get_list():
                child.parent = None
        else:
            self._parent = parent

    @parent.setter
    def parent(self, parent):
        """Set the parent for the object and its children."""
        self._parent = parent
        for child in self.disordered_get_list():
            child.parent = parent

    def disordered_has_id(self, id):
        """True if there is an object present associated with this id."""
        return (id in self._child_dict)

    def disordered_select(self, id):
        """Select the object with given id as the currently active object.

        Uncaught method calls are forwarded to the selected child object.
        """
        self.selected_child = self._child_dict[id]

    def disordered_add(self, child):
        """This is implemented by DisorderedAtom and DisorderedResidue."""
        raise NotImplementedError

    def disordered_get_id_list(self):
        """Return a list of id's."""
        # sort id list alphabetically
        return sorted(self._child_dict)

    def disordered_get(self, id=None):
        """Get the child object associated with id.

        If id is None, the currently selected child is returned.
        """
        if id is None:
            return self.selected_child
        return self._child_dict[id]

    def disordered_get_list(self):
        """Return list of children."""
        return list(self._child_dict.values())
