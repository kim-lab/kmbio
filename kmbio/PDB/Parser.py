from abc import ABC, abstractmethod


class Parser(ABC):

    @abstractmethod
    def get_structure(self, structure_id, filename):
        raise NotImplementedError

    @abstractmethod
    def get_biological_assembly(self, structure_id, filename, biological_assembly_index):
        structure = self.get_structure(structure_id, filename)
