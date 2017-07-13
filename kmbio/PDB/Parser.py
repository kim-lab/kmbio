from abc import ABC, abstractmethod


class Parser(ABC):

    @abstractmethod
    def get_structure(self, filename, structure_id=None, biological_assembly_id=None):
        raise NotImplementedError
