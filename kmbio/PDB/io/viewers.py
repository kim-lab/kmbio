import io

import kmbio.PDB
from kmbio.PDB import Structure


def structure_to_ngl(structure: Structure):
    import nglview

    class NGLStructure(nglview.Structure):
        def __init__(self, entity, ext="pdb", params={}) -> None:
            super().__init__()
            self.path = ""
            self.ext = ext
            self.params = params
            self._entity = entity

        def get_structure_string(self) -> str:
            io_str = io.StringIO()
            kmbio.PDB.save(self._entity, io_str)
            return io_str.getvalue()

    return NGLStructure(structure)


def view_structure(structure: Structure, **kwargs):
    import nglview

    structure_ngl = structure_to_ngl(structure)
    return nglview.NGLWidget(structure_ngl, **kwargs)
