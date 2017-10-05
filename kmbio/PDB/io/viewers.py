import io

import nglview

from . import save
from .. import Entity


class KMBioStructure(nglview.Structure):

    def __init__(self, entity, ext='pdb', params={}):
        super().__init__()
        self.path = ''
        self.ext = ext
        self.params = params
        self._entity = entity

    def get_structure_string(self):
        io_str = io.StringIO()
        save(self._entity, io_str)
        return io_str.getvalue()


def show(entity, **kwargs):
    """Veiw structure (or another entity) inside an NGLViewer."""
    return nglview.NGLWidget(KMBioStructure(entity), **kwargs)


Entity.show = lambda self, **kwargs: show(self, **kwargs)
Entity.to_ngl = lambda self: KMBioStructure(self)
