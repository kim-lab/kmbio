import abc
import io

from .savers import save


def view(entity, **kwargs):
    """Veiw structure (or another entity) inside an NGLViewer."""
    return NGLViewer(entity).get_widget(**kwargs)


class Viewer(abc.ABC):

    @abc.abstractmethod
    @property
    def _structure(self):
        raise NotImplementedError

    @abc.abstractmethod
    def get_widget(self, **kwargs):
        raise NotImplementedError


class NGLViewer(Viewer):

    def __init__(self, entity):
        self.entity = entity

    @property
    def _structure(self):
        from nglview import Structure

        class NGLStructure(Structure):
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

        return NGLStructure(self.entity)

    def get_widget(self, **kwargs):
        from nglview import NGLWidget
        return NGLWidget(self._structure, **kwargs)
