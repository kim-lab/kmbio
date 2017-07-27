"""Test doctests

Can't run doctests using ``pytest --doctest-modules`` because this adds the root path
 to `sys.path`.
 """
import doctest
import importlib
import logging
import os
import os.path as op
import pkgutil
import tempfile

import numpy as np
import pytest

import kmbio

logger = logging.getLogger(__name__)

DOCTEST_OPTIONFLAGS = (
    doctest.NORMALIZE_WHITESPACE | doctest.ELLIPSIS | doctest.IGNORE_EXCEPTION_DETAIL)

DOCTEST_EXTRAGLOBS = {'os': os, 'op': op, 'tempfile': tempfile, 'np': np}


def iter_submodules(package):
    """Import all submodules of a module, recursively, including subpackages.

    Adapted from https://stackoverflow.com/a/25562415/2063031
    """
    yield package.__name__, package
    for loader, name, ispkg in pkgutil.walk_packages(package.__path__):
        module = importlib.import_module(package.__name__ + '.' + name)
        if ispkg:
            yield from iter_submodules(module)
        else:
            yield module.__name__, module


@pytest.mark.skip(reason="TODO: Not working because of missing dependencies")
@pytest.mark.parametrize("module_name, module", iter_submodules(kmbio))
def test_doctest(module_name, module):
    failure_count, test_count = doctest.testmod(
        module, optionflags=DOCTEST_OPTIONFLAGS, extraglobs=DOCTEST_EXTRAGLOBS)
    assert failure_count == 0
