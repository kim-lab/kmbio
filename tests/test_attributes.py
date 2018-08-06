import pytest

import kmbio


@pytest.mark.parametrize("attribute", ["__version__"])
def test_attribute(attribute):
    assert getattr(kmbio, attribute)


def test_main():
    import kmbio

    assert kmbio
