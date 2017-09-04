from collections import OrderedDict

import pytest

from kmbio.PDB.utils import sort_ordered_dict


@pytest.mark.parametrize(
    "unsorted_dict, sorted_dict",
    [
        (OrderedDict([(1, 10), (3, 30), (2, 20)]), OrderedDict([(1, 10), (2, 20), (3, 30)])),
        (OrderedDict([('b', 'bye bye'), ('a', 'hello')]),
         OrderedDict([('a', 'hello'), ('b', 'bye bye')])),
    ])
def test_sort_ordered_dict(unsorted_dict, sorted_dict):
    assert unsorted_dict != sorted_dict
    print(unsorted_dict, len(sorted_dict))
    sort_ordered_dict(unsorted_dict)
    assert unsorted_dict == sorted_dict
