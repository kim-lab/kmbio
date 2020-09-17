import logging

import pytest

logger = logging.getLogger(__name__)


def parametrize(arg_string, arg_list):
    """
    Args:
        arg_string: Comma-separated string of arguments (e.g. 'pdb_id, pdb_type').
        arg_list: List of arguments or argument dictionaries.
    """
    logger.info("arg_string: %s", arg_string)
    logger.info("arg_list: %s", arg_list)
    if "," in arg_string:
        keys = arg_string.replace(" ", "").split(",")
        args = [tuple(r[k] for k in keys) for r in arg_list]
    else:
        key = arg_string.replace(" ", "")
        if isinstance(arg_list[0], dict) and key in arg_list[0]:
            args = [r[key] for r in arg_list]
        else:
            args = arg_list
    logger.info("args: %s", args)
    return pytest.mark.parametrize(arg_string, args)
