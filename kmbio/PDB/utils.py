from kmbio.PDB import Atom


def allequal(s1, s2):
    if type(s1) != type(s2):
        raise Exception
    if isinstance(s1, Atom):
        return s1 == s2
    equal = (
        len(s1) == len(s2) and
        all(allequal(so1, so2) for (so1, so2) in zip(s1.values(), s2.values()))
    )
    return equal
