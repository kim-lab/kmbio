"""
Turn an mmCIF file into a dictionary.

"""
from cpython cimport bool


cdef list split(unicode line):

    cdef list tokens = []
    cdef bool is_double_quote = False
    cdef bool is_single_quote = False
    cdef unicode cur_token = u''
    cdef unicode whitespace = u' \t\r\n'
    cdef unicode DQ = u'"'
    cdef unicode SQ = u"'"
    cdef unicode c

    def flush():
        nonlocal tokens, cur_token
        if cur_token:
            tokens.append(cur_token)
            cur_token = u''

    for c in line:
        if is_double_quote:
            if c == DQ:
                flush()
                is_double_quote = False
            else:
                cur_token += c
        elif is_single_quote:
            if c == SQ:
                flush()
                is_single_quote = False
            else:
                cur_token += c
        elif c == DQ:
            flush()
            is_double_quote = True
        elif c == SQ:
            flush()
            is_single_quote = True
        elif c in whitespace:
            flush()
        else:
            cur_token += c
    flush()

    return tokens


def tokenize(handle):
    cdef unicode line
    cdef unicode token

    for line in handle:
        if line.startswith(u"#"):
            yield line.strip()
        elif line.startswith(u";"):
            # Multi-line string are enclosed in ';'
            token = line[1:].strip()
            for line in handle:
                line = line.strip()
                if line == u';':
                    break
                token += line
            yield token
        elif '"' in line or "'" in line:
            for token in split(line.strip()):
                yield token
        else:
            for token in line.strip().split():
                yield token


def process_tokens(tokens):
    cdef dict mmcif_dict = {}
    cdef list keys = []
    cdef unicode token = next(tokens)
    cdef bool loop_flag = False
    cdef unicode key = None

    mmcif_dict[token[0:5]] = token[5:]
    i = 0
    n = 0
    for token in tokens:
        # Parse loops
        if token == u"loop_":
            # Start a new loop block
            loop_flag = True
            keys = []
            i = 0  # Index of the current column
            n = 0  # Number of columns in this "DataFrame"
            continue
        elif loop_flag:
            if token == '#':
                loop_flag = False
                continue
            elif token.startswith(u"_") and i == 0:
                # Parse column names into `keys`
                mmcif_dict[token] = []
                keys.append(token)
                n += 1
                continue
            else:
                # Add elements one column at a time
                mmcif_dict[keys[i % n]].append(token)
                i += 1
                continue
        if token == '#':
            continue 
        # Key-value pairs
        if key is None:
            key = token
        else:
            mmcif_dict[key] = token
            key = None
    return mmcif_dict


def mmcif2dict(file):
    """Parse a mmCIF file and return a dictionary.

    Parameters
    ----------
    file : str
        Name of the mmCIF file OR an open filehandle.
    """
    cdef dict mmcif_dict
    cdef bool close_handle

    if isinstance(file, str):
        handle = open(file, 'rt')
        close_handle = True
    else:
        handle = file
        close_handle = False

    try:
        tokens = tokenize(handle)
        mmcif_dict = process_tokens(tokens)
    finally:
        if close_handle:
            handle.close()

    return mmcif_dict
