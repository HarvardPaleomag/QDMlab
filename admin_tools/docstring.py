import os
import admin_tools
from admin_tools import arguments

def get_docstring_indices(reload = False):
    """
    Walks through all functions and gets the index for the docstrings
    """
    code = admin_tools.preload_code(reload = reload)
    params = arguments.get_in_out_params(reload = reload)
    comment_block_index = {}
    for f in params:
        for package, files in code.items():
            for file, c in files.items():
                strt_idx = [i for i, l in enumerate(c) if l.startswith('function') and f in l]
                if len(strt_idx) != 1:
                    # file does not include the function
                    continue
                else:
                    strt_idx = strt_idx[0]
                comment_idx = []
                for i, l in enumerate(c[strt_idx+1:]):
                    if '%' in l:
                        comment_idx.append(i)
                    else:
                        break
                comment_block_index[f] = comment_idx
                break
    return comment_block_index


def get_docstring(reload = False):
    """
    walks through all functions and gets the docstring as a list
    """
    code = admin_tools.preload_code(reload = reload)
    params = arguments.get_in_out_params(reload = reload)
    comment_block_index = get_docstring_indices(reload = reload)
    docstring = {}
    for f, idx in comment_block_index.items():
        for package, files in code.items():
            for file, c in files.items():
                strt_idx = [i for i, l in enumerate(c) if l.startswith('function') and f in l]
                if len(strt_idx) != 1:
                    # file does not include the function
                    continue
                else:
                    docstring[f] = [c[i] for i in idx[1:]]
    return docstring