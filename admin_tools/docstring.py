import os
import admin_tools
from admin_tools import arguments
import numpy as np


def get_docstring_indices(reload=False, debug_fname=''):
    """
    Walks through all functions and gets the index for the docstrings
    """
    code = admin_tools.preload_code(reload=reload)
    params = arguments.get_in_out_params(reload=reload)
    comment_block_index = {}
    for f in params:
        for package, files in code.items():
            for file, c in files.items():
                strt_idx = [i for i, l in enumerate(c) if l.startswith('function') and f' {f}(' in l.replace('=', ' ')]
                if len(strt_idx) != 1:
                    # file does not include the function
                    continue
                else:
                    strt_idx = strt_idx[0]
                comment_idx = []
                idx = np.array([i for i, l in enumerate(c[strt_idx:]) if l.rstrip().startswith('%')])
                df = np.diff(idx, prepend=True)

                try:
                    until_idx = np.where(df > 1)[0][0]
                except IndexError:
                    until_idx = len(idx)

                comment_idx = idx[:until_idx]

                comment_block_index[f] = comment_idx

                if f == debug_fname:
                    print(f, file)
                    print(strt_idx, c[strt_idx])
                    print('comment indices: ', comment_idx)
                    print()
                    print('extracted comments')
                    print('------------------')
                    for l in np.array(c)[comment_idx]:
                        print(l, end='')

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
                strt_idx = [i for i, l in enumerate(c) if l.startswith('function') and f' {f}(' in l.replace('=', ' ')]
                if len(strt_idx) != 1:
                    # file does not include the function
                    continue
                else:
                    docstring[f] = [c[i] for i in idx[1:]]
    return docstring


def make_params_doc(p):
    lines = ['', 'Parameters', '----------']
    for i in p['parameters']:
        lines.append(f'  {i}:')
    for k, v in p['kwargs'].items():
        lines.append(f'  {k}: ({v})')

    lines.extend(['', 'Returns', '----------'])
    for i in p['returns']:
        lines.append(f'  {i}:')

    lines = [f'% {i}' for i in lines]

    for l in lines:
        print(l)

    # return lines
    # return ''.join(lines)
    # return '\n\r'.join(lines)