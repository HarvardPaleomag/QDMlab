import os
import admin_tools


def get_in_out_params(reload = False, save = False):
    """
    INTERNAL
    """
    func_files = {}
    code = admin_tools.preload_code(reload=reload)

    # iterate over all files/folders
    for root in code:
        for f in code[root]:
            if 'GPUfit_MATLAB' in root:
                continue
            # skip non matlab files
            if not f.endswith('.m') :
                continue
            if f in func_files:
                continue

            # skip all tests
            if 'test' in f:
                continue

            lines = code[root][f]

            for i,l in enumerate(lines):

                if l.lstrip().startswith('function'):
                    indent = detect_indent(l)

                l = l.lstrip()

                if l.startswith('function'):
                    fName = split_funcname(l)

                    if '=' in l:
                        out = l[9:].split('=')[0].replace(' ','').replace('[', '').replace(']', '').split(',')
                        out = [elem.split('%')[0] for elem in out] # remove comments
                    else:
                        out = ''

                    if fName in ['gpufit_cuda_available', 'gpufit_version','gpufit']:
                        continue

                    if '(' in l:
                        inputs = l[9:].split('(')[1].replace(' ','').replace(')','').replace('\n','').split(',')

                    func_files[fName] = {}
                    func_files[fName]['parameters'] = inputs
                    func_files[fName]['returns'] = out
                    func_files[fName]['funcPath'] = os.path.join(root, f)
                    func_files[fName]['line'] = l
                    func_files[fName]['indent'] = indent
                    if i > 1:
                        func_files[fName]['internal'] = split_funcname(lines[0])
                if l.startswith('arguments'):
                    save = True
                    continue
                if l.startswith('end'):
                    save = False
                    continue
                if l.startswith('%'):
                    continue
                if save:
                    # arguments
                    arg = l.replace("\n", "").split(' ')
                    arg = [n for n in arg if n]
                    if not arg:
                        continue

                    arg = arg[0].split('.')

                    if len(arg) == 1:
                        arg = arg[0]
                        kwargs = False
                    else:
                        kwarg = arg[0]
                        arg = arg[1]
                        # remove kwargs name from args
                        if kwarg in func_files[fName]['parameters']:
                            func_files[fName]['parameters'].remove(kwarg)

                    # add kwargs and default values for kwargs
                    if '=' in l:
                        if not 'kwargs' in func_files[fName]:
                            func_files[fName]['kwargs'] = {}
                        default = l.split('=')[1].replace("\n", "").strip().replace(';','')
                        func_files[fName]['kwargs'][arg] = default

    return func_files

def split_funcname(l):
    """
    INTERNAL
    """
    if '=' in l:
        fName = l[9:].split('=')[1].split('(')[0].strip(' ')
    else:
        fName = l[9:].split('(')[0].strip(' ')
    return fName

def detect_indent(l):
    """
    INTERNAL
    """
    out = ''
    for i,c in enumerate(l):
        if c in [' ', '\t']:
            out += c
        else:
            break
    return out

def first_line_comments(reload = False, save = False):
    """
    Function to automatically add the first line comment. Adds all possible parameters to it.
    """
    func_files = get_in_out_params(reload = reload)

    for func in sorted(func_files):
        if 'GPUfit_MATLAB' in func_files[func]['funcPath']:
            continue

        p = ', '.join([k for k in func_files[func]['parameters'] if not k in ['kwargs', 'filter', 'cline_idx', 'funcPath']])
        if 'kwargs' in func_files[func]:
            args_list = [k for k in func_files[func]['kwargs'].keys()]
        else:
            args_list = []

        line = func_files[func]['indent']+'%'
        if func_files[func]['returns']:
            line += '[' + ', '.join(func_files[func]['returns']) + '] = '
        line += f'{func}('
        if p:
            line += f'{p}'
            if args_list:
                line += '; '
        if args_list:
            line += '\''+ '\', \''.join(args_list) + '\''
        line += ')\n'

        with open(func_files[func]['funcPath'], 'r+', encoding='latin-1') as f:
            lines = f.readlines()
            for i, l in enumerate(lines):
                if func_files[func]['line'] in l:
                    print(func_files[func]['funcPath'])
                    print('-'*100)
                    write = True
                    if func not in lines[i+1]:
                        l0 = line.replace('\n', '')
                        print(f'ADDING:  {l0}')
                        l1 = lines[i].replace('\n', '\\n')
                        print(f'BETWEEN: {l1}')
                        l2 = lines[i+1].replace('\n', '\\n')
                        print(f'AND:     {l2}')
                        lines.insert(i+1,line)  # inserts line
                    else:
                        if lines[i+1].rstrip() == line.rstrip():
                            print('NO CHANGE')
                            write = False
                        else:
                            lold = lines[i+1].replace('\n', '\\n')
                            print(f'REPLACING: {lold}')
                            lnew = line.replace('\n', '\\n')
                            print(f'WITH:      {lnew}')
                            lines[i+1] = line
                    print('-'*100)
                    print()
                    break
            if save:
                if write:
                    print(func)
                    f.truncate(0)         # truncates the file
                    f.seek(0)             # moves the pointer to the start of the file
                    f.writelines(lines)   # write the new data to the file
