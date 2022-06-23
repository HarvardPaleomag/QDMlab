import pandas as pd
import os
from pandas import ExcelWriter
import admin_tools

def preload_code(reload = False, debug = False):
    if admin_tools.QDMlab_code is not None and not reload:
        return admin_tools.QDMlab_code

    QDMlab_code = {}
    for root, dirs, files in os.walk(admin_tools.install_directory):
        if 'GPUfit_MATLAB' in root:
            continue
        for file in files:
            if not file.endswith('.m'):
                continue
            if debug:
                print(f'reading {file}')
            with open(os.path.join(root, file), encoding='latin-1') as f:
                fdata = f.readlines()
                if not root in QDMlab_code:
                    QDMlab_code[root] = {}

                QDMlab_code[root][file] = fdata
    admin_tools.QDMlab_code = QDMlab_code
    return QDMlab_code

def called_in(save = False, reload = False, debug = False):

    code = pd.DataFrame()
    QDMlab_code = preload_code(reload = reload, debug = debug)

    for root, files in QDMlab_code.items():
        for file, fdata in files.items():
            for line in fdata:
                if line.startswith('function'):
                    funcName = line
                    if '=' in funcName:
                        funcName = funcName.split('=')[1]
                    funcName = funcName.split('(')[0]
                    funcName = funcName.replace('function', '').rstrip().lstrip()
                    code.loc[funcName,'root'] = root
                    code.loc[funcName,'file'] = file

    for funcName in code.index:
        n = 1
        for root, files in QDMlab_code.items():
            for file, fdata in files.items():
                    if f'{funcName}(' in ''.join(fdata):
                        code.loc[funcName, n] = file
                        n += 1
                        break


    # for file in code:
    #     if not 'called in' in code[file]:
    #          code[file]['called in'] = 'none'
    if save:
        with ExcelWriter(os.path.join(r'C:\Users\micha\OneDrive\Desktop', 'dependencies.xlsx'), engine="openpyxl") as writer:
            code.to_excel(writer, sheet_name='called in')

    return code

def get_dependencies(save = False, reload = False, debug = False):
    deps = pd.DataFrame()
    calledIn = called_in()

    for item in [i for i in set(calledIn.iloc[:,1:].values.flatten()) if not 'nan' in str(i)]:
        n = 1
        for r, row in calledIn.iterrows():
            row = row.values
            if item in row:
                deps.loc[item,n] = r
                n += 1
                continue

    if save:
        with ExcelWriter(os.path.join(r'C:\Users\micha\OneDrive\Desktop', 'dependencies.xlsx'), engine="openpyxl", mode='r+') as writer:
            deps.to_excel(writer, sheet_name='dependencies')

    return deps

if __name__ == '__main__':
    calls = called_in(save=True)
    deps = get_dependencies(save=True)
    print(deps)