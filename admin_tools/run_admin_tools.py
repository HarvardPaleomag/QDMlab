import os, sys
path = os.path.dirname(__file__)
path = os.path.abspath(os.path.join(os.path.dirname(path), '..'))

sys.path.append(path)

import admin_tools
from admin_tools import arguments, dependencies

def main():
	arguments.first_line_comments(save=True)
	dependencies.get_dependencies(save = True, reload = True, debug = True)

if __name__ == '__main__':
	main()