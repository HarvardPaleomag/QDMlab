import os, sys
path = os.path.dirname(__file__)
path = os.path.abspath(os.path.join(os.path.dirname(path), '..'))

sys.path.append(path)

import admin_tools
from admin_tools import arguments

def main():
	arguments.first_line_comments(save=True)

if __name__ == '__main__':
	main()