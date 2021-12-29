import os, sys
path = os.path.dirname(__file__)

sys.path.append(os.path.split(path)[0])

import admin_tools
from admin_tools import arguments

def main():
	arguments.first_line_comments(save=True)

if __name__ == '__main__':
	main()