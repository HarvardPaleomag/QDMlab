import os
from os.path import dirname, abspath
from admin_tools.arguments import get_in_out_params, first_line_comments
from admin_tools.dependencies import called_in, get_dependencies, preload_code

install_directory = dirname(dirname(abspath(__file__)))
QDMlab_code = None
