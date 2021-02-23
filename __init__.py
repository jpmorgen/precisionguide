"""Module to enable precise telescope pointing"""

import sys

if sys.platform == 'win32':
    from win32.py import *

from pgdata import *

