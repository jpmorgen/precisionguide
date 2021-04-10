"""Module to enable precise telescope pointing"""

import sys

from .pgdata import *

if sys.platform == 'win32':
    from .win32 import *

