#   nls/__init__.py
#   (c) Daniel Bershatsky, 2016
#   See LICENSE for details

from __future__ import print_function

from .animation import *
from .model import *
from .pumping import *

__all__ = [
    'animation',
    'model',
    'native'
    'pumping',
    'solver',
]

MAJOR = 0
MINOR = 1
PATCH = 0


def version():
    return  MAJOR, MINOR, PATCH