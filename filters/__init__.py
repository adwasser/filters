import os
_root = os.path.abspath(os.path.dirname(__file__)) + '/'

from .version import __version__
from .filters import Filter, FilterSet


