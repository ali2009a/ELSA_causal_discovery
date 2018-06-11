
import os
from os.path import basename
from constants import *
import numpy as np
currentPath = os.path.dirname(__file__)
print os.listdir(currentPath)
print currentPath
__all__ = [os.path.splitext(basename(f))[0] for f in os.listdir(currentPath)]
print "salammmmmmmmm"
print __all__