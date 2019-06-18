from constants import *
import numpy as np

def harmonize(value):
	if ( (value == REFUSAL) or (value == DONT_KNOW) or (value== SCHD_NOT_APPLICABLE)):
		return np.nan
	elif (value == NOT_APPLICABLE):
		return 0
	else:
		return value

def binarize(value):
	if (value==0):
		return 0
	elif (value>3):
		return 1
	else:
		return 2
