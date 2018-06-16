import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
from munkres import Munkres
import scipy.stats
import re
import time
import subprocess
import os
from features import *


trtmntVar = set(["scfrdm","heacta", "heactb","heactc", "scorg03","scorg06","scorg05","scorg07","scako","heskb"])
confoundersVar = set(["scfrda","scfrdg","indager", "hehelf","dhsex","totwq10_bu_s"])
targetVar = set(["memtotb"])
allVar = trtmntVar|confoundersVar|targetVar


basePath = "/home/ali/Downloads/UKDA-5050-stata (2)/stata/stata13_se"
REFUSAL=-9
DONT_KNOW=-8
NOT_APPLICABLE=-1
SCHD_NOT_APPLICABLE=-2

NOT_ASKED=-3

NOT_IMPUTED = -999.0
NON_SAMPLE = -998.0
INST_RESPONDENT=  -995.0


def report(df, var):
	for i in [3,4,5]:
		print "wave",i
		print "min", df["{}_{}".format(var, i)].min()
		print "max", df["{}_{}".format(var, i)].max()
		print "mean", df["{}_{}".format(var, i)].mean()
		print "std", df["{}_{}".format(var, i)].std()

def harmonizeData(df):
	print allVar
	for var in allVar:
		df[var] = df[var].apply((globals()[var].harmonize))
	
	# for var in ["heacta", "heactb", "heactc", "scako", "heskb"]:
	#     df[var] = df[var].apply((globals()[var].binarize))

	return df



def binarizeData(df):
	pattern = r"[a-zA-Z0-9_]*_n$"
	cols = list(df.columns)
	cols.remove('idauniq')
	for var in cols:
		if not re.match(pattern, var):  
		    col_bin = var + '_b'
		    df[col_bin] = df[var].apply((globals()[var].binarize))
	return df


def normalizeData(df):
	cols = list(df.columns)
	cols.remove('idauniq')
	for col in cols:
	    col_norm = col + '_n'
	    df[col_norm] = (df[col] - df[col].min())/(df[col].max()- df[col].min())
	return df





def readWave2Data(basePath):
	waveNumber=2
	w3Core = pd.read_stata("{}/wave_{}_elsa_data.dta".format(basePath, waveNumber),convert_categoricals=False)
	w3Drv =  pd.read_stata("{}/wave_{}_ifs_derived_variables.dta".format(basePath, waveNumber),convert_categoricals=False)
	w3FinDrv = pd.read_stata('{}/wave_{}_financial_derived_variables.dta'.format(basePath, waveNumber),convert_categoricals=False)
	

	s1 = pd.merge(w3Core, w3Drv, how='inner', on=['idauniq'])
	df = pd.merge(s1, w3FinDrv, how='inner', on=['idauniq'])

	col_list = ["idauniq","heacta","heactb","heactc", "scorg03", "scorg06", "scorg05", "scorg07", "hegenh",
				 "scfrda" , "scfrdg","scako", "heskb", "indager", "dhsex" , "scfrdm", "memtotb","totwq10_bu_s" ]
	df = df [col_list] 

	df = df.rename(columns = {'hegenh':'hehelf'})
	df = harmonizeData(df)
	df = normalizeData(df)
	df = binarizeData(df)
	# df = df.ix[0:50,:]
	return df


def readWave3Data(basePath):
	waveNumber=3
	w3Core = pd.read_stata("{}/wave_{}_elsa_data.dta".format(basePath, waveNumber),convert_categoricals=False)
	w3Drv =  pd.read_stata("{}/wave_{}_ifs_derived_variables.dta".format(basePath, waveNumber),convert_categoricals=False)
	w3FinDrv = pd.read_stata('{}/wave_{}_financial_derived_variables.dta'.format(basePath, waveNumber),convert_categoricals=False)
	

	s1 = pd.merge(w3Core, w3Drv, how='inner', on=['idauniq'])
	df = pd.merge(s1, w3FinDrv, how='inner', on=['idauniq'])

	col_list = ["idauniq","heacta","heactb","heactc", "scorg03", "scorg06", "scorg05", "scorg07", "hegenh",
				 "scfrda" , "scfrdg","scako", "heskb", "indager", "dhsex" , "scfrdm", "memtotb","totwq10_bu_s" ]
	df = df [col_list] 

	df = df.rename(columns = {'hegenh':'hehelf'})
	df = harmonizeData(df)
	df = normalizeData(df)
	df = binarizeData(df)
	# df = df.ix[0:50,:]
	return df



def readWave4Data(basePath):
	waveNumber=4
	w3Core = pd.read_stata("{}/wave_{}_elsa_data.dta".format(basePath, waveNumber),convert_categoricals=False)
	w3Drv =  pd.read_stata("{}/wave_{}_ifs_derived_variables.dta".format(basePath, waveNumber),convert_categoricals=False)
	w3FinDrv = pd.read_stata('{}/wave_{}_financial_derived_variables.dta'.format(basePath, waveNumber),convert_categoricals=False)
	
	s1 = pd.merge(w3Core, w3Drv, how='inner', on=['idauniq'])
	df = pd.merge(s1, w3FinDrv, how='inner', on=['idauniq'])

	col_list = ["idauniq","heacta","heactb","heactc", "scorg03", "scorg06", "scorg05", "scorg07", "hehelf",
				 "scfrda" , "scfrdg","scako", "heskb", "indager", "dhsex" , "scfrdm", "memtotb","totwq10_bu_s" ]
	df = df [col_list] 

	df = harmonizeData(df)
	df = normalizeData(df)
	df = binarizeData(df)
	# df = df.ix[0:50,:]
	return df


def readWave5Data(basePath):
	waveNumber=5
	w3Core = pd.read_stata("{}/wave_{}_elsa_data.dta".format(basePath, waveNumber),convert_categoricals=False)
	w3Drv =  pd.read_stata("{}/wave_{}_ifs_derived_variables.dta".format(basePath, waveNumber),convert_categoricals=False)
	w3FinDrv = pd.read_stata('{}/wave_{}_financial_derived_variables.dta'.format(basePath, waveNumber),convert_categoricals=False)
	
	s1 = pd.merge(w3Core, w3Drv, how='inner', on=['idauniq'])
	df = pd.merge(s1, w3FinDrv, how='inner', on=['idauniq'])

	col_list = ["idauniq","heacta","heactb","heactc", "scorg03", "scorg06", "scorg05", "scorg07", "hehelf",
				 "scfrda" , "scfrdg","scako", "heskb", "indager", "dhsex" , "scfrdm", "memtotb","totwq10_bu_s" ]
	df = df [col_list] 

	df = harmonizeData(df)
	df = normalizeData(df)
	df = binarizeData(df)
	# df = df.ix[0:50,:]
	return df


def readData():
	df3 = readWave3Data(basePath)
	df4 = readWave4Data(basePath)
	df5 = readWave5Data(basePath)
	df34 = pd.merge(df3, df4, how='inner', on=['idauniq'],suffixes=('_3', ''))
	df345 = pd.merge(df34, df5, how='inner', on=['idauniq'],suffixes=('_4', '_5'))

	return df345


def computeMemIndexChange(row, waveNumber):
	memtotVarCur = "memtotb_{}".format(waveNumber) 
	memtotVarPrev = "memtotb_{}".format(waveNumber-1)
	return row[memtotVarCur] - row[memtotVarPrev]





def computeDistance(row1,row2):
	diff  = row1 - row2
	diff = diff[~np.isnan(diff)]
	return np.linalg.norm(diff)


def preProcessData(df):
	# df= df.dropna(axis=0, how="any")
	df= df.dropna(subset=["memtotb_3","memtotb_4","memtotb_5"])
	df["memtotChangeW4"] = df.apply(computeMemIndexChange,waveNumber=4,axis=1)
	df["memtotChangeW5"] = df.apply(computeMemIndexChange,waveNumber=5,axis=1)
	df["memtotb_n_4"] = df["memtotb_n_3"]
	df["memtotb_n_5"] = df["memtotb_n_4"]
	return df


def getTreatmentGroups(df, indVariable, waveNumber):
	varCurrWave = "{}_b_{}".format(indVariable, waveNumber)
	varPrevWave = "{}_b_{}".format(indVariable, waveNumber-1)

	treatmentIndexes = df.index[df[varCurrWave] == 1].tolist()
	controlIndexes = df.index[df[varCurrWave] == 0].tolist()	
	
	for i in treatmentIndexes:
		if (df.loc[i][varPrevWave]==1) or (df.loc[i][varPrevWave]==np.nan):
			treatmentIndexes.remove(i)

	for i in controlIndexes:
		if df.loc[i][varPrevWave]==1 or (df.loc[i][varPrevWave]==np.nan):
			controlIndexes.remove(i)
	print "G"
	print len(controlIndexes)
	print len(treatmentIndexes)
	return [controlIndexes, treatmentIndexes]


def ComputeCostMatrix(df, treatmentGroups, indVariable, waveNumber):
	controlIndexes = treatmentGroups[0]
	treatmentIndexes = treatmentGroups[1]

	cols = df.columns.tolist()
	cols.remove('idauniq')
	pattern = r"[a-zA-Z0-9]*_n_{}$".format(waveNumber)
	confounders = []
	for colName in cols:
		if (re.match(pattern, colName) and not (indVariable in colName)):
			confounders.append(colName)

	confDF = df[confounders]

	numTreat = len(treatmentIndexes)
	numControl = len(controlIndexes)
	C = np.zeros(shape = (numTreat, numControl))
	for i in range(numTreat):
		for j in range(numControl):
			C[i,j] = computeDistance(confDF.loc[treatmentIndexes[i]].values, confDF.loc[controlIndexes[j]].values)

	return C



def run_cmd(cmd, working_directory=None):
	if working_directory!= None:
		try:
			output = subprocess.check_output(cmd,shell=True,cwd=working_directory)
			print "output:"+output
		except:
			print "failed:"+cmd
			# pass
	else:
		try:
			output = subprocess.check_output(cmd,shell=True)
			print(output)
		except:
			print "failed:"+cmd
			# pass


def performMatching(C):

	r,c = C.shape
	with open('matrix.txt', 'w') as f:
		f.write("{} {}\n".format(r,c))
		for i in range(0,r):
			for j in range(0,c):
				f.write( "{} ".format(C[i][j]))
			f.write("\n")

	command = "hungarian/test"
	run_cmd(command)
	with open('matching.txt', 'r') as f:
		indexes = []
		for line in f:
			words = line.rstrip('\n').split(',')
			L = int(words[0])
			R = int(words[1])
			if R!= -1:
				pair = (L,R)
				indexes.append(pair)

	# m = Munkres()
	# indexes = m.compute(C)
	return indexes

def getTargetValues(df, treatmentGroups, indexes, waveNumber):
	memTotChangeVar = "memtotChangeW{}".format(waveNumber)
	controlIndexes = treatmentGroups[0]
	treatmentIndexes = treatmentGroups[1]
	memtotT = [  df.loc[treatmentIndexes[i[0]]][memTotChangeVar]  for i in indexes]
	memtotC = [  df.loc[controlIndexes[i[1]]][memTotChangeVar]  for i in indexes]
	return [memtotC, memtotT]

def computePValue(X,Y):
	res= scipy.stats.wilcoxon(X,Y,"wilcox")
	pVal = res[1]
	return pVal


def f():
	start_time = time.time()

	df = readData()
	df = preProcessData(df)
	indVariables = ["heacta", "heactb", "heactc", "scorg03", "scorg05", "scorg06","scorg07",
					"scako","heskb"]

	# indVariables = ["heacta"]

	# indVariables = ["scako"]

	# indVariables=indVariables[0:3]
	# targetValues = []
	pVals = {}
	for indVariable in indVariables:
		pVals[indVariable] = []
	
	for indVariable in indVariables:
		s =time.time()
		print indVariable
		for waveNumber in [4,5]:
			print waveNumber
			treatmentGroups = getTreatmentGroups(df,indVariable, waveNumber)
			C= ComputeCostMatrix(df, treatmentGroups, indVariable, waveNumber)
			matchedPairs = performMatching(C)
			targetValues = getTargetValues(df,treatmentGroups, matchedPairs, waveNumber)

			pval = computePValue(targetValues[0], targetValues[1])
			print pval
			pVals[indVariable].append(pval)	
		elapsedTime = time.time()-s
		print "processing time:", elapsedTime/60		

	return pVals

if __name__ == "__main__":
	print "a"