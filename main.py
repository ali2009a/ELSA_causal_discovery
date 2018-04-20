import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
from munkres import Munkres
import scipy
import re

basePath = "/home/ali/Downloads/UKDA-5050-stata (2)/stata/stata13_se"
REFUSAL=-9
DONT_KNOW=-8
NOT_APPLICABLE=-1
SCHD_NOT_APPLICABLE=-2

NOT_ASKED=-3

NOT_IMPUTED = -999.0
NON_SAMPLE = -998.0
INST_RESPONDENT=  -995.0
def preProcessPhysicalActivity(value):
	if ( (value == REFUSAL) or (value == DONT_KNOW) or 
		(value ==  NOT_APPLICABLE)):
		return np.nan
	elif ( (value == 3) or (value == 4)):
		return 0
	elif ( (value == 1) or (value == 2)):
		return 1


def preProcessBinaryVariables(value):
	if ( (value == REFUSAL) or (value == DONT_KNOW) or 
		(value ==  NOT_APPLICABLE) or (value== SCHD_NOT_APPLICABLE)):
		return np.nan
	return value

def preProcessCovariates(value):
	if ( (value == REFUSAL) or (value == DONT_KNOW) or 
		(value ==  NOT_APPLICABLE) or (value== SCHD_NOT_APPLICABLE)):
		return np.nan
	return value

def preProcessWealthDecile(value):
	if ( (value == NOT_IMPUTED) or (value == NON_SAMPLE) or 
		(value ==  INST_RESPONDENT)):
		return np.nan
	return value


def preProcessToboccoUse(value):
	if ( (value == REFUSAL) or (value == DONT_KNOW) or (value== SCHD_NOT_APPLICABLE)):
		return np.nan
	elif (value == NOT_APPLICABLE):
		return 0
	else:
		return value

def preProcessMemIndex(value):
	if ( (value == REFUSAL) or (value == NOT_ASKED) or (value== SCHD_NOT_APPLICABLE) or 
		(value == NOT_APPLICABLE) or (value == -7)):
		return np.nan
	else:
		return value

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



	df["heacta"]=df["heacta"].apply(preProcessPhysicalActivity)
	df["heactb"]=df["heactb"].apply(preProcessPhysicalActivity)
	df["heactc"]=df["heactc"].apply(preProcessPhysicalActivity)

	df["scorg03"] = df["scorg03"].apply(preProcessBinaryVariables)
	df["scorg05"] = df["scorg05"].apply(preProcessBinaryVariables)
	df["scorg06"] = df["scorg06"].apply(preProcessBinaryVariables)
	df["scorg07"] = df["scorg07"].apply(preProcessBinaryVariables)

	df["hegenh"] = df["hegenh"].apply(preProcessCovariates)
	df["scfrda"] = df["scfrda"].apply(preProcessCovariates)
	df["scfrdg"] = df["scfrdg"].apply(preProcessCovariates)
	df["scako"] = df["scako"].apply(preProcessCovariates)
	df["heskb"] = df["heskb"].apply(preProcessToboccoUse)
	df["indager"] = df["indager"].apply(preProcessCovariates)
	df["dhsex"] = df["dhsex"].apply(preProcessCovariates)
	df["scfrdm"] = df["scfrdm"].apply(preProcessCovariates)
	

	df["totwq10_bu_s"] = df["totwq10_bu_s"].apply(preProcessWealthDecile)

	df["memtotb"] = df["memtotb"].apply(preProcessMemIndex)


	df = df.rename(columns = {'hegenh':'hehelf'})
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



	df["heacta"]=df["heacta"].apply(preProcessPhysicalActivity)
	df["heactb"]=df["heactb"].apply(preProcessPhysicalActivity)
	df["heactc"]=df["heactc"].apply(preProcessPhysicalActivity)

	df["scorg03"] = df["scorg03"].apply(preProcessBinaryVariables)
	df["scorg05"] = df["scorg05"].apply(preProcessBinaryVariables)
	df["scorg06"] = df["scorg06"].apply(preProcessBinaryVariables)
	df["scorg07"] = df["scorg07"].apply(preProcessBinaryVariables)

	df["hehelf"] = df["hehelf"].apply(preProcessCovariates)
	df["scfrda"] = df["scfrda"].apply(preProcessCovariates)
	df["scfrdg"] = df["scfrdg"].apply(preProcessCovariates)
	df["scako"] = df["scako"].apply(preProcessCovariates)
	df["heskb"] = df["heskb"].apply(preProcessToboccoUse)
	df["indager"] = df["indager"].apply(preProcessCovariates)
	df["dhsex"] = df["dhsex"].apply(preProcessCovariates)
	df["scfrdm"] = df["scfrdm"].apply(preProcessCovariates)
	

	df["totwq10_bu_s"] = df["totwq10_bu_s"].apply(preProcessWealthDecile)

	df["memtotb"] = df["memtotb"].apply(preProcessMemIndex)

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



	df["heacta"]=df["heacta"].apply(preProcessPhysicalActivity)
	df["heactb"]=df["heactb"].apply(preProcessPhysicalActivity)
	df["heactc"]=df["heactc"].apply(preProcessPhysicalActivity)

	df["scorg03"] = df["scorg03"].apply(preProcessBinaryVariables)
	df["scorg05"] = df["scorg05"].apply(preProcessBinaryVariables)
	df["scorg06"] = df["scorg06"].apply(preProcessBinaryVariables)
	df["scorg07"] = df["scorg07"].apply(preProcessBinaryVariables)

	df["hehelf"] = df["hehelf"].apply(preProcessCovariates)
	df["scfrda"] = df["scfrda"].apply(preProcessCovariates)
	df["scfrdg"] = df["scfrdg"].apply(preProcessCovariates)
	df["scako"] = df["scako"].apply(preProcessCovariates)
	df["heskb"] = df["heskb"].apply(preProcessToboccoUse)
	df["indager"] = df["indager"].apply(preProcessCovariates)
	df["dhsex"] = df["dhsex"].apply(preProcessCovariates)
	df["scfrdm"] = df["scfrdm"].apply(preProcessCovariates)
	

	df["totwq10_bu_s"] = df["totwq10_bu_s"].apply(preProcessWealthDecile)

	df["memtotb"] = df["memtotb"].apply(preProcessMemIndex)

	# df = df.ix[0:50,:]
	return df


def readData():
	df3 = readWave3Data(basePath)
	df4 = readWave4Data(basePath)
	df5 = readWave5Data(basePath)

	df34 = pd.merge(df3, df4, how='inner', on=['idauniq'],suffixes=('_3', ''))
	df345 = pd.merge(df34, df5, how='inner', on=['idauniq'],suffixes=('_4', '_5'))

	return df345


def computeMemIndexChange(row):
	return row["memtotb_4"] - row["memtotb_3"]


def normalizeData(df):
	cols = list(df.columns)
	cols.remove('idauniq')
	for col in cols:
	    col_norm = col + '_n'
	    df[col_norm] = (df[col] - df[col].min())/(df[col].max()- df[col].min())
	return df

def computeDistance(row1,row2,waveNumber, indVariable):
	cols = row1.keys().tolist()
	cols.remove('idauniq')

	pattern = r"[a-zA-Z0-9]*_{}_n$".format(waveNumber)
	confounders = []
	for colName in cols:
		if (re.match(pattern, colName) and not (indVariable in colName)):
			confounders.append(colName)

	cost = 0;
	for colName in confounders:
		if (not np.isnan(row1[colName])) and (not np.isnan(row2[colName])):
			# print row1[colName], row2[colName]
			cost+= abs(row1[colName] - row2[colName])
	return cost

def ComputeCostMatrix():
	df = readData()
	df = normalizeData(df)
	df= df.dropna(axis=0, how="any")
	
	df["memtotChangeW4"] = df.apply(computeMemIndexChange,axis=1)
	treatmentIndexes = df.index[df["heacta_3"] == 1].tolist()
	controlIndexes = df.index[df["heacta_3"] == 0].tolist()	

	treatmentIndexes = treatmentIndexes[0:100]
	controlIndexes = controlIndexes[0:200]

	numTreat = len(treatmentIndexes)
	numControl = len(controlIndexes)
	C = np.zeros(shape = (numTreat, numControl))
	for i in range(numTreat):
		for j in range(numControl):
			C[i,j] = computeDistance(df.loc[treatmentIndexes[i]], df.loc[controlIndexes[j]],4,"heacta")

	m = Munkres()
	indexes = m.compute(C)
	memtotT = [  df.loc[treatmentIndexes[i[0]]]["memtotChangeW4"]  for i in indexes]
	memtotC = [  df.loc[controlIndexes[i[1]]]["memtotChangeW4"]  for i in indexes]


	res= scipy.stats.wilcoxon(T[0],T[1],zero_method-"wilcox")
	pVal = res[1]

	return pVal


if __name__ == "__main__":
	print "a"