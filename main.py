import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt

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
	print "value", value
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
		(value == NOT_APPLICABLE)):
		return np.nan
	else:
		return value

def readWaveData(basePath, waveNumber):
	w3Core = pd.read_stata("{}/wave_{}_elsa_data.dta".format(basePath, waveNumber),convert_categoricals=False)
	w3Drv =  pd.read_stata("{}/wave_{}_ifs_derived_variables.dta".format(basePath, waveNumber),convert_categoricals=False)
	w3FinDrv = pd.read_stata('{}/wave_{}_financial_derived_variables.dta'.format(basePath, waveNumber),convert_categoricals=False)
	
	df = w3Core[["heacta","heactb","heactc", "scorg03", "scorg06", "scorg05", "scorg07", "hegenh",
				 "scfrda" , "scfrdg","scako", "heskb", "indager", "dhsex" , "scfrdm" ]]
	df["memtotb"] = w3Drv.ix[:,"memtotb"]
	df["totwq10_bu_s"] = w3FinDrv.ix[:, "totwq10_bu_s"]



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

	df = df.ix[0:50,:]
	return df

def readData():
	dfw3 = readWaveData(basePath , 3)
	dfw4 = readWaveData(basePath , 4)
	dfw5 = readWaveData(basePath , 5)


if __name__ == "__main__":
	
	print a
	# print w3Core["heskb"]
	# print w3Core["indager"].cat.codes

	# print w3Core['scorg06']
	# print w3Drv["memtotb"]

	# print w3FinDrv["totwq10_bu_s"]
	# print w3FinDrv["nfwq10_bu_s"]