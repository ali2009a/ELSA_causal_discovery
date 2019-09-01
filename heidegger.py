import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
# from munkres import Munkres
import scipy.stats
import re
import time
import subprocess
import os
from features import *
import scipy.signal as ss
import numpy as np
import math
from numpy import diff
from sklearn.cluster import KMeans
from tqdm import tqdm
dfPath = "df.pkl"
nanLabelPath= "nanLabel.pkl"
import pickle
from sklearn.cluster import AgglomerativeClustering

#"scako" was removed because wave 1 had different scale
trtmntVar = set(["scfrda","scfrdg","scfrdm","heacta", "heactb","heactc", "scorg03","scorg06","scorg05","scorg07","heskb"]) #11
confoundersVar = set(["indager", "hehelf","dhsex","totwq10_bu_s"])   #4
binaryVariables = set(["scorg03","scorg06","scorg05","scorg07","dhsex"])
targetVar = set(["memIndex"])  #1
auxVar = set(["cfdscr","cflisen", "cflisd","cfdatd"])
drvVar = set(["memIndexChange", "baseMemIndex"])  #1
allVar = trtmntVar|confoundersVar|targetVar
WAVE_NUM=7
weights = {"scfrda":1,"scfrdg":1,"scfrdm":1,"heacta":1, "heactb":1,"heactc":1, "scorg03":1,"scorg06":1,"scorg05":1,"scorg07":1,"heskb":1,"indager":2, "hehelf":1,"dhsex":1,"totwq10_bu_s":1,"baseMemIndex":1}
reverseVars = ["heacta", "heactb", "heactc", "scfrdg", "scfrda", "hehelf"]
# continousVars = ["heskb, scfrdm"]

basePath = "/home/ali/Downloads/UKDA-5050-stata (2)/stata/stata13_se"
REFUSAL=-9
DONT_KNOW=-8
NOT_APPLICABLE=-1
SCHD_NOT_APPLICABLE=-2

NOT_ASKED=-3

NOT_IMPUTED = -999.0
NON_SAMPLE = -998.0
INST_RESPONDENT=  -995.0


def returnSample(df):
    df = df.sample(20)
    varNames = ["heactb","scfrdm","scorg06"]
    col_list = []
    
        # for num in range(1,8):
            

    for num in range(1,8):
        
        for var in varNames:
            col_list.append("{}_b_{}".format(var,num))
        col_list.append( "memIndex_{}".format(num))     
    return df[col_list]

def removeMedianValues(df):
    for i in range(1,8):
        columnName = "scfrdm_b_{}".format(i)
        df= df.drop( df[ (df[columnName]==2)].index)
    return df

def report(df, var):
    for i in [3,4,5]:
        print "wave",i
        print "min", df["{}_{}".format(var, i)].min()
        print "max", df["{}_{}".format(var, i)].max()
        print "mean", df["{}_{}".format(var, i)].mean()
        print "std", df["{}_{}".format(var, i)].std()

def harmonizeData(df):
    # print allVar
    for var in (trtmntVar|confoundersVar):
        df[var] = df[var].apply((globals()[var].harmonize))
    
    # for var in ["heacta", "heactb", "heactc", "scako", "heskb"]:
    #     df[var] = df[var].apply((globals()[var].binarize))

    return df



def binarizeData(df):
    # pattern = r"[a-zA-Z0-9_]*_n$"
    # for var in trtmntVar:
    #   if not re.match(pattern, var):  
    #       col_bin = var + '_b'
    #       df[col_bin] = df[var].apply((globals()[var].binarize))
    return df


def normalizeData(df):
    # cols = list(df.columns)
    # cols.remove('idauniq')
    # for col in cols:
    #     col_norm = col + '_n'
    #     df[col_norm] = (df[col] - df[col].min())/(df[col].max()- df[col].min())

    for var in (trtmntVar|confoundersVar|targetVar):
        dfs=[]
        for i in range(1,8):
            col = "{}_{}".format(var,i)
            dfs.append(pd.DataFrame( {var: df[col]}))
        mergedDf = pd.concat(dfs)
        nonZeros = mergedDf[mergedDf[var]!=0][var]
        NZMean = nonZeros.mean()
        mean= mergedDf[var].mean()
        std = mergedDf[var].std()
        minValue = mergedDf[var].min()
        maxValue = mergedDf[var].max()
        mergedDf.to_csv("rawValues_{}".format(var), index=False)
        for i in range(1,8):
            col = "{}_{}".format(var,i)
            col_norm = "{}_n_{}".format(var,i)
            if(var == "heskb"):
                df[col_norm] = 1/(1+ np.exp(-(df[col]-(NZMean/2))*(0.5)))
            elif var == "scfrdm":
                df[col_norm] = 1/(1+ np.exp(-(df[col]-(NZMean/2))*(3))) 
            else:
                df[col_norm] = (df[col] - minValue)/(maxValue - minValue)
        if not var in targetVar:
            for i in range(1,8):
                # print "removed"
                col = "{}_{}".format(var,i)
                # print col
                df = df.drop(columns=col)
    return df



def readWave1Data(basePath):
    waveNumber=1
    Core = pd.read_stata("{}/wave_1_core_data_v3.dta".format(basePath, waveNumber),convert_categoricals=False)
    Drv =  pd.read_stata("{}/wave_1_ifs_derived_variables.dta".format(basePath, waveNumber),convert_categoricals=False)
    FinDrv = pd.read_stata('{}/wave_1_financial_derived_variables.dta'.format(basePath, waveNumber),convert_categoricals=False)
    
    s1 = pd.merge(Core, Drv, how='inner', on=['idauniq'])
    df = pd.merge(s1, FinDrv, how='inner', on=['idauniq'])
    
    df = df.rename(columns = {'scorg3':'scorg03'})
    df = df.rename(columns = {'scorg5':'scorg05'})
    df = df.rename(columns = {'scorg6':'scorg06'})
    df = df.rename(columns = {'scorg7':'scorg07'})
    
    col_list = ["idauniq"] + list(trtmntVar) + list(confoundersVar) + list(auxVar); 
    df = df [col_list] 

    df = addMemIndex(df)
    df = removeAuxVars(df)

    df = harmonizeData(df)
    # df = normalizeData(df)
    df = binarizeData(df)
    # df= addSuffix(df,1)
    # df = df.ix[0:50,:]
    return df

def addSuffix(df, num):
    for var in (trtmntVar|confoundersVar|targetVar):
        newName = "{}_{}".format(var,num)
        df = df.rename(columns = {var:newName})
    return df


def readWave2Data(basePath):
    waveNumber=2
    Core = pd.read_stata("{}/wave_2_core_data_v4.dta".format(basePath, waveNumber),convert_categoricals=False)
    Drv =  pd.read_stata("{}/wave_2_derived_variables.dta".format(basePath, waveNumber),convert_categoricals=False)
    FinDrv = pd.read_stata('{}/wave_{}_financial_derived_variables.dta'.format(basePath, waveNumber),convert_categoricals=False)
    

    s1 = pd.merge(Core, Drv, how='inner', on=['idauniq'])
    df = pd.merge(s1, FinDrv, how='inner', on=['idauniq'])

    df = df.rename(columns = {'HeActa':'heacta'})
    df = df.rename(columns = {'HeActb':'heactb'})
    df = df.rename(columns = {'HeActc':'heactc'})
    df = df.rename(columns = {'Hehelf':'hehelf'})
    df = df.rename(columns = {'HeSkb':'heskb'})
    df = df.rename(columns = {'DhSex':'dhsex'})

    df = df.rename(columns = {'CfDScr':'cfdscr'})
    df = df.rename(columns = {'CfLisEn':'cflisen'})
    df = df.rename(columns = {'CfLisD':'cflisd'})
    df = df.rename(columns = {'CfDatD':'cfdatd'})

    col_list = ["idauniq"] + list(trtmntVar) + list(confoundersVar) + list(auxVar); 
    df = df [col_list] 

    df = addMemIndex(df)
    df = removeAuxVars(df)

    df = harmonizeData(df)
    # df = normalizeData(df)
    df = binarizeData(df)
    # df = df.ix[0:50,:]
    return df


def readWave3Data(basePath):
    waveNumber=3
    Core = pd.read_stata("{}/wave_{}_elsa_data.dta".format(basePath, waveNumber),convert_categoricals=False)
    Drv =  pd.read_stata("{}/wave_{}_ifs_derived_variables.dta".format(basePath, waveNumber),convert_categoricals=False)
    FinDrv = pd.read_stata('{}/wave_{}_financial_derived_variables.dta'.format(basePath, waveNumber),convert_categoricals=False)
    

    s1 = pd.merge(Core, Drv, how='inner', on=['idauniq'])
    df = pd.merge(s1, FinDrv, how='inner', on=['idauniq'])

    df = df.rename(columns = {'hegenh':'hehelf'})
    col_list = ["idauniq"] + list(trtmntVar) + list(confoundersVar) + list(auxVar); 
    df = df [col_list] 

    df = addMemIndex(df)
    df = removeAuxVars(df)

    df = harmonizeData(df)
    # df = normalizeData(df)
    df = binarizeData(df)
    # df = df.ix[0:50,:]
    return df



def readWave4Data(basePath):
    waveNumber=4
    Core = pd.read_stata("{}/wave_{}_elsa_data.dta".format(basePath, waveNumber),convert_categoricals=False)
    Drv =  pd.read_stata("{}/wave_{}_ifs_derived_variables.dta".format(basePath, waveNumber),convert_categoricals=False)
    FinDrv = pd.read_stata('{}/wave_{}_financial_derived_variables.dta'.format(basePath, waveNumber),convert_categoricals=False)
    
    s1 = pd.merge(Core, Drv, how='inner', on=['idauniq'])
    df = pd.merge(s1, FinDrv, how='inner', on=['idauniq'])

    col_list = ["idauniq"] + list(trtmntVar) + list(confoundersVar) + list(auxVar); 
    df = df [col_list] 

    df = addMemIndex(df)
    df = removeAuxVars(df)

    df = harmonizeData(df)
    # df = normalizeData(df)
    df = binarizeData(df)
    # df = df.ix[0:50,:]
    return df


def readWave5Data(basePath):
    waveNumber=5
    Core = pd.read_stata("{}/wave_{}_elsa_data.dta".format(basePath, waveNumber),convert_categoricals=False)
    Drv =  pd.read_stata("{}/wave_{}_ifs_derived_variables.dta".format(basePath, waveNumber),convert_categoricals=False)
    FinDrv = pd.read_stata('{}/wave_{}_financial_derived_variables.dta'.format(basePath, waveNumber),convert_categoricals=False)

    s1 = pd.merge(Core, Drv, how='inner', on=['idauniq'])
    df = pd.merge(s1, FinDrv, how='inner', on=['idauniq'])

    col_list = ["idauniq"] + list(trtmntVar) + list(confoundersVar) + list(auxVar); 
    df = df [col_list] 

    df = addMemIndex(df)
    df = removeAuxVars(df)

    df = harmonizeData(df)
    # df = normalizeData(df)
    df = binarizeData(df)
    # df = df.ix[0:50,:]
    return df


def readWave6Data(basePath):
    waveNumber=6
    w6Core = pd.read_stata("{}/wave_{}_elsa_data_v2.dta".format(basePath, waveNumber),convert_categoricals=False)
    w6Drv =  pd.read_stata("{}/wave_{}_ifs_derived_variables.dta".format(basePath, waveNumber),convert_categoricals=False)
    w6FinDrv = pd.read_stata('{}/wave_{}_financial_derived_variables.dta'.format(basePath, waveNumber),convert_categoricals=False)
    

    s1 = pd.merge(w6Core, w6Drv, how='inner', on=['idauniq'])
    df = pd.merge(s1, w6FinDrv, how='inner', on=['idauniq'])

    # df = df.rename(columns = {'hegenh':'hehelf'})
    df = df.rename(columns = {'HeActa':'heacta'})
    df = df.rename(columns = {'HeActb':'heactb'})
    df = df.rename(columns = {'HeActc':'heactc'})
    df = df.rename(columns = {'Hehelf':'hehelf'})
    df = df.rename(columns = {'HeSkb':'heskb'})
    df = df.rename(columns = {'DhSex':'dhsex'})

    df = df.rename(columns = {'CfDScr':'cfdscr'})
    df = df.rename(columns = {'CfLisEn':'cflisen'})
    df = df.rename(columns = {'CfLisD':'cflisd'})
    df = df.rename(columns = {'CfDatD':'cfdatd'})



    # col_list = ["idauniq","heacta","heactb","heactc", "scorg03", "scorg06", "scorg05", "scorg07", "hehelf",
    #            "scfrda" , "scfrdg","scako", "heskb", "indager", "dhsex" , "scfrdm", "memtotb","totwq10_bu_s",  ]
    
    col_list = ["idauniq"] + list(trtmntVar) + list(confoundersVar) + list(auxVar); 
    df = df [col_list] 

    df = addMemIndex(df)
    df = removeAuxVars(df)

    df = harmonizeData(df)
    # df = normalizeData(df)
    df = binarizeData(df)
    # df = df.ix[0:50,:]
    return df




def readWave7Data(basePath):
    waveNumber=7
    w6Core = pd.read_stata("{}/wave_{}_elsa_data.dta".format(basePath, waveNumber),convert_categoricals=False)
    w6Drv =  pd.read_stata("{}/wave_{}_ifs_derived_variables.dta".format(basePath, waveNumber),convert_categoricals=False)
    w6FinDrv = pd.read_stata('{}/wave_7_financial_derived_variables.dta'.format(basePath, waveNumber),convert_categoricals=False)
    

    s1 = pd.merge(w6Core, w6Drv, how='inner', on=['idauniq'])
    df = pd.merge(s1, w6FinDrv, how='inner', on=['idauniq'])

    df = df.rename(columns = {'HeActa':'heacta'})
    df = df.rename(columns = {'HeActb':'heactb'})
    df = df.rename(columns = {'HeActc':'heactc'})
    df = df.rename(columns = {'Hehelf':'hehelf'})
    df = df.rename(columns = {'HeSkb':'heskb'})
    df = df.rename(columns = {'DhSex':'dhsex'})
    df = df.rename(columns = {'scfrdl':'scfrdm'})



    df = df.rename(columns = {'CfDScr':'cfdscr'})
    df = df.rename(columns = {'CfLisEn':'cflisen'})
    df = df.rename(columns = {'CfLisD':'cflisd'})
    df = df.rename(columns = {'CfDatD':'cfdatd'})



    # col_list = ["idauniq","heacta","heactb","heactc", "scorg03", "scorg06", "scorg05", "scorg07", "hehelf",
    #            "scfrda" , "scfrdg","scako", "heskb", "indager", "dhsex" , "scfrdm", "memtotb","totwq10_bu_s",  ]
    
    col_list = ["idauniq"] + list(trtmntVar) + list(confoundersVar) + list(auxVar); 
    df = df [col_list] 

    df = addMemIndex(df)
    df = removeAuxVars(df)

    df = harmonizeData(df)
    # df = normalizeData(df)
    df = binarizeData(df)
    # df = df.ix[0:50,:]
    return df


def removeAuxVars(df):
    df= df.drop(list(auxVar),axis=1)
    return df


def readData(mergeMethod="outer"):
    df1 = readWave1Data(basePath)
    df2 = readWave2Data(basePath)
    df3 = readWave3Data(basePath)
    df4 = readWave4Data(basePath)
    df5 = readWave5Data(basePath)
    df6 = readWave6Data(basePath)
    df7 = readWave7Data(basePath)


    df12 = pd.merge(df1,  df2, how=mergeMethod, on=['idauniq'],suffixes=('_1', ''))
    df13 = pd.merge(df12, df3, how=mergeMethod, on=['idauniq'],suffixes=('_2', ''))
    df14 = pd.merge(df13, df4, how=mergeMethod, on=['idauniq'],suffixes=('_3', ''))
    df15 = pd.merge(df14, df5, how=mergeMethod, on=['idauniq'],suffixes=('_4', ''))
    df16 = pd.merge(df15, df6, how=mergeMethod, on=['idauniq'],suffixes=('_5', ''))
    df17 = pd.merge(df16, df7, how=mergeMethod, on=['idauniq'],suffixes=('_6', '_7'))

    # df34 = pd.merge(df3, df4, how='inner', on=['idauniq'],suffixes=('_3', ''))
    # df35 = pd.merge(df34, df5, how='inner', on=['idauniq'],suffixes=('_4', '_5'))

    return df17

    # return [df1, df2, df3, df4, df5, df6 , df7]

def addMemIndex(df):
    df["memIndex"] = df.apply(computeMemIndex, axis=1)
    df= df.dropna(subset=["memIndex"])
    return df


def computeMemIndexChange(row, waveNumber):
    memtotVarCur = "memIndex_{}".format(waveNumber) 
    memtotVarPrev = "memIndex_{}".format(waveNumber-1)
    return row[memtotVarCur] - row[memtotVarPrev]



def computeMemIndex(row):
    if row["cfdatd"] == REFUSAL:
        return np.nan
    if row ["cflisd"] == DONT_KNOW:
        row ["cflisd"] = 0

    if row ["cflisen"] == DONT_KNOW:
        row ["cflisen"] = 0

    if (row ["cfdscr"]<0) or (row ["cflisd"]<0) or (row ["cflisen"]<0):
        return np.nan
    else:
        return row["cfdscr"] + row["cflisd"] + row["cflisen"]


def computeDistance(row1,row2, weights_local):
    diff  = row1 - row2
    diff = diff*weights_local
    diff = diff[~np.isnan(diff)]
    return np.linalg.norm(diff)


def preProcessData(df):
    # df= df.dropna(axis=0, how="any")
    # df= df.dropna(subset=["memIndex_1", "memIndex_", "memIndex", "memIndex","memIndex","memIndex","memIndex"])
    
    # df["memIndexChange_1"]=  np.nan
    # df["memIndexChange_2"] = df.apply(computeMemIndexChange,waveNumber=4,axis=1)
    # df["memIndexChange_3"] = df.apply(computeMemIndexChange,waveNumber=5,axis=1)
    # df["memIndexChange_4"] = df.apply(computeMemIndexChange,waveNumber=5,axis=1)
    # df["memIndexChange_5"] = df.apply(computeMemIndexChange,waveNumber=5,axis=1)
    # df["memIndexChange_6"] = df.apply(computeMemIndexChange,waveNumber=5,axis=1)
    # df["memIndexChange_7"] = df.apply(computeMemIndexChange,waveNumber=5,axis=1)


    # df["baseMemIndex_7"] = df["memIndex_6"]
    # df["baseMemIndex_6"] = df["memIndex_5"]
    # df["baseMemIndex_5"] = df["memIndex_4"]
    # df["baseMemIndex_4"] = df["memIndex_3"]
    # df["baseMemIndex_3"] = df["memIndex_2"]
    # df["baseMemIndex_2"] = df["memIndex_1"]
    # df["baseMemIndex_1"] = np.nan


    df = normalizeData(df)
    
    
    return df


def getTreatmentGroups(df, indVariable, waveNumber):
    varCurrWave = "{}_b_{}".format(indVariable, waveNumber)
    varPrevWave = "{}_b_{}".format(indVariable, waveNumber-1)
    memIndexChangeVar = "memIndexChange_{}".format(waveNumber)

    treatmentIndexes = df.index[df[varCurrWave] == 1].tolist()
    controlIndexes = df.index[df[varCurrWave] == 0].tolist()    
    
    for i in treatmentIndexes:
        if (df.loc[i][varPrevWave]==1) or (df.loc[i][varPrevWave]==np.nan) or (df.loc[i][memIndexChangeVar]==np.nan):
            treatmentIndexes.remove(i)

    for i in controlIndexes:
        if df.loc[i][varPrevWave]==1 or (df.loc[i][varPrevWave]==np.nan) or (df.loc[i][memIndexChangeVar]==np.nan):
            controlIndexes.remove(i)
    print "Group size:"
    print len(controlIndexes)
    print len(treatmentIndexes)
    return [controlIndexes, treatmentIndexes]


def ComputeCostMatrix(df, treatmentGroups, indVariable, waveNumber):
    controlIndexes = treatmentGroups[0]
    treatmentIndexes = treatmentGroups[1]

    # cols = df.columns.tolist()
    # cols.remove('idauniq')
    # pattern = r"[a-zA-Z0-9]*_n_{}$".format(waveNumber)
    confounders = []
    # for colName in cols:
    #   if (re.match(pattern, colName) and not (indVariable in colName)):
    #       confounders.append(colName)

    weights_local = []

    for var in ((trtmntVar| confoundersVar | set(["baseMemIndex"]))- set([indVariable])):
        if var in binaryVariables:
            colName = "{}_{}".format(var,waveNumber)
        else:
            colName= "{}_n_{}".format(var,waveNumber)
        confounders.append(colName) 
        if var == "indager":
            weights_local.append(3)
        elif var == "baseMemIndex":
            weights_local.append(3)
        else:
            weights_local.append(1)


    confDF = df[confounders]
    # print confounders

    numTreat = len(treatmentIndexes)
    numControl = len(controlIndexes)
    C = np.zeros(shape = (numTreat, numControl))
    for i in range(numTreat):
        for j in range(numControl):
            # print confDF.loc[treatmentIndexes[i]]
            # print confDF.loc[treatmentIndexes[i]].values
            # print weights_local
            C[i,j] = computeDistance(confDF.loc[treatmentIndexes[i]].values, confDF.loc[controlIndexes[j]].values,weights_local)

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
        costs = []
        for line in f:
            words = line.rstrip('\n').split(',')
            L = int(words[0])
            R = int(words[1])
            if R!= -1:
                pair = (L,R)
                indexes.append(pair)
                costs.append(C[L,R])

    costs = np.array(costs)
    passedPairs = [pair for idx, pair in enumerate(indexes) if costs[idx]< 0.01 ]           
    return passedPairs


    # L=3
    # costs = np.array(costs)
    # print costs.min(), costs.max()
    # logC= np.log(C.flatten())
    # threshold = (np.median(logC) - L*MAD(logC))
    # passThreshold = np.log(costs) < threshold
    # print threshold
    # passedPairs = [pair for idx, pair in enumerate(indexes) if passThreshold[idx] ]
    # return passedPairs


def doRandomMatching(k, labels, trtNUM):
    indices = np.where(labels==k)[0]
    trt = indices[np.where(indices<trtNUM)[0]]
    ctrl = indices[np.where(indices>=trtNUM)[0]]
    if (len(trt)<len(ctrl)):
        ctrl = np.random.choice(ctrl, len(trt), replace=False)
    else:
        trt = np.random.choice(trt, len(ctrl), replace=False)

    return zip(trt,ctrl)
    

def performMatching_RBD(C, trtNUM):
    CLUSTER_NUM=22
    model = AgglomerativeClustering(affinity='precomputed', n_clusters=CLUSTER_NUM, linkage='complete').fit(C)
    labels = model.labels_
    print labels
    finalPairs = []
    for k in range(0, CLUSTER_NUM):
        pairs= doRandomMatching(k, labels, trtNUM)
        finalPairs = finalPairs+ pairs
    
    print len(finalPairs) 


    costs = []
    for pair in finalPairs:
        costs.append(C[pair[0], pair[1]])
    print costs
    costs = np.array(costs)
    passedPairs = [pair for idx, pair in enumerate(finalPairs) if costs[idx]< 0.01 ]           
    print len(passedPairs)
    return passedPairs



def getTargetValues(df, treatmentGroups, indexes, waveNumber):
    memTotChangeVar = "memIndexChange_{}".format(waveNumber)
    controlIndexes = treatmentGroups[0]
    treatmentIndexes = treatmentGroups[1]
    memtotT = [  df.loc[treatmentIndexes[i[0]]][memTotChangeVar]  for i in indexes]
    memtotC = [  df.loc[controlIndexes[i[1]]][memTotChangeVar]  for i in indexes]
    return [memtotC, memtotT]

def computePValue(X,Y):
    res= scipy.stats.wilcoxon(X,Y,"wilcox")
    pVal = res[1]
    return pVal


def getVariableData(df, var):
    Xvars=[]
    for i in range(1,8):
        Xvar = "{}_{}".format(var, i)
        Xvars.append(Xvar)
    # print Xvars

    newDF = df[Xvars]

    # minVal= float("-inf")
    # maxVal = float("inf")
    # for var in Xvars:
    #   # print var
    #   if newDF[var].min()<minVal:
    #       minVal = newDF[var].min()
    #   if newDF[var].max()>maxVal:
    #       maxVal = newDF[var].max()
    
    # for var in Xvars:
    #   newDF[var] =  (newDF[var] - newDF[var].min())/(newDF[var].max()- newDF[var].min())  
            
    return newDF



def exportVariables(df):

    for var in (trtmntVar|targetVar):
        D = getVariableData(df, var)
        D.to_csv("{}_data.csv".format(var))

    return



    X = getVariableData(df, X)

def getVariableDataBinary(df, var):
    Xvars=[]
    for i in range(1,8):
        Xvar = "{}_b_{}".format(var, i)
        Xvars.append(Xvar)
    # print Xvars

    newDF = df[Xvars]

    minVal= float("-inf")
    maxVal = float("inf")
    for var in Xvars:
        # print var
        if newDF[var].min()<minVal:
            minVal = newDF[var].min()
        if newDF[var].max()>maxVal:
            maxVal = newDF[var].max()
    
    for var in Xvars:
        newDF[var] =  (newDF[var] - newDF[var].min())/(newDF[var].max()- newDF[var].min())  
            
    return newDF

def detectLag(a,b):
    # print a
    # print b
    # # result = ss.correlate(Xvec, Yvec, method="direct")
    # print "a", len(a)
    # print "b", len(b)

    # result= ss.correlate(a - np.mean(a), b - np.mean(b), method='direct')/(np.std(a)*np.std(b)*len(a))
    result= ss.correlate(a, b, method='direct')
    return result

def findWaveVars(inputList):
    pattern = r"[a-zA-Z0-9_]*_1$"
    for var in inputList:
        if re.match(pattern, var):
            print var 


def computeLag(df, X,Y):
    X = getVariableData(df, X)
    Y = getVariableData(df, Y)

    X = X.interpolate()

    res = []

    # lags = {}
    # counter= {}
    # for i in range(-7,8):
        # lags[i]=0
        # counter[i] = 0


    for i in range(0,len(Y)):
        lag = detectLag(diff(X.loc[i]), diff(Y.loc[i]))
        res.append(lag)
        # print X.loc[i]
        # print Y.loc[i]
        # print lag
        # if  not math.isnan(lag[1]):
        #   # print type(lag[1])
        #   # print lag[1]
        #   lags[lag[0]]+= lag[1]
        #   counter[lag[0]]+= 1

    return (X,Y,res)

    # for i in range(-7,8):
    #   print "lag: {} , sum: {:.2f}".format(i, lags[i]) 
    #   if counter[i]:
    #       print "\t avg: {0:.2f}".format(lags[i]/counter[i])  

    # return lags, counter


def computeLagForAllVars(df):
    for var in trtmntVar:
        print "examining ", var
        res= computeLag(df,var, list(targetVar)[0])


def f():
    # start_time = time.time()

    df = readData()
    df = preProcessData(df)



    pVals = {}
    for indVariable in trtmntVar:
        pVals[indVariable] = []
    
    for indVariable in trtmntVar:
    # for indVariable in ["heactb"]:
        s =time.time()
        print indVariable
        for waveNumber in [2,3,4,5,6,7]:
        # for waveNumber in [5]:
            print waveNumber
            treatmentGroups = getTreatmentGroups(df,indVariable, waveNumber)
            C= ComputeCostMatrix(df, treatmentGroups, indVariable, waveNumber)
            matchedPairs = performMatching(C)
            targetValues = getTargetValues(df,treatmentGroups, matchedPairs, waveNumber)

            pval = computePValue(targetValues[0], targetValues[1])
            print "pval", pval
            pVals[indVariable].append(pval) 
        elapsedTime = time.time()-s
        print "processing time:", elapsedTime/60        

    return pVals



def expandData(df):
    newDF =  pd.DataFrame( {'idauniq':df["idauniq"]})
    for var in (trtmntVar|confoundersVar|targetVar):
        for i in range(1,8):
            col = "{}_n_{}".format(var,i)
            newDF[col]=np.nan

    for i in range(1,8):
        for var in targetVar:
            col = "{}_{}".format(var,i)
            newDF[col]= np.nan

    offset=7
    for var in (trtmntVar|confoundersVar|targetVar):
        for i in range(1,8):
            col = "{}_n_{}".format(var,i)
            newCol = "{}_n_{}".format(var,i+offset)
            df = df.rename(columns={col: newCol})

    for i in range(1,8):
        for var in targetVar:
            col = "{}_{}".format(var,i)
            newCol = "{}_{}".format(var,i+offset)       
            df = df.rename(columns={col: newCol})

    finalDF = pd.merge(newDF,  df, how="inner", on=['idauniq'])     
    return finalDF


def interpolate(df):
    for index, row in df.iterrows():
        for var in (trtmntVar|confoundersVar|targetVar):
            cols= []
            for i in range(1,15):
                col = "{}_n_{}".format(var,i)
                cols.append(col)
            seq = df.loc[index, cols]
            f2  = pd.DataFrame(np.array([seq]), columns=cols)
            f2  = f2.interpolate(method ='linear', axis=1, limit_direction="both" )
            df.loc[index, cols] = f2.loc[0, cols]
    return df


def getTreatmentSignal():
    signal = np.array( [0, 1, 1])
    weights =  np.array( [1, 1, 1]) 

    return (signal, weights)

def getControlSignal():
    signal = np.array([0, 0, 0])
    weights = np.array( [1, 1, 1])  

    return (signal, weights)


MIExcludingPoints = 1  #must be equal or larger than 1(to exclude the MI for current wave)
def getMatchingWeights():
    weights = np.array( [1, 1, 1, 1, 1, 1, 1])
    return weights

def computeDistance(seq1, seq2, seq1Label, seq2Label=None, weights=None):
    # print "seq", seq
    # print "label", seqLabel
    alpha = 0.3
#   if weights is None:
#       F_index = 0.75
#       L = len(seq1)
#       weights= np.zeros(L)
#       for i in range(0,L):
#           weights[i] = F_index ** ( (L-1)-i);
    weights = weights[-len(seq1):]
    penalty =0.5
    sumDiff=0
    for i in range(len(seq1)):
        if seq1Label[i]:
            diff =np.abs(seq1[i]- seq2[i])
            diff =  (1-2*alpha)*diff+alpha
        elif (not seq2Label is None) and seq2Label[i]:
            diff =np.abs(seq1[i]- seq2[i])
            diff =  (1-2*alpha)*diff+alpha
        else:
            diff =np.abs(seq1[i]- seq2[i])
        if np.isnan(diff):
            diff= 1
        sumDiff = sumDiff + diff*weights[i]
    avgDiff = sumDiff/np.sum(weights)
    if np.isnan(avgDiff):
        print "found Nan"
    return avgDiff

def measureSimilarity(var, signal, df, nanLabel):
    [samplesNum, columnsNum] = df.shape
    distanceValues = np.empty((samplesNum*WAVE_NUM,3))   #to be fixed
    # distanceValues = np.empty((7*1000,3))
    distanceValues[:] = np.nan


    windowLength = len(signal[0])
    counter= 0
    for index in tqdm(range(0, len(df))):   #to be fixed
    # for index in tqdm(range(0,1000)):
        # print "index", index  
        for w in range(8,15):
        # for w in [11]:

            # print "w", w
            cols= []
            for i in range(w-(windowLength-1),w+1):

                col = "{}_n_{}".format(var,i)
                cols.append(col)

            seq = np.array(df.ix[index, cols])
            seqLabel = np.array(nanLabel.ix[index, cols])
            diff = computeDistance(seq, signal[0], seqLabel, weights=signal[1])
            distanceValues[counter,:]=  [int(index), int(w), diff]
            counter = counter+1     
    return distanceValues


def measureSimilarity2(var, signal, df, nanLabel):
    [samplesNum, columnsNum] = df.shape
    distanceValues = np.empty((samplesNum*WAVE_NUM,3))
    distanceValues[:] = np.nan
    winLen = len(signal[0])
    signalSeq= signal[0]
    weights = signal[1]
    seqLabel = np.zeros((winLen,), dtype=int)
    
    S = np.zeros(shape = (samplesNum*WAVE_NUM, winLen))
    SL = np.zeros(shape = (samplesNum*WAVE_NUM, winLen))
    counter= 0
    for index in tqdm(range(0, len(df))):
        for w in range(8,15):
            seqs= extractSeq(df, nanLabel, var, index, w, False, length = winLen)
            S[counter,:]= seqs[0]
            SL[counter,:]= seqs[1]
            distanceValues[counter,:]=  [int(index), int(w), np.nan]
            counter = counter+1

    alpha =0.3
    S = S.reshape(S.shape[0], 1, S.shape[1]) 
    SL = SL.reshape(SL.shape[0], 1, SL.shape[1])

    diff = np.abs(S- signalSeq )
    isNan = np.logical_or(SL, seqLabel)
    penalizedDiff = (1-isNan)*diff + (isNan)*((1-2*alpha)*diff+alpha)
    penalizedDiff =  np.isnan(penalizedDiff)*np.ones(penalizedDiff.shape, dtype=int) + (1 - np.isnan(penalizedDiff))*penalizedDiff
    weightedDiff  = penalizedDiff * weights  # to get the right most weights 
    costSum =np.sum(diff,axis=2)
    varCost = costSum / np.sum(getMatchingWeights())
    distanceValues[:,2] = varCost.reshape(len(distanceValues))
    return distanceValues




def checkIsNan(x):
   return np.isnan(x)


def identifyTreatmentGroup(var, signal, df, nanLabel):
    D = measureSimilarity(var, signal, df, nanLabel)


    

def MAD(data, axis=None):
    res= np.median(np.absolute(data - np.median(data, axis)), axis)
    return res/0.67


def fixReverseDir(df):
    for var in (reverseVars):
        dfs=[]
        for i in range(1,8):
            col = "{}_{}".format(var,i)
            dfs.append(pd.DataFrame( {var: df[col]}))
        mergedDf = pd.concat(dfs)
        maxValue = mergedDf[var].max()
        for i in range(1,8):
            col = "{}_{}".format(var,i)
            df[col] = maxValue - df[col]
    return df

def preprocess(df):
    df = fixReverseDir(df)
    df= normalizeData(df)
    df= expandData(df)
    nanLabel = df.apply(checkIsNan) 
    df = interpolate(df)
    return df, nanLabel


def detectOutliers(distanceInfo, nanLabel, var, string, L=3):
    D = distanceInfo[:,2].copy()
    outliersIndex, L =  tuneL(D, nanLabel, distanceInfo, var, string)
    return outliersIndex


def removeUnknownMI(outliersIndex, distanceInfo, nanLabel):
    filteredIndex = []
    for i in outliersIndex:
        index = int(distanceInfo[i,0])
        w = int(distanceInfo[i,1])
        col= "memIndex_{}".format(w)
        if (not nanLabel.loc[index,col]):
            filteredIndex.append(i)
    return filteredIndex


def computeAvgDistance(df, nanLabel, outliersIndexT, outliersIndexC, distanceInfoT, distanceInfoC, varSet, isTargetVar, i, j):
        sumDiff = 0
        for var in (varSet):
            seqs_T = extractSeq(df, nanLabel, var, distanceInfoT[outliersIndexT[i],0], distanceInfoT[outliersIndexT[i],1], isTargetVar, length = len(getMatchingWeights()))
            seqs_C = extractSeq(df, nanLabel, var, distanceInfoC[outliersIndexC[j],0], distanceInfoC[outliersIndexC[j],1], isTargetVar, length = len(getMatchingWeights()))
            diff = computeDistance(seqs_T[0], seqs_C[0], seqs_T[1], seqs_C[1], weights= getMatchingWeights())
            sumDiff = sumDiff + diff
        return sumDiff/len(varSet)  


def computeDistanceMatrix(df, nanLabel, trtVariable, outliersIndexT, outliersIndexC, distanceInfoT, distanceInfoC):
    C = np.zeros(shape = (len(outliersIndexT), len(outliersIndexC)))
    
    for i in tqdm(range(0,len(outliersIndexT))):
        print i
        for j in range(0, len(outliersIndexC)):
            distances = []
            trtDist= computeAvgDistance(df, nanLabel, outliersIndexT, outliersIndexC, distanceInfoT, distanceInfoC, trtmntVar-set([trtVariable]), False, i, j)
            #distances.append(d)
            confDist= computeAvgDistance(df, nanLabel, outliersIndexT, outliersIndexC, distanceInfoT, distanceInfoC, confoundersVar, False, i, j)
            #distances.append(d)
            targetDist=computeAvgDistance(df, nanLabel, outliersIndexT, outliersIndexC, distanceInfoT, distanceInfoC, targetVar, True, i, j)
            #distances.append(d)
            C[i,j]= (trtDist + confDist + 50*targetDist)/52
    return C    




def computeAvgDistance2(df, nanLabel, outliersIndexT, outliersIndexC, distanceInfoT, distanceInfoC, varSet, isTargetVar, trtSeq):
        alpha = 0.3
        anchorPoint = (np.where(trtSeq==1)[0][0]-1)
        if (anchorPoint<0):
            anchorPoint=0
        anchorDist = len(trtSeq)- anchorPoint

        anchorPoint2 = (np.where(trtSeq!=2)[0][0])
        anchorDist2 = len(trtSeq)- anchorPoint2   
        
        if (isTargetVar):
            extractLen = anchorDist
            # winLen= 1
            winLen= extractLen-1
            # effectiveWeights = np.ones(winLen)
            effectiveWeights = np.geomspace(1, (0.5)**(winLen-1) , num=winLen)
        else:
            extractLen = anchorDist2 
            winLen = extractLen
            effectiveWeights = np.ones(winLen)
        
        costSum = np.zeros(shape = (len(outliersIndexT), len(outliersIndexC)))

        for var in (varSet):
            T = np.zeros(shape = (len(outliersIndexT), winLen))
            C = np.zeros(shape = (len(outliersIndexC), winLen))
            TL = np.zeros(shape = (len(outliersIndexT), winLen))
            CL = np.zeros(shape = (len(outliersIndexC), winLen))
            
            for i in range(0,len(outliersIndexT)):
                seqs =extractSeq(df, nanLabel, var, distanceInfoT[outliersIndexT[i],0], distanceInfoT[outliersIndexT[i],1], isTargetVar, length = extractLen)
                T[i,:] = seqs[0][:winLen] # to discard memIndex for the current weight (right most one)
                TL[i,:] = seqs[1][:winLen]

            for j in range(0,len(outliersIndexC)):
                seqs =extractSeq(df, nanLabel, var, distanceInfoC[outliersIndexC[j],0], distanceInfoC[outliersIndexC[j],1], isTargetVar, length = extractLen)
                C[j,:] = seqs[0][:winLen]
                CL[j,:] = seqs[1][:winLen]
            
            T = T.reshape(T.shape[0], 1, T.shape[1])
            TL = TL.reshape(TL.shape[0], 1, TL.shape[1])
            diff = np.abs(T-C)  
            isNan = np.logical_or(TL,CL).astype(int)
            penalizedDiff = (1-isNan)*diff + (isNan)*((1-2*alpha)*diff+alpha)           
            penalizedDiff =  np.isnan(penalizedDiff).astype(int)*np.ones(penalizedDiff.shape, dtype=int) + (1 - np.isnan(penalizedDiff).astype(int))*np.nan_to_num(penalizedDiff)
            weightedDiff  = penalizedDiff * effectiveWeights
            aggregatedCost =np.sum(weightedDiff,axis=2)
            varCost = aggregatedCost / np.sum(effectiveWeights)
            costSum = costSum + varCost

        avgCost= costSum/float(len(varSet))
        return avgCost
def computeDistanceMatrix2(df, nanLabel, trtVariable, outliersIndexT, outliersIndexC, distanceInfoT, distanceInfoC, trtSeq):
    trtDist= computeAvgDistance2(df, nanLabel, outliersIndexT, outliersIndexC, distanceInfoT, distanceInfoC, trtmntVar-set([trtVariable]), False, trtSeq)
    confDist= computeAvgDistance2(df, nanLabel, outliersIndexT, outliersIndexC, distanceInfoT, distanceInfoC, confoundersVar, False, trtSeq)
    targetDist=computeAvgDistance2(df, nanLabel, outliersIndexT, outliersIndexC, distanceInfoT, distanceInfoC, targetVar, True, trtSeq)
    C= (trtDist + confDist + 198*targetDist)/200
    np.savetxt('X_Dist.csv', trtDist, delimiter=',')
    np.savetxt('Z_Dist.csv', confDist, delimiter=',')
    np.savetxt('Y_Dist.csv', targetDist, delimiter=',')

    return C





def computeAvgDistance2_RBD(df, nanLabel, outliersIndexT, outliersIndexC, distanceInfoT, distanceInfoC, varSet, isTargetVar, trtSeq):
        alpha = 0.3
        anchorPoint = (np.where(trtSeq==1)[0][0]-1)
        if (anchorPoint<0):
            anchorPoint=0
        anchorDist = len(trtSeq)- anchorPoint

        anchorPoint2 = (np.where(trtSeq!=2)[0][0])
        anchorDist2 = len(trtSeq)- anchorPoint2   
        
        if (isTargetVar):
            extractLen = anchorDist
            winLen= 1
            effectiveWeights = np.ones(winLen)
        else:
            extractLen = anchorDist2 
            winLen = extractLen
            effectiveWeights = np.ones(winLen)
        
        totalNum = len(outliersIndexT)+len(outliersIndexC)
        costSum = np.zeros(shape = (totalNum,totalNum))
        print costSum.shape
        for var in (varSet):
            T = np.zeros(shape = (totalNum, winLen))
            C = np.zeros(shape = (totalNum, winLen))
            TL = np.zeros(shape = (totalNum, winLen))
            CL = np.zeros(shape = (totalNum, winLen))
            
            ids = []
            for i in range(0,len(outliersIndexT)):
                ids.append((distanceInfoT[outliersIndexT[i],0], distanceInfoT[outliersIndexT[i],1]))

            for j in range(0,len(outliersIndexC)):
                ids.append((distanceInfoC[outliersIndexC[j],0], distanceInfoC[outliersIndexC[j],1]))


            for i, pair in enumerate(ids):
                seqs =extractSeq(df, nanLabel, var, pair[0], pair[1], isTargetVar, length = extractLen)
                T[i,:] = seqs[0][:winLen] # to discard memIndex for the current weight (right most one)
                TL[i,:] = seqs[1][:winLen]
            C = T.copy()
            CL = TL.copy()
            
            T = T.reshape(T.shape[0], 1, T.shape[1])
            TL = TL.reshape(TL.shape[0], 1, TL.shape[1])
            diff = np.abs(T-C)  
            isNan = np.logical_or(TL,CL).astype(int)
            penalizedDiff = (1-isNan)*diff + (isNan)*((1-2*alpha)*diff+alpha)           
            penalizedDiff =  np.isnan(penalizedDiff).astype(int)*np.ones(penalizedDiff.shape, dtype=int) + (1 - np.isnan(penalizedDiff).astype(int))*np.nan_to_num(penalizedDiff)
            weightedDiff  = penalizedDiff * effectiveWeights
            aggregatedCost =np.sum(weightedDiff,axis=2)
            varCost = aggregatedCost / np.sum(effectiveWeights)
            print varCost.shape
            costSum = costSum + varCost

        avgCost= costSum/float(len(varSet))
        return avgCost

def computeDistanceMatrix2_RBD(df, nanLabel, trtVariable, outliersIndexT, outliersIndexC, distanceInfoT, distanceInfoC, trtSeq):

    trtDist= computeAvgDistance2_RBD(df, nanLabel, outliersIndexT, outliersIndexC, distanceInfoT, distanceInfoC, trtmntVar-set([trtVariable]), False, trtSeq)
    confDist= computeAvgDistance2_RBD(df, nanLabel, outliersIndexT, outliersIndexC, distanceInfoT, distanceInfoC, confoundersVar, False, trtSeq)
    targetDist=computeAvgDistance2_RBD(df, nanLabel, outliersIndexT, outliersIndexC, distanceInfoT, distanceInfoC, targetVar, True, trtSeq)
    C= (trtDist + confDist + 198*targetDist)/200
    np.savetxt('X_Dist.csv', trtDist, delimiter=',')
    np.savetxt('Z_Dist.csv', confDist, delimiter=',')
    np.savetxt('Y_Dist.csv', targetDist, delimiter=',')

    return C



def dump2DMatrix(C, var, outputPrefix): 
    r,c = C.shape
    with open('{}_{}.txt'.format(outputPrefix, var), 'w') as f:
        f.write("{} {}\n".format(r,c))
        for i in range(0,r):
            for j in range(0,c):
                f.write( "{} ".format(C[i][j]))
            f.write("\n")


def heidegger():
    if (os.path.isfile(dfPath) and os.path.isfile(nanLabelPath)):
        df = pd.read_pickle(dfPath)
        nanLabel = pd.read_pickle(nanLabelPath)
    else:
        df = readData()
        df, nanLabel = preprocess(df)

    f = open("result.txt","w")

    for var in trtmntVar:
    #for var in ["heactb"]:
        print "evaluting {}".format(var)
        distanceInfoT = measureSimilarity(var, getTreatmentSignal(), df, nanLabel)
        distanceInfoC = measureSimilarity(var, getControlSignal(), df, nanLabel)
        
        #distanceInfoT[:,2].tofile("Treatment_Distance_{}.csv".format(var),sep=',') ## for debug
        #distanceInfoC[:,2].tofile("Control_Distance_{}.csv".format(var),sep=',')   ## for debug
        
        outliersIndexT = detectOutliers(distanceInfoT, nanLabel, var, "Treatment")
        outliersIndexC = detectOutliers(distanceInfoC, nanLabel, var, "Control")
        C = computeDistanceMatrix2(df, nanLabel, var, outliersIndexT, outliersIndexC, distanceInfoT, distanceInfoC)
        
        C.tofile("FlattenCostMatrix_{}.csv".format(var),sep=',')  ## for debug
        dump2DMatrix(C, var, "CostMatrix")                          ## for debug

        matchedPairs = performMatching(C)
        targetValues = extractTargetValues(df, matchedPairs, outliersIndexT, outliersIndexC,distanceInfoT, distanceInfoC, var)
        pval = computePValue(targetValues[0], targetValues[1])
        f.write("{}:{}\n".format(var, pval))
        f.flush()
        os.fsync(f.fileno())
    f.close()

from numpy import genfromtxt
def runHyps():
        sigs = genfromtxt('LowH_Hyp.csv', delimiter=',')
        sigs = sigs[1:,:-1]
        print sigs
        if (os.path.isfile(dfPath) and os.path.isfile(nanLabelPath)):
            df = pd.read_pickle(dfPath)
            nanLabel = pd.read_pickle(nanLabelPath)
        else:
            df = readData()
            df, nanLabel = preprocess(df)

        f = open("result.txt","w")
        #for var in trtmntVar:
        var = "heactb"
        # ctrlSignal = ([0,0,0,0,0,0,0],[1,1,1,1,1,1,1])
     

        for i in range(0,len(sigs)):
            trtSeq = sigs[i,:]
            if (not isValidHyp(trtSeq)):
                continue
            weights = (~(trtSeq==2)).astype(int)
            trtSignal = (trtSeq,weights)                
            print "evaluating:"
            print trtSignal
            distanceInfoT = measureSimilarity(var, trtSignal, df, nanLabel)
            outliersIndexT = detectOutliers(distanceInfoT, nanLabel, var, "Treatment")

            anchorPoint = (np.where(trtSeq!=2)[0][0])
            anchorDist = len(trtSeq)- anchorPoint
            ctrlSeq = np.zeros(anchorDist)
            ctrlWeights = np.ones(anchorDist)
            ctrlSignal = (ctrlSeq, ctrlWeights)

            distanceInfoC = measureSimilarity(var, ctrlSignal, df, nanLabel)
            outliersIndexC = detectOutliers(distanceInfoC, nanLabel, var, "Control")

            if (len(outliersIndexC)==0 or  len(outliersIndexT)==0 ):
                f.write("{0} pattern: {1} , {2}\n".format(var, trtSeq.astype(int), "NA - outlier detection returned zero samples"))
                f.flush()
                os.fsync(f.fileno())
                continue

            C = computeDistanceMatrix2(df, nanLabel, var, outliersIndexT, outliersIndexC, distanceInfoT, distanceInfoC, trtSeq)

            C.tofile("FlattenCostMatrix_{}.csv".format(var),sep=',')  ## for debug
            dump2DMatrix(C, var, "CostMatrix")                                                      ## for debug

            matchedPairs = performMatching(C)
            if (len(matchedPairs)<4):
                f.write("{0} pattern: {1} , {2}\n".format(var, trtSeq.astype(int), "NA - matching returned less than four samples"))
                f.flush()
                os.fsync(f.fileno())
                continue                




            [isBiased, meanVals] = isDCBiased(df, matchedPairs, outliersIndexT, outliersIndexC,distanceInfoT, distanceInfoC, var, trtSeq)
            print "in runHyps"
            print meanVals
             

            targetValues = extractTargetValues(df, matchedPairs, outliersIndexT, outliersIndexC,distanceInfoT, distanceInfoC, var)
            pval = computePValue(targetValues[0], targetValues[1])


            # if(isBiased):
            #     # f.write("{0} pattern: {1} , {2}\n".format(var, trtSeq.astype(int), "Dont Care terms are biased"))
                
            #     f.write("{0} pattern: {1} , pval={2:} , ACE={4: .2f} , n={3:d}; dont care terms are biased, meanVals={4:} \n".format(var, trtSeq.astype(int), pval, len(matchedPairs), np.mean(targetValues[1])- np.mean(targetValues[0]), str(meanVals)))
            #     f.flush()
            #     os.fsync(f.fileno())
            #     continue    
            print meanVals
            meanValsStr = str(meanVals)
            print "meanValsStr:"
            print meanValsStr
            print "pval={0:.5f} , n={1:d}\n".format(pval, len(matchedPairs))
            f.write("{0} pattern: {1}, pval={2:}, ACE={4: .2f}, n={3:d}, DCT Mean={5}\n".format(var, trtSeq.astype(int), pval, len(matchedPairs), np.mean(targetValues[1])- np.mean(targetValues[0]),meanValsStr))
            f.flush()
            os.fsync(f.fileno())


        f.close()
    
def extractTargetValues(df, matchedPairs, outliersIndexT, outliersIndexC,distanceInfoT, distanceInfoC, var, anchorDist):
    print anchorDist  
    memtotT = []
    memtotC = []
    memtotT_prev = []
    memtotC_prev = []


    T_ids= []
    C_ids = []

    matchedPairs_fixed = []
    for pair in  matchedPairs:
        matchedPairs_fixed.append((pair[0], pair[1]-len(outliersIndexT)))

    for pair in matchedPairs_fixed:
        index, w  = distanceInfoT[outliersIndexT[pair[0]],0], distanceInfoT[outliersIndexT[pair[0]], 1]
        w=int(w)
        index=  int(index)
        T_ids.append((index,w))
        col= "memIndex_{}".format(w)
        col_prev = "memIndex_{}".format(w-anchorDist+1)
        memtotT.append( df.loc[index, col]- df.loc[index, col_prev]                                                                  )
        memtotT_prev.append(df.loc[index, col_prev])

    for pair in matchedPairs_fixed:
        index, w  = distanceInfoC[outliersIndexC[pair[1]],0],distanceInfoC[outliersIndexC[pair[1]], 1]
        w=int(w)
        index=  int(index)
        C_ids.append((index,w))
        col= "memIndex_{}".format(w)
        col_prev = "memIndex_{}".format(w-anchorDist+1)
        memtotC.append( df.loc[index, col]-df.loc[index, col_prev])
        memtotC_prev.append( df.loc[index, col_prev])



    dump2DMatrix(np.array(T_ids), var, "Treatment_Ids")
    dump2DMatrix(np.array(C_ids), var, "Control_Ids")
    np.array(memtotT).tofile("memIndexValues_Treatment_{}.csv".format(var), sep=',')
    np.array(memtotC).tofile("memIndexValues_Control_{}.csv".format(var), sep=',')
    
    np.array(memtotT_prev).tofile("memIndexValues_base_Treatment_{}.csv".format(var), sep=',')
    np.array(memtotC_prev).tofile("memIndexValues_base_Control_{}.csv".format(var), sep=',')
    print "memtotC:{}, memtotT:{}, memtotC_prev:{}, memtotT_prev:{}".format(np.mean(memtotC), np.mean(memtotT), np.mean(memtotC_prev), np.mean(memtotT_prev))
    return [memtotC, memtotT, memtotC_prev, memtotT_prev]   



def isValidHyp(h):
    return (h==1).any()
def produceHistograms(df, nanLabel):
    for var in trtmntVar:
        distanceInfoT = measureSimilarity(var, getTreatmentSignal(), df, nanLabel)
        distanceInfoC = measureSimilarity(var, getControlSignal(), df, nanLabel)
        DT= distanceInfoT[:,2]
        DC= distanceInfoC[:,2]

        plt.figure()
        plt.hist(DT, bins=20) 
        title = "{} treatment".format(var)
        plt.title(title)
        plt.savefig("{}.png".format(title))
        DT.tofile("{}.csv".format(title),sep=',')

        plt.figure()
        plt.hist(DC, bins=20) 
        title = "{} control".format(var)
        plt.title(title)
        plt.savefig("{}.png".format(title))
        DC.tofile("{}.csv".format(title),sep=',')



def getprevWaveMIs(indexes, distanceInfo):
    memtot = []
    for pair in matchedPairs:
        index, w  = distanceInfo[indexes[i],0] ,distanceInfoT[indexes[i],1]
        w=int(w)
        index=  int(index)
        col= "memIndex_{}".format(w)
        col_prev = "memIndex_{}".format(w-1)
        memtot.append( df.loc[index, col_prev])

    return memtot


def isDCBiased(df, matchedPairs, outliersIndexT, outliersIndexC,distanceInfoT, distanceInfoC, var, trtSeq):
    isBiased=False
    meanValList = []
    for i in trtSeq:
        if i==2:
            resPair = haveSameDCDistrubtion(df, matchedPairs, outliersIndexT, outliersIndexC,distanceInfoT, distanceInfoC, var, len(trtSeq)-(i+1))
            resValue = resPair[0]
            meanVal = resPair[1]
            print "meanVal:{}".format(meanVal)          
            isBiased = isBiased or not resValue
            meanValList.append(round(meanVal,2))
    print "meanValList"
    print meanValList
    return [isBiased, meanValList]


def haveSameDCDistrubtion(df, matchedPairs, outliersIndexT, outliersIndexC,distanceInfoT, distanceInfoC, var, offset):
    trtVal = []
    ctrlVal = []

    for pair in matchedPairs:
        index, w  = distanceInfoT[outliersIndexT[pair[0]],0] ,distanceInfoT[outliersIndexT[pair[0]],1]
        w=int(w-offset)
        index=  int(index)
        col= "{}_n_{}".format(var, w)
        trtVal.append( df.loc[index, col])

    for pair in matchedPairs:
        index, w  = distanceInfoC[outliersIndexC[pair[1]],0] ,distanceInfoC[outliersIndexC[pair[1]],1]
        w=int(w-offset)
        index=  int(index)
        col= "{}_n_{}".format(var, w)
        ctrlVal.append( df.loc[index, col])


    print "Treatment Mean:{}".format(np.mean(trtVal))
    # pval = computePValue(ctrlVal, trtVal)
    # print pval
    return [np.mean(trtVal) < 0.9 and np.mean(trtVal) >0.1, np.mean(trtVal)]


def extractSeq(df, nanLabel, var, index, w, isTargetVar, length):
    s = int(w-(length-1))
    e = int(w+1)
    # if isTargetVar:
    #   e = int(w)
    cols= []
    for k in range(s,e):
        col = "{}_n_{}".format(var,k)
        cols.append(col)

    seq = np.array(df.loc[index, cols])
    seqLabel = np.array(nanLabel.loc[index, cols])
    return (seq, seqLabel)


def evaluateLValues(distances, df, nanLabel):
    for pair in distances:
        print "\n"
        # distanceInfoT = measureSimilarity(var, getTreatmentSignal(), df, nanLabel)
        # distanceInfoC = measureSimilarity(var, getControlSignal(), df, nanLabel)
        
        # DT= distanceInfoT[:,2]
        # DC= distanceInfoC[:,2]
        DT =  pair[0]
        DC = pair[1]
        outliersIndexT, LT = tuneL(DT)
        outliersIndexC, LC = tuneL(DC)

        print "Control:{}".format(len(outliersIndexC))
        print "Treatment:{}".format(len(outliersIndexT))

def tuneL(Data, nanLabel, distanceInfo, var, string):
    isKnown = []
    for i in range(0,len(Data)):
        index = int(distanceInfo[i,0])
        w = int(distanceInfo[i,1])
        col= "memIndex_{}".format(w)
        isKnown.append(not nanLabel.loc[index,col])
    isKnown = np.array(isKnown)
    

    UPPER_LIMIT=500  #to be fixed
    LOWER_LIMIT=100
    Data[Data==0] = np.partition(np.unique(Data),1)[1]/10.0
    Data = np.log(Data)
    Data.tofile("{}_distance_{}.csv".format(string, var), sep=',')
    isKnown.tofile("{}_isMIKnown_{}.csv".format(string, var), sep=',')
    L=1
    Median = np.median(Data[np.where(isKnown)])
    MADVal = MAD(Data[np.where(isKnown)])
    if MADVal==0:
       MADVal=0.01
    print "MAD:{}".format(MADVal)
    #print "1:{}".format(L)
    while (True):
        Outliers = (Data < (Median - L*MADVal))
        Outliers =  np.logical_and(Outliers, isKnown) 
        outliersIndex = np.where(Outliers)[0]
        size= len(outliersIndex)
        # print "size:{}, L:{}".format(size,L)
        if(size<LOWER_LIMIT): #modified - needs to be checked. Prevs: UPPER_LIMIT
            break
        L=L+1
        #print "2:{}".format(L)
    while (True):
        Outliers = (Data < (Median - L*MADVal))
        Outliers =  np.logical_and(Outliers, isKnown) 
        outliersIndex = np.where(Outliers)[0]
        size= len(outliersIndex)
        # print "size:{}, L:{}".format(size,L)
        if(size>LOWER_LIMIT):
            break
        L=L-1
        #print "3:{}".format(L)
    if (len(outliersIndex)>UPPER_LIMIT):
        idx = np.argpartition(Data[outliersIndex], UPPER_LIMIT)[:UPPER_LIMIT]
        outliersIndex = outliersIndex[idx]
        #outliersIndex = np.random.choice(outliersIndex, UPPER_LIMIT, replace=False)
    passedIndexes = [index for  index in outliersIndex if Data[index]< -2 ]
    print "Size before pruning:{}".format(len(outliersIndex))
    print "Final size after pruning:{}".format(len(passedIndexes))
    #print Data[outliersIndex]
    return (passedIndexes, L)


def getVarDistances(df, nanLabel):
    result= []
    for var in trtmntVar:
        distanceInfoT = measureSimilarity(var, getTreatmentSignal(), df, nanLabel)
        distanceInfoC = measureSimilarity(var, getControlSignal(), df, nanLabel)
        
        DT= distanceInfoT[:,2]
        DC= distanceInfoC[:,2]

        result.append((DT,DC))
    return result



def drawKmeanDiagram(Data):
    Data= Data.reshape(-1, 1)
    Sum_of_squared_distances = []
    K = range(1,15)
    for k in K:
        km = KMeans(n_clusters=k)
        km = km.fit(Data)
        Sum_of_squared_distances.append(km.inertia_)

    plt.plot(K, Sum_of_squared_distances, 'bx-')
    plt.xlabel('k')
    plt.ylabel('Sum_of_squared_distances')
    plt.title('Elbow Method For Optimal k')
    plt.show()


def computeDiff(ids, df):
    diffs = []
    for pair in ids:
        index, w= pair.rstrip().split(" ")
        index = int(index)
        w= int(w)
        col = "memIndex_{}".format(w)
        col_prev = "memIndex_{}".format(w-1)
        diffs.append(df.loc[index,col]-df.loc[index,col_prev])
    return diffs



def fetchLEHyps(path):
    hyps = genfromtxt(path, delimiter=',')
    hyps = hyps[1:,:-1] 
    LEHyps = set()
    for i in range(len(hyps)):
        seq = hyps[i,:].astype(int)
        print seq
        if (seq==1).any():
            strSeq= str(seq)
            pattern = re.sub("[^0-9]", "",strSeq)
            LEHyps.add(pattern)
    print LEHyps
    return LEHyps

def splitStr(word): 
    return [int(char) for char in word]

def array2id(array):
    strSeq= str(array.astype(int))
    pattern = re.sub("[^0-9]", "",strSeq)
    return pattern

def id2array(id):
    charList = splitStr(id)
    return np.array(charList)

def getNeighbours(h):
    neighbours= set()
    for i in range(0,len(h)):
        for j in [0,1,2]:
            if int(h[i])!= j:
                t= id2array(h)
                t[i]=j
                neighbours.add(array2id(t))
    return neighbours

# def evaluate(var, trtSeq):
#     pvals = {
# "001": 0.0940618851493,
# "010": 0.220666932262,
# "011": 4.56538888534e-06,
# "012": 1.36081090715e-05,
# "021": 5.20254263439e-08,
# "100": 0.599608704212,
# "101": 0.000422794877288,
# "102": 0.0282634922079,
# "110": 0.00457976459746,
# "111": 0.000241111600374,
# "112": 2.50396924891e-05,
# "120": 0.212761805018,
# "121": 0.000613897909022,
# "122": 0.00123168476547,
# "201": 0.00616095632461,
# "210": 0.00219695225111,
# "211": 5.26449229881e-06,
# "212": 1.31416673629e-05,
# "221": 0.00144933730433,
#     }
#     return pvals[trtSeq]


cache={}

def evaluate_RBD(var, trtSeq):
    if trtSeq in cache:
        return cache[trtSeq]

    if (os.path.isfile(dfPath) and os.path.isfile(nanLabelPath)):
        df = pd.read_pickle(dfPath)
        nanLabel = pd.read_pickle(nanLabelPath)
    else:
        df = readData()
        df, nanLabel = preprocess(df)

    trtSeq = id2array(trtSeq)
    weights = (~(trtSeq==2)).astype(int)
    trtSignal = (trtSeq,weights)                
    distanceInfoT = measureSimilarity(var, trtSignal, df, nanLabel)
    outliersIndexT = detectOutliers(distanceInfoT, nanLabel, var, "Treatment")

    anchorPoint = (np.where(trtSeq!=2)[0][0])
    anchorDist = len(trtSeq)- anchorPoint
    ctrlSeq = np.zeros(anchorDist)
    ctrlWeights = np.ones(anchorDist)
    ctrlSignal = (ctrlSeq, ctrlWeights)

    distanceInfoC = measureSimilarity(var, ctrlSignal, df, nanLabel)
    outliersIndexC = detectOutliers(distanceInfoC, nanLabel, var, "Control")

    if (len(outliersIndexC)==0 or  len(outliersIndexT)==0 ):
        with open("searchResult.txt","a") as f:
            f.write("{0} pattern: {1} , {2}\n".format(var, trtSeq.astype(int), "NA - outlier detection returned zero samples"))
        cache[array2id(trtSeq)] = np.nan
        return np.nan

    #C = computeDistanceMatrix2(df, nanLabel, var, outliersIndexT, outliersIndexC, distanceInfoT, distanceInfoC, trtSeq)
    C = computeDistanceMatrix2_RBD(df, nanLabel, var, outliersIndexT, outliersIndexC, distanceInfoT, distanceInfoC, trtSeq)
    matchedPairs = performMatching_RBD(C, len(outliersIndexT))
    if (len(matchedPairs)<4):
        with open("searchResult.txt", "a") as f:
            f.write("{0} pattern: {1} , {2}\n".format(var, trtSeq.astype(int), "NA - matching returned less than four samples"))
        cache[array2id(trtSeq)] = np.nan
        return np.nan             

    #[isBiased, meanVals] = isDCBiased(df, matchedPairs, outliersIndexT, outliersIndexC,distanceInfoT, distanceInfoC, var, trtSeq)
    
    anchorPoint = (np.where(trtSeq==1)[0][0]-1)
    if (anchorPoint<0):
        anchorPoint=0
    anchorDist = len(trtSeq)- anchorPoint

    targetValues = extractTargetValues(df, matchedPairs, outliersIndexT, outliersIndexC,distanceInfoT, distanceInfoC, var, anchorDist)
    pval = computePValue(targetValues[0], targetValues[1])
    with open("searchResult.txt","a") as f:
        meanValsStr = str(meanVals)
        f.write("{0} pattern: {1}, pval={2:}, ACE={4: .2f}, n={3:d}, DCT Mean={5}\n".format(var, trtSeq.astype(int), pval, len(matchedPairs), np.mean(targetValues[1])- np.mean(targetValues[0]),meanValsStr))
    
    
    # if (isBiased):
    #     cache[array2id(trtSeq)]= np.nan
    #     return np.nan
    # else:
    cache[array2id(trtSeq)]= pval
    return pval



def evaluate(var, trtSeq):
    if trtSeq in cache:
        return cache[trtSeq]

    if (os.path.isfile(dfPath) and os.path.isfile(nanLabelPath)):
        df = pd.read_pickle(dfPath)
        nanLabel = pd.read_pickle(nanLabelPath)
    else:
        df = readData()
        df, nanLabel = preprocess(df)

    trtSeq = id2array(trtSeq)
    weights = (~(trtSeq==2)).astype(int)
    trtSignal = (trtSeq,weights)                
    distanceInfoT = measureSimilarity(var, trtSignal, df, nanLabel)
    outliersIndexT = detectOutliers(distanceInfoT, nanLabel, var, "Treatment")

    anchorPoint = (np.where(trtSeq!=2)[0][0])
    anchorDist = len(trtSeq)- anchorPoint
    ctrlSeq = np.zeros(anchorDist)
    ctrlWeights = np.ones(anchorDist)
    ctrlSignal = (ctrlSeq, ctrlWeights)

    distanceInfoC = measureSimilarity(var, ctrlSignal, df, nanLabel)
    outliersIndexC = detectOutliers(distanceInfoC, nanLabel, var, "Control")

    if (len(outliersIndexC)==0 or  len(outliersIndexT)==0 ):
        with open("searchResult.txt","a") as f:
            f.write("{0} pattern: {1} , {2}\n".format(var, trtSeq.astype(int), "NA - outlier detection returned zero samples"))
        cache[array2id(trtSeq)] = np.nan
        return np.nan

    C = computeDistanceMatrix2(df, nanLabel, var, outliersIndexT, outliersIndexC, distanceInfoT, distanceInfoC, trtSeq)
    matchedPairs = performMatching(C)
    if (len(matchedPairs)<4):
        with open("searchResult.txt", "a") as f:
            f.write("{0} pattern: {1} , {2}\n".format(var, trtSeq.astype(int), "NA - matching returned less than four samples"))
        cache[array2id(trtSeq)] = np.nan
        return np.nan             

    [isBiased, meanVals] = isDCBiased(df, matchedPairs, outliersIndexT, outliersIndexC,distanceInfoT, distanceInfoC, var, trtSeq)
    targetValues = extractTargetValues(df, matchedPairs, outliersIndexT, outliersIndexC,distanceInfoT, distanceInfoC, var)
    pval = computePValue(targetValues[0], targetValues[1])
    with open("searchResult.txt","a") as f:
        meanValsStr = str(meanVals)
        f.write("{0} pattern: {1}, pval={2:}, ACE={4: .2f}, n={3:d}, DCT Mean={5}\n".format(var, trtSeq.astype(int), pval, len(matchedPairs), np.mean(targetValues[1])- np.mean(targetValues[0]),meanValsStr))
    
    
    if (isBiased):
        cache[array2id(trtSeq)]= np.nan
        return np.nan
    else:
        cache[array2id(trtSeq)]= pval
        return pval



def findMin(U, pVals):
    minVal = 2
    minID = -1
    for hypID in U:
        if(pVals[hypID]<minVal):
            minVal = pVals[hypID]
            minID = hypID
    return (minID, minVal)  


def search(var, s, LowE_Path):
    # var="heactb"
    # s = "212"
    pVals={}
    prev = {}
    U = fetchLEHyps(LowE_Path)
    print U
    for hypId in U:
        pVals[hypId] = 2
        prev[hypId] = "nan"

    pVals[s]=evaluate_RBD(var, s)
    
    bestSoFarVal =  pVals[s]
    bestSoFarID = s

    print pVals

    print bestSoFarID
    print bestSoFarVal


    counter= 0
    while(len(U)>0):
        print "it : {}".format(counter)
        minID, minVal=findMin(U, pVals)
        
        
        if (minVal >1):
            break
        print "min ID:{}, min value:{}, bestSoFarVal:{}".format(minID, minVal, bestSoFarVal)
        if (minVal<=bestSoFarVal):
            print "best so far changed"
            bestSoFarID= minID
            bestSoFarVal = minVal
        else:
            print "goint to break"
            break

        for v in getNeighbours(minID):
            if (v in U):
                pVals[v] = evaluate(var, v)
                prev[v] = minID
        U.remove(minID)
        counter=counter+1

    print bestSoFarID
    print bestSoFarVal
    print pVals
    print prev


    printPath(bestSoFarID, prev)
    with open('pVals.txt', 'w') as file:
        file.write(pickle.dumps(pVals)) 
    
    with open('prev.txt', 'w') as file:
        file.write(pickle.dumps(prev)) 

    return (pVals, prev, bestSoFarID)


def printPath(node, prev):
    print "prev[node]:",prev[node]
    if ( prev[node]=="nan"):
        print node,
    else:
        printPath(prev[node], prev)
        print "->",
        print node,

def getPvalStats(pVals):
    calclulated=0
    pruned = 0
    notReached=0
    for node, pval in pVals.items():
        if (pval==np.nan):
            pruned=pruned+1
        elif(pval==2):
            notReached=notReached+1
        else:
            calclulated= calclulated+1

    return (calclulated, pruned, notReached)   


if __name__ == "__main__":
    print "a"


