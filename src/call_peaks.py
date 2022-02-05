import os, sys
import pandas as pd
import numpy as np


def uniquePeak (df):
    idxList = []
    peaks = []
    pNumber = 0
    window = 1e6
    for idx, dfgrp in df.groupby('chr'):
        pNumber += 1
        peakList = []
        dfgrp.sort_values(by = 'ps', inplace=True)
        idxList += list(dfgrp.index)
       # print('number of sig snps in chr', idx, 'is ', len(dfgrp))
        sigSNPList = list(dfgrp.ps)
        if len(sigSNPList) == 1: 
            peaks.append(pNumber)
        elif len(sigSNPList) == 2: 
            if sigSNPList[1] - sigSNPList[0] > window: 
                peaks.append(pNumber)
                pNumber +=1 
                peaks.append(pNumber)
            else: peaks += [pNumber, pNumber]
        else:
            lastpos = sigSNPList[0]
            for x in sigSNPList: 
                if x - lastpos <= window:
                    #print('chr', idx, x, pNumber)
                    lastpos = x
                    peaks.append(pNumber)
                else: 
                    pNumber +=1
                    peaks.append(pNumber)
                    lastpos = x
    df = df.loc[idxList]
    df['peaks'] = peaks
    return(df)
 
indDir = sys.argv[1]
peakDir = sys.argv[2] 
    
SigSNP =[]

icList = [x.split('.')[0] for x in os.listdir('ICA/')] 
#eFList = [os.path.join(indDir,IC + '_transformed_signal.csv') for IC in icList]

dfconsolidate = pd.DataFrame()

cutoff = -np.log10(0.05/12191984)

for i, ic in enumerate(icList):
    F_signal = os.path.join(indDir, ic+ '_GWAS_signal.csv')
    if not os.path.exists(F_signal) : continue
    outName = peakDir + '/' + ic + '.summary.csv'
    print(i, ic)
    df = pd.read_csv(F_signal, index_col = 0)
    pcol = list(df)[-1]
    df = df.loc[df[pcol] >=cutoff]
    if len(df) == 0: continue
    df = uniquePeak(df)
    #if len(df.peaks.unique())>5: continue
    #print(df.head())
    dfpeak = pd.DataFrame() 
    for idx, grp in df.groupby('peaks'):
        st = grp['ps'].min()
        ed = grp['ps'].max()
        distance = ed - st
        num = len(grp)
        grp = grp.loc[[grp[pcol].idxmax()]]
        grp['SnpNumer'] = [num]
        grp['pStart'] = st
        grp['pStop'] = ed
        grp['pLength'] = distance  
        #grp.index = [x]
        #print(grp)
        dfpeak = dfpeak.append(grp)
    SigSNP += list(dfpeak['rs'].values)
    dfpeak.to_csv(outName, index= False)
    dfpeak['PC'] = ic
    dfconsolidate = pd.concat([dfconsolidate, dfpeak])

dfconsolidate.reset_index(inplace=True)
df = uniquePeak(dfconsolidate)

df = df.loc[df['SnpNumer']>=5]
df = df.iloc[:, 1:]


df.to_csv(peakDir +'/consolidatedPeaks_peaks.txt', index = False)
#dfconsolidate.to_csv('consolidatedPeaks.csv')






