import os
import pandas as pd
import numpy as np
dftas = pd.read_csv('sigGWAS_RES/consolidatedPeaks_peaks.txt')
dftas_peak = pd.DataFrame()
for afile in os.listdir('ICA_GWAS_Signal/'):
    if afile.endswith('.summary.csv'):
        IC = afile.split('.')[0]
        dfsignal = pd.read_csv('ICA_GWAS_Signal/{0}_GWAS_signal.csv'.format(IC))
        dfpeak = pd.read_csv('ICA_GWAS_Signal/{0}'.format(afile), index_col = 'peaks')
        dfpeak = dfpeak.loc[dfpeak['SnpNumer'] >=3]
        for apeak in dfpeak.index:
            pChr = dfpeak.loc[apeak,'chr']
            pST = dfpeak.loc[apeak,'pStart']
            pED = dfpeak.loc[apeak, 'pStop']
            dfsig = dfsignal.loc[dfsignal.chr == pChr]
            dfsig = dfsig.loc[(dfsig.ps >= pST) & (dfsig.ps <= pED)]
            dfsig['IC'] = IC
            dfsig['Peak'] = apeak
            dftas_peak = dftas_peak.append(dfsig)
dftas_peak.iloc[:, 1:].to_csv('GoodPeakSummary.csv') #outputing all of the significant SNPs within the consolidated peaks. 
dftas_peak.rs.to_csv('TAS_In_GoodPeaks.list', header = None, index = False)
