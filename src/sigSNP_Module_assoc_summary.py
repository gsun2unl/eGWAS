from scipy.stats import fisher_exact,chi2_contingency
import os
import pandas as pd
import numpy as np

dftas_peak = pd.read_csv('GoodPeakSummary.csv', index_col = 'IC')
ICList = dftas_peak.index.unique()
dfecis = pd.read_csv('TAS_cis.out',sep = '\t', index_col = 0)
dfetrans = pd.read_csv('TAS_Trans.out',sep = '\t', index_col = 0)
dfetas = pd.concat([dfecis, dfetrans])
dfetas['logP'] = -np.log(dfetas['p-value'])
dfetas = dfetas.loc[dfetas['logP']>=5]
nameIC, pVal, pNum = [], [], []
fpout = open('eTraits_peaks_enrichement.txt','w')
headList = ['leadSNP', 'eTraitNumber','ModuleSize', 'EnrichmentFactor' ,'chi2_Pval', 'IC', 'Peak']
fpout.write('\t'.join(headList) + '\n')

for IC in ICList:
    dfTAS_IC = dftas_peak.loc[IC]
    IC_Gene = open('ICA/{0}.list'.format(IC)).readlines()[1:]
    IC_Gene = [agene.strip() for agene in IC_Gene]
    dfTAS_IC.set_index('Peak', inplace=True)
    for apeak in dfTAS_IC.index.unique():
        dfTAS_IC_Peak = dfTAS_IC.loc[apeak].sort_values(by='logPvalue')
        peakSNP = dfTAS_IC_Peak['rs']
        leadSNP = map(str,list(dfTAS_IC_Peak.iloc[0,[1,3]].values))
        leadSNP = ':'.join(leadSNP)
        rsList = [rs for rs in peakSNP if rs in dfetas.index] #note that the lead SNP might not have significant associations with expression of genes while the SNPs in the peak with lower logP than the
        #lead SNP could have more associations with expression of genes. This lead to the SNP you labeled in this graph not the same as the lead SNP in the peak. 
        if len(rsList)==0:continue
        dfrs = dfetas.loc[set(rsList)].reset_index().set_index('gene')
        meGene = [agene for agene in dfrs.index if agene in IC_Gene]
        meGene = list(set(meGene))
        me = len(meGene) 
        mne = len(IC_Gene) - me
        nme = len(dfrs.index.unique()) - me
        nmne = 19585 - len(IC_Gene) - nme
        #print(me, mne, nme, nmne)
        chi2, newp, dof, ex = chi2_contingency([[me,mne],[nme, nmne]])
        if newp <= (0.01/len(dfTAS_IC.index)): 
            em = me/len(IC_Gene)
            eb = len(dfrs.index.unique())/19565
            enrichment = em/eb
            if enrichment > 1:
                leadSNP = leadSNP + '(' + IC+ ')'
                peakName = 'Peak'+str(apeak)
                nameIC.append(leadSNP)
                pVal.append(newp)
                pNum.append(me)
                print('etrait enriched in ', leadSNP + '(' + IC+ ')' , 'Peak'+str(apeak))
                print(len(IC_Gene), 'module genes; ', me, 'etraits in module gene',  )
                #fout = open(IC + '_peak'+str(apeak)+'_gene.txt','w')
                lineList = [leadSNP, str(me), str(len(IC_Gene)), str(enrichment), str(newp), IC, str(apeak)]
                fpout.write('\t'.join(lineList) + '\n')
                #for agene in meGene:
                #    fout.write(agene + '\n')
                #fout.close()
fpout.close()               

