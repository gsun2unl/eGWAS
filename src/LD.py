import pandas as pd
import os

dftas = pd.read_csv('sigGWAS_RES/consolidatedPeaks_peaks.txt')
dfsnp = pd.read_csv('sigGWAS_RES/V1_transformed.assoc.txt', sep = '\t', usecols = ['chr','rs','ps'], index_col = 'rs')

sigList = dftas.rs.unique()
sList = []

dfsig = dfsnp.loc[sigList]
dfsig.to_csv('sigSites.list', sep = '\t', index=False)
for asit in sigList:
    print(asit)
    ssList = []
    sChr = dfsnp.loc[asit]['chr']
    sPos = dfsnp.loc[asit]['ps']
    sst = sPos - 1000000
    ssd = sPos + 1000000 
    df = dfsnp.loc[dfsnp['chr'] == sChr] 
    df = df.loc[df['ps'] >= sst]
    df = df.loc[df['ps'] <= ssd]
    ssList = list(df.index)
    sList += ssList 

sList = list(set(sList))
sList.sort()
print(len(sList))

fout = open('Allsites.list','w')
for x in sList:
    fout.write(x + '\n')
    
fout.close()


#vcftools --gzvcf /home/gsun2unl/Documents/ePscoR_Schnablelab/Phenotype/GEMMA/e340.vcf.gz --snps Allsites.list --recode --stdout | bgzip -c > AllSites.vcf.gz
#vcftools --gzvcf AllSites.vcf.gz --geno-r2-positions sigSites.list --out AllSiteLD
