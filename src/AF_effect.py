import os
import subprocess as sp
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
%matplotlib inline


dfePlot = pd.read_csv('AF_effect.csv')
fig, axes = plt.subplots(1,3, figsize=(10,3), tight_layout = True)
colorDict = {'Cis':'lightblue', 'Trans':'orange'}
for etype, grp in dfePlot.groupby('Type'):
    color = colorDict[etype]
    print(etype, color)
    ax1 = axes[0]
    sns.distplot(grp['maf'], bins=30, hist = False, kde=True, color = color, kde_kws = {'shade':True, 'linewidth':2}, label = 'eTraitType', ax = ax1)
    ax1.set_xlim(0,0.5)
    ax1.legend(['Cis eQTL', 'Trans eQTL'])
    ax1.set_xlabel('Minor allele frequency')
    ax1.set_ylabel('Density')
    
    ax2 = axes[1]
    sns.regplot(x = 'bin', y='R2', data = grp, x_estimator=np.mean, color = color, ax = ax2)
    ax2.set_xlabel('Minor allele frequency')
    ax2.set_ylabel('R2')
    
    ax3 = axes[2]
    sns.regplot(x = 'bin', y='beta', data = grp, x_estimator=np.mean, color = color, ax = ax3)
    ax3.set_xlabel('Minor allele frequency')
    ax3.set_ylabel('Effect size')
    #ax2.legend(['Cis eQTL', 'Trans eQTL'])
    #ax3.legend(['Cis eQTL', 'Trans eQTL'])
plt.savefig('eQTL_AlleFreq_Effect.png', dpi = 300)
plt.savefig('eQTL_AlleFreq_Effect.svg')
