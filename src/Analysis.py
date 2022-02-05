#!/usr/bin/env python
# coding: utf-8


import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from collections import Counter


dfeSum = pd.read_csv('e340_eSummary.csv',index_col=0)
dfeSum = dfeSum.loc[dfeSum.index.value_counts()<=10]
dfeSum = dfeSum.loc[dfeSum.index.unique().sort_values()]

dfHb = pd.read_csv("eTrait_Hb_PVE.csv",index_col = 0)

####### Panel B

fig, ax = plt.subplots(figsize = (4,4), tight_layout=True)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.tick_params(axis='both', which='minor', labelsize=8)
#sckw = {'s':8,'alpha':0.7, }
#line_kws={"color": "red", 'alpha':0.7}

ax.scatter(dfHb.Hb,dfHb.pve,color = 'black', s=5)
x = np.linspace(0,1,100)
y = np.linspace(0,1,100)
ax.plot(x,y, color = 'red', linewidth=1.5, alpha=0.5)
plt.xlabel('Broad sense heritability', fontsize=12)
plt.ylabel('eQTL PVE', fontsize=12)
plt.savefig('../Plots/Hb_PVE_scatter.png', dpi=300)
plt.savefig('../Plots/Hb_PVE_scatter.svg')



####### Panel A #########

Color = {'Cis':'orange', 'Trans':'blue'}
fig, ax = plt.subplots(figsize= (3,4), tight_layout=True)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.tick_params(axis='both', which='minor', labelsize=8)

for idx, grp in dfeSum.groupby('Type'):
    print (idx)
    color = Color[idx]
    #ne,xe,_e = ax.hist(grp['PVE'],bins = 100, color=color, density = True, alpha = 0.5)
    #bin_centers_e = 0.5*(xe[1:]+xe[:-1])
    #ax.plot(bin_centers_e,ne)
    sns.distplot(grp['pve'], hist = False, kde=True, kde_kws = {'shade':True, 'linewidth':2}, label = 'eType')
    ax.set_xlim(0,1)
    
ax.set_xlabel('eQTL PVE', fontsize=12)
ax.set_ylabel('Density', fontsize=12)
plt.legend(prop = {'size':12}, labels = ['Cis eQTLs','Trans eQTLs'])
plt.savefig('../Plots/PVE_DensityPlot.png', dpi = 200)
plt.savefig('../Plots/PVE_DensityPlot.svg')


########## Panel C ###############

fig, ax = plt.subplots(figsize = (4,4), tight_layout=True)
ax.tick_params(axis='both', which='major', labelsize=12)
#ax.tick_params(axis='both', which='minor', labelsize=8)

dfPlot = dfHb.loc[dfHb.prop<=1]
dfPlot['prop'].hist(bins = 20, edgecolor = 'black', linewidth=1.2, color = 'grey', ax = ax, grid = False)
plt.xlabel('Fraction of heritability explained',fontsize=12)
plt.ylabel('Number of genes', fontsize=12)
plt.savefig('../Plots/FractionOfHbExplained.png', dpi = 300)
plt.savefig('../Plots/FractionOfHbExplained.svg')







