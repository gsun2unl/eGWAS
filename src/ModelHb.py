import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

dfhb = pd.read_csv('ModuleHeritability.csv', index_col = 0)
dfhb.sort_values(by = 'Heritability', inplace=True, ascending=False)

from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.size'] = 12

fig, ax =plt.subplots(figsize=(4,3), tight_layout = True)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=12,direction='out', length=4, width=2)

xpos = list(range(1,len(dfhb)+1))
ax.bar(xpos, dfhb.Heritability, color = 'purple', alpha=0.6, width = 1)
ax.set_xlim(0,165)
ax.set_ylabel('Broad sense heritability')
ax.set_xlabel('Models')
plt.savefig('ModelHeritability.svg')
plt.savefig('ModelHeritability.png', dpi = 500)
