import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

dftrans340 = pd.read_csv('e340trans_afexp.csv', index_col = 0 )
dftrans140 = pd.read_csv('e140trans_afexp.csv', index_col = 0 )


x340 = np.array(dftrans340.index)
y340 = np.array(dftrans340['e340mean'])
y340_perm = np.array(dftrans340['e340meanPerm'])

x140 = np.array(dftrans140.index)
y140 = np.array(dftrans140['e140mean'])
y140_perm = np.array(dftrans140['e140meanPerm'])

fig, (ax1, ax2) = plt.subplots(1,2, figsize=(8,4), tight_layout=True)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

ax1.set_title('340 inbred lines')
ax2.set_title('200 inbred lines')

ax1.plot(x340, y340, '.',color = 'b', label='Before_Permutation') 
ax1.plot(x340,y340_perm,'.', color = 'r', label='After_Permutation')

ax2.plot(x140, y140, '.' , color = 'b', label = 'Before permutation') 
ax2.plot(x140, y140_perm, '.' , color = 'r', label = 'After permutation')

ax1.set_xlabel('Low      Expression rank      High')
ax2.set_xlabel('Low      Expression rank      High')
ax1.set_ylabel('Trans rare allele count')
ax1.legend(loc='best', fontsize = 10)
ax2.legend(loc='best', fontsize = 10)
#plt.savefig('../../eGWAS340_Plots/Top5000GeneRank.png', dpi=300)
