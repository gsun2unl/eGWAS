import os
import pandas as pd
import numpy as np
dftasnum = pd.read_csv('TAS_In_Good_Peaks.numeric.xmat',sep = '\t', index_col = 0)
dftasnum.transpose().to_csv('TAS_In_GoodPeaks.numeric.txt',sep = '\t')
dftasmap = pd.read_csv('TAS_In_Good_Peaks.numeric.map',sep = '\t', header = None, usecols = [1,0,3])
dftasmap = dftasmap[[1,0,3]]
dftasmap.to_csv('TAS_In_GoodPeaks_snpsloc.txt', sep = '\t', index = False, header = None)
dfcvrt = pd.read_csv('e340.Covariates',sep = '\t', header = None, index_col = 0)
dfcvrt.transpose().to_csv('e340.covar',sep='\t')
