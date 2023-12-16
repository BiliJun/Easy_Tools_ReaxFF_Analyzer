# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 20:54:08 2023

@author: zh
"""
import pandas as pd

pfn = pd.DataFrame()

import os

# f = r'D:\zh\JP-10\REAXFF\JP-10-2200\bond\bondanalysis\\'

# pth = [p for p in os.listdir(f) if 'out' not in p]

# for p in pth :
#     pf = pd.read_csv(f+p+'\\Default.spec.tab_sps.csv')
#     pfn['t'+p] = pf['Time']
#     pfn[p] = pf['C10H16']
    
# pfn = pfn.dropna()

# pfn.to_csv('D:\zh\JP-10\结果\\C10H16.csv')

# ---------------------

import numpy as np
import pandas as pd
import math  
    
timstep = 10**(-12)   # ps -> s


df = pd.read_csv(r'D:\ActiveEnergy\CH4.csv' , index_col=0 )         
ppf = pd.DataFrame()  # 平均值列表
rateK = pd.DataFrame()  # ln(k) 列表
df = df.iloc[1:].astype(float)   
nm = df.columns

for i in range(0, len(nm), 2) :
    pass
    a = pd.DataFrame()
    a['0'] = df[nm[i]]
    a['1'] = df[nm[i+1]]
    a = a.dropna(axis=0,how='all')
    
    # ------
    k = (a['1'].apply(np.log) - math.log(a.loc[1,'1']) ) / a['0']/timstep
    k = k.dropna(axis=0,how='all').reset_index(drop=True)
    kave = abs(k.mean())
    rateK[nm[i+1]] = [kave, math.log(kave) ]
    ppf[nm[i+1]] = k
rateK.to_csv(r'D:\ActiveEnergy\rate2.csv')


