# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 19:46:51 2023

@author: Administrator
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem import Draw

pth = r'D:\纤维素裂解\数据处理\产物\产物\\'
f = r'2800.spec.tab_sps.csv'
df = pd.read_csv(pth+f,index_col=0)
time = df['Time']
df = df.loc[:,~df.columns.isin(["Time"])]

# df = df.loc[:, lambda d: d.columns.str.contains('C')]
# df.filter(regex='C')

def replaing(val) :
    try:
        num = int(val.replace('C',''))
    except :
        num = 1
    return num

def varing(sm) :
    try :
        mol = Chem.MolFromSmiles(sm)
        formula = CalcMolFormula(mol)
        num=-100
    except:
        num = 100
    return num

C_list = []

for sps in df.columns :

    if 'C' in sps : 
        if varing(df.loc[ 0 ,sps]) != 100 :
            
        # sps =   'CO'
    # CO
            if 'O' in sps and ('H' not in sps) :
                c_num = sps.split('O')
                C_list.append(replaing(c_num[0]))
    # CH/CHO    
            else :
                c_num = sps.split('H')
                C_list.append(replaing(c_num[0]))
        else :
            C_list.append(100)
    # 不含碳        
    else :
        C_list.append(0)


df.loc[len(df.index)] = C_list
df = df.loc[2:,]
C_gas = df[df.columns[df.iloc[-1] == 0 ]].astype(int) # 1-4 
other_gas = df[df.columns[df.iloc[-1] > 0 ] & df.columns[df.iloc[-1] < 5 ]].astype(int)      # 1-4 
light = df[df.columns[df.iloc[-1] > 4 ] & df.columns[df.iloc[-1] < 14 ]  ].astype(int)       # 5-13
heavy = df[df.columns[df.iloc[-1] > 13 ] & df.columns[df.iloc[-1] < 32 ]  ].astype(int)      # 14-31
coke = df[df.columns[df.iloc[-1] > 31 ]                                 ].astype(int)        # > 32

catagary = pd.DataFrame()
catagary['Time'] = time[2:]
catagary['C_gas'] = C_gas.apply(lambda x: x.sum(), axis=1)
catagary['other_gas'] = other_gas.apply(lambda x: x.sum(), axis=1)
catagary['light'] = light.apply(lambda x: x.sum(), axis=1)
catagary['heavy'] = heavy.apply(lambda x: x.sum(), axis=1)
catagary['coke'] = coke.apply(lambda x: x.sum(), axis=1)     

catagary.to_csv(pth+f.replace('.spec.tab_sps.csv','_')+'category.csv')

        
        


