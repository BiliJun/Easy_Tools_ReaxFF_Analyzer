# -*- coding: utf-8 -*-
"""
Created on Sun Nov  5 10:31:23 2023

@author: zhouw
"""


fpth = r'D:\zh\cell-pe\数据\2400\\'
fnms = 'Default.spec.tab_sps.csv'


from rdkit import Chem
import pandas as pd
import tqdm
# Analysing bond information
def bond_information (smi) : 
    # smi = '[H][H]'
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.RemoveHs(mol)
    if mol is not None:        
        '''# 获取分子中的原子的相邻原子的数目以及序号
        # 调用 GetNumAtoms 函数, #获取分子中的原子数目
        atom_num = mol.GetNumAtoms()
        nei_atom = []
        nei_atomEle = []
        for i in range(atom_num):
           nei_atom.append([(i.GetSmarts(),i.GetIdx()) for i in mol.GetAtomWithIdx(i).GetNeighbors()])
           nei_atomEle.append([i.GetSmarts() for i in mol.GetAtomWithIdx(i).GetNeighbors()])
        '''
        #获取分子中键的相关类型以及键的特征，以及成键原子序号
        # nei_bond = []
        nei_ba = []
        #获取分子中的键数目
        bond_num = mol.GetNumBonds()
        for i in range(bond_num):
           # nei_bond.append((mol.GetBondWithIdx(i).GetBondType().name,mol.GetBondWithIdx(i).GetBeginAtomIdx(),mol.GetBondWithIdx(i).GetEndAtomIdx()))
           nei_ba.append((mol.GetBondWithIdx(i).GetBondType().name,mol.GetBondWithIdx(i).GetBeginAtom().GetSmarts(),mol.GetBondWithIdx(i).GetEndAtom().GetSmarts()))
        
        bdlst = pd.DataFrame(nei_ba,columns=['B_type','A1','A2']) 
        bdlst.replace('-SINGLE-','-', inplace = True)
        bdlst.replace('-DOUBLE-','=', inplace = True)
        bdlst.replace('(H\d)|\[|\]|@','', regex = True, inplace = True)
        bdlst.replace('H','', regex = True, inplace = True)
        bdlst['FULLtype'] = bdlst['A1'].str.upper() +'-'+ bdlst['B_type'] + '-'+ bdlst['A2'].str.upper()
        x = bdlst['FULLtype'].value_counts()
        x = pd.DataFrame(x)
        x.columns = [smi]
        x = x.to_dict('dict')
        y = list( ii for i in x.keys() for ii in x[i] )
    else :
        x = [], y = []
    return x,y



pf = pd.read_csv(fpth+fnms, header = 1,index_col = 0)  # read spes.csv
pf = pf.reset_index(drop=True) 
sps_lst = [ii for ii in list(pf.columns) ]   # !!  要检索的产物  spes simles 
bd_type, bond_lst =[], {}   # Summary of sps bond type and number list

for i_sps in sps_lst :  # 产物遍历
    try :
        bd_sps, b_nms = bond_information (i_sps)
        bond_lst.update( bd_sps )
        bd_type.extend(b_nms)
    except :
        print(i_sps)
        
bond_lst = pd.DataFrame(bond_lst)
bond_lst = bond_lst.fillna(0)            # 产物的键数列表
bd_type = list(set(bd_type))             # 键类型
BOND_DATA = pd.DataFrame ( )
BOND_DATA['Time'] = pf['0.0']


for ib in bd_type :
    a = [0 for i in range(len(pf))]
    for isp in sps_lst :
        try :
            b = bond_lst[isp][ib]*pf[isp]
            b.tolist() 
            a = a + b
        except :
            print('drop this col in .csv file', isp)
            
    BOND_DATA [ib] = a
BOND_DATA.to_csv(fpth+fnms.replace('.csv','_bond_type.csv'))
        
   
    

