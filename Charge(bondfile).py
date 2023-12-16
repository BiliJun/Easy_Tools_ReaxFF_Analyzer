#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 22:39:09 2022

@author: zh
"""


# 原子编号-关联原子列表，判断：列表中的原子是否为H,O，进而判断与原子形成的键。归纳到键列表NBOND中，去除重复的项。



# bonding atom list generate

file_names ='bonds.2000.out'    
ele_n = 3
Atom_type = ['C','H','O']
timestep = 10000 # 0.1 fs
skip = 50

#from numba import jit
import pandas as pd 
from collections import Counter 



''' 
### record chapter lines
'''


def chpter (lns) :
    # 记录段落首行号，numbr of partcl是chpter行数
    print('chpter')
    chp_st, chp_len,ids = [], [],[]
    for idx, val in enumerate(lns) :
        if '# id type nb id_1...id_nb mol bo_1...bo_nb abo nlp q' in val :
            chp_st.append((idx+1))
        if '# Number of particles' in val :
            aa = val.strip()
            aa = aa.split(' ')
            chp_len.append(int(aa[-1]))
    
    idx = [{x,x+y} for x,y in zip(chp_st, chp_len)]
    
    for i, x in enumerate(idx) : 
        
        #print(f"{int(i)/int(tl)*100:.3f}")
        if i %skip == 0 :
            ids.append (x)
            

    return chp_st, chp_len, ids



''' # >>>>>>>>>> Main pro >>>>>>>>>>>>>>>>>>> '''

with open (file_names,'r') as f :   # 文件名
    lns = f.readlines()
print ('read completed')

chp_st, chp_len,idx = chpter (lns)
Q_lst = pd.DataFrame(Atom_type)
lst = pd.DataFrame()

print('main')
     
for val_cs, val_cl in idx : 
    print (val_cl,val_cs)

    Time = lns[val_cl-7].strip(' ').replace("\n", "").split(' ')                    # 时间行
    Time = int(Time[2])/timestep
    
    chp_q=[]                            # q sum list in each chpter
    for val in lns[val_cl : val_cs] :                    # enumerate each line to give the charge to the corspingding elemental llist
        data = val.strip(' ').replace("\n", "").split(' ')       # d[i].append(data[-1]), when at_type=2, adding at_type2's q in d[2] list
                                                                 # 去除列表中所有空格
        data = [x for x in data if x!='' ]
        chp_q.append([data[1], float(data[-1])])
    
    pd_q = pd.DataFrame(chp_q,columns=['at_type',str(Time)])
    pd_avq = pd_q.groupby('at_type').mean().T       # 以at_type分组,求平均，以T转职
    lst = pd.concat([lst,pd_avq])

    # ouput progressing ~
    print(f"progressing:  {val_cs}/{chp_st[-1]}    {int(val_cs)/int(chp_st[-1])*100:.3f}%")
lst.columns = Atom_type

lst.to_csv('bond_q.csv')






# pandas function: groupby, 
# grouped.sum()),list(grouped.size()



