# -*- coding: utf-8 -*-
"""


"""
# Input file folder
# fils =  [str(x) for x in range(1,9)]
# fils.extend(['62','64','1+'])

fils =  ['64']









''' >>>>>>>> Main program <<<<<<<<<< '''

import pandas as pd 
import collections
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import numpy as np 




def ring_cal (itx) :
    
    chsmiles, chformule, ty_sm_n, TolRing   = [], [], [], []
    Lo, Hi = [],[]
    for ii, x in enumerate(itx) :
        # pass
        # smiles line
        chsmiles.append(x)
        # smiles -> chemical formular
        mol = Chem.MolFromSmiles(x)
        formula = CalcMolFormula(mol)
        chformule.append(str(formula))
  
        # total rings number
        r_num = mol.GetRingInfo()
        tl_r = r_num.NumRings()
        TolRing.append(tl_r)
        
        # 统计最小环的个数
        ssr = Chem.GetSymmSSSR(mol)
        # 几个环
        n_type = len(ssr) 
        gpb = []
        for ring in ssr:
                                        # 列出ssr的所有环及其每个环的原子组成
            # 判断环的类型，数原子个数
            r_type = len(list(ring))    #  ssr : 5-ring [1,2,3,4,5]
            gpb.append(r_type)          #        6-ring [1,2,3,4,5,6]
        
        # 确定环的最大、最小类型
        try:
            Lo.append(min(gpb) )
            Hi.append(max(gpb) )
        except :
            a=1

        
        # ring type for each sps
        ty_sm_n.append(gpb)
        
    lo = min(Lo)
    hi = max(Hi) 
    R_np = []
    for i, x in enumerate(ty_sm_n) :

        if x == [] :

            r_np = np.zeros((hi-lo+1), dtype = int) 
            R_np.append(r_np)
        
        else :
            molecu = collections.Counter(x) 
            # molecu_k = molecu.most_common(len(molecu))
            
            r_np = np.zeros((hi-lo+1), dtype = int)    # molecu
            
            for i in molecu :
                r_np[i-lo] = molecu[i]
            
            R_np.append(r_np) 
        
    R_np = np.mat(R_np)
    R_names = [str(na)+'-Ring' for na in range(lo,hi)]
   
    return chsmiles, chformule, R_np, R_names, TolRing

                        

def speslst (tar_file):
    # 读取产物列表
    # tar_file = tar_folder+oun+'.spec.tab'
    # skip chem_SMILES lines in the large file
    # head charcter is ';'
    # smiles to fomular moudle
    
    with open (tar_file,'r') as ff :

        line = ff.readlines ()
        print ('read complete')
        
        for i, im in enumerate (line) : # im = line[1] 
            if i == 0 :
                continue            
            if im[0] == ';' :
                itx = im.strip(' ')
                itx = itx.split(';')
                
                chsmiles, chformule, R_np, R_names, TolRing  = ring_cal(itx)
                print ('ring calculation')    
                    
            if im[0] != ';' :
                skip_lns = i
                break
            
    #print (chsmiles)
    t = [x.strip('\n') for x in line[skip_lns+1:]]
    Sn = line[0].strip('\n').split(';')
    Chm = pd.DataFrame(chsmiles).T
    Chm.to_csv(tar_folder+out+'smiles.csv')
    
    
    sp = pd.DataFrame(t) 
    sp = sp[0].str.split(';', expand = True)

    sp.set_index([0], inplace=True)
    sp.columns = chformule[1:]


    ''' 环列表 '''
    print('Building Ring type list...')
    lst=[]
    
    for i, inx in enumerate(R_names)  :
        c= []                           # 1 ring num list
        fa = [fx for ff in R_np[1:,i].tolist() for fx in ff]
        for xi in range(sp.shape[0]) :          # time
            a = [] # 1 time sum
            for ii, ifa in zip(sp.iloc[xi],fa) : # ele
                a.append(int(ii)*int(ifa))        
            c.append(sum(a))                    # 1 time sum
        print(f'{inx} list')                                        # 1 ring num list
        lst.append(c)
        
    lst=pd.DataFrame(lst)
    m=lst.T
    m.columns=R_names
    
    m.index = list(sp.index)
    
    m.to_csv('./Ring_number/'+out+'Ring_num.csv')
    print(f'Ring list of .csv file complete. Ring type : {R_names}')
    of = open (oun+'name.txt','w+') 
    print(chsmiles, file=of)
                                                               
    return sp    

                                              

def order_lst (lst) :
# 排序，提取前10产物
    decrib = sp.describe()
    df = pd.concat([sp,decrib],axis=0)
    lst_order = df.sort_values(by=['unique'], axis=1, ascending=False)
    
    return lst_order



for fi in fils :
    
    oun = fi
    out = 'rp'+oun
    tar_folder = './rp'+oun+'/'

    print(f'Folder : {tar_folder}')
    sp = speslst (tar_folder+oun+'.spec.tab')
    
    # sort by mean
    lst_order = order_lst (sp)
    print ('sps list construction')
    # 输出

    # CSV output
    # lst_order.to_csv (tar_folder+out+'bd_sps.csv' )





