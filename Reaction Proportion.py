# -*- coding: utf-8 -*-


import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx
import os
file_name = r'D:\zh\JP-10\REAXFF\JP-10-2200\bond\bondanalysis\\'
pth= [p for p in os.listdir(file_name) if ('2200' in p) and ('csv'not in p) ]
Summ = pd.DataFrame()
for p in pth :
#     pass
    file_names =file_name+p+'\Default.rate.tab'
    df=pd.read_csv(file_names,sep=';',engine='c',low_memory=False)
    N = df[df.index%2==0]['N'].values -df [df.index%2==1]['N'].values
    k = df[df.index%2==0]['k'].values -df [df.index%2==1]['k'].values
    
    df.loc[df.index%2==0,'N_ave'] = N
    df.loc[df.index%2==1,'N_ave'] = -N
    df.loc[df.index%2==0,'k_ave'] = abs(k)
    df.loc[df.index%2==1,'k_ave'] = abs(k)
      
    
    
    df[['R','P']] = df["Formula's"].str.split('->',expand=True)       # 分列
    df[['RS','PS']] = df["SMILES"].str.split(':',expand=True)   
    # 奇数行，正反应
    df=df[df.index%2==0]
    df = df.reset_index(drop=True)
    
    for index, row in df.iterrows() :
        I = row['R'].replace(' ','').split('+')
        J = row['P'].replace(' ','').split('+')
        df.loc[index,'R1'] = len(I)
        df.loc[index,'R2'] = len(J)
        
        if (len(I) == 1) and (len(J) == 1) :
            df.loc[index,'Rtype'] = 'iso'
        if (len(I) == 1) and (len(J) > 1) :
            df.loc[index,'Rtype'] = 'crack2'
        if len(I) > 1 :
            df.loc[index,'Rtype'] = 'add'
        if (len(I) == 1) and (len(J) == 2) :
            df.loc[index,'Rtype'] = 'crack1'
            if 'H' in J :
                df.loc[index,'crackType'] = 'H-abtract'
            else :
                df.loc[index,'crackType'] = 'C-C_crack'
                
            
    
    '''
    寻找反应列表
    '''
    a = df[df['RS'].str.contains(r"C1C[C@H]2[C@H](C1)[C@H]1C[C@H]2CC1")]
    a = df[df['R'].str.contains("C10H16")]
    # a.to_csv(r'D:\zh\n-alkane\结果\\C7_'+p+'_reactList.csv')
    # b = df[df['R'].str.contains('C16H34')]
    # c = df[df['R'].str.contains('C7H14')]
    
    dftol1 = df.groupby(by=['Rtype'])['N_ave'].sum()
    dftol2 = df.groupby(by=['crackType'])['N_ave'].sum()
    dfa = a.groupby(by=['Rtype'])['N_ave'].sum()
    # dfb = b.groupby(by=['Rtype'])['N'].sum()
    # dfc = c.groupby(by=['Rtype'])['N'].sum()
    dfaH = a.groupby(by=['crackType'])['N_ave'].sum()
    # dfbH = b.groupby(by=['crackType'])['N'].sum()
    # dfcH = c.groupby(by=['crackType'])['N'].sum()
    # print(dfa)
    # print(dfaH)
    Summ[p] = pd.concat([dfa,dfaH])
Summ.to_csv(r'D:\zh\JP-10\REAXFF\JP-10-2200\bond\bondanalysis\\'+p+'TolreacType.csv')
print(Summ) 
# f = pd.concat([dfa, dfb, dfc],keys=['C12H26', 'C16H34', 'C7H14'], axis=1)
# f["总和"] =f.apply(lambda x:x.sum(),axis =1)
# f.loc[len(f),:] =f.apply(lambda x:x.sum(),axis =0)
# ff = pd.concat([dftol1,dftol2])
# print(f)