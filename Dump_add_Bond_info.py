# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 15:05:17 2023

@author: Administrator
"""
import pandas as pd

# 电荷云图
STEP = [s for s in range(0,1000,200)]

def yuntu(data1,f ) :
    import matplotlib.pyplot as plt 
    import numpy as np 
    
    # 原始坐标 x,y 和对应的状态变量值 a
    from scipy import interpolate
    T = data1
    #盒子边界
    T['xs'] = T['xs'].astype(float)*(f[1]-f[0])+f[0]
    T['ys'] = T['ys'].astype(float)*(f[3]-f[2])+f[2]
    T['zs'] = T['zs'].astype(float)*(f[5]-f[4])+f[4]
    # 筛选纵坐标
    T1 = T[T['zs'] < 10]
    T1 = T1[T1['type'] == 3]

    
    h = 1000
    X = np.linspace(min(T1['xs']),max(T1['xs']), h)
    Y = np.linspace(min(T1['ys']),max(T1['ys']), h)
    #生成二维数据坐标点，可以想象成围棋棋盘上的一个个落子点
    X1, Y1 = np.meshgrid(X, Y)
    
    Z = interpolate.griddata((T1['xs'], T1['ys']), T1['q'], (X1, Y1), method='cubic')
    
    plt.figure()
    plt.subplot(111)
    extent=(min(T1['xs']), max(T1['xs']),min(T1['ys']), max(T1['ys'])) #任意设置显示的坐标范围
    plt.imshow(Z,cmap='RdBu_r',extent=extent, vmin=min(T1['q']),vmax=max(T1['q']))
    # plt.title('$\Omega$') #'plt.cm.rainbow
    plt.colorbar()
    plt.show()
    # ,vmin=min(T['q']),vmax=max(T['q'])
    
    
import os
for p in list(os.listdir(r'D:\ECH4_MoSNi\Data_CH4Elec\PartE_Ni-2\\')):

    pth = r'D:\ECH4_MoSNi\Data_CH4Elec\PartE_Ni-2\\'+p+r'\result\\'
    foupt = open (pth+'dump.new.lammps','w')
    print(f'processing {p}')
    with open (pth+'dump.reax.lammpstraj') as ff1, open(pth+'bonds.out')as ff2 :
        f1 = ff1.readlines()
        print('Read dump file completely')
        f2 = ff2.readlines()
        print('Read bond file completely')
        box = [ float(ii) for i in f1[5:8] for ii in i.strip().split()  ] # 盒子边界
        row_No1 = [i+1 for i, t in enumerate(f1) if 'TIMESTEP' in t]
        row_No2 = [i for i, t in enumerate(f2) if '# Timestep' in t]
        row_No1 = row_No1[abs(len(row_No1)-len(row_No2)):]  # bond与lammps文件的步数差值
        from tqdm import tqdm
        for i in tqdm(range(len(row_No1)-1)) :
            pass
            dat1 = pd.DataFrame ([ l.strip().split() for l in f1[row_No1[i]+8:row_No1[i+1]-1]]) # dump
            dat1.columns = ['ITEM: ATOMS id', 'type', 'xs', 'ys', 'zs']
            dat1 = dat1.astype(float)
            dat2 = pd.DataFrame ([[ l.strip().split()[0], l.strip().split()[1], l.strip().split()[-1] ] for l in f2[row_No2[i]+7:row_No2[i+1]-1]]) # bond
            dat2.columns = ['id', 'type','q']
            dat2 = dat2.astype(float)
            dat1 = dat1.sort_values(by='ITEM: ATOMS id', inplace=False)
            dat2 = dat2.sort_values(by='id', inplace=False)
            dat1['q'] = dat2['q']
            #设置value的显示长度为200，默认为50
            pd.set_option('max_colwidth',200)
            #显示所有列，把行显示设置成最大
            pd.set_option('display.max_columns', None)
            #显示所有行，把列显示设置成最大
            pd.set_option('display.max_rows', None)
            pd.set_option('display.unicode.east_asian_width', True)
            dat1['ITEM: ATOMS id'] = dat1['ITEM: ATOMS id'].astype(int)
            dat1['type'] = dat1['type'].astype(int)
            
            if i in STEP :
                yuntu(dat1,box) 
                
                
            print(f'ITEM: TIMESTEP\n{i+1}',file = foupt)
            for ff in f1[2:9] :
                print(ff,end='',file = foupt)
            print(dat1.to_string(header =False,index=False),file = foupt)
    
    foupt.close()
            