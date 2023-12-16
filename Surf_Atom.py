# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 12:33:54 2023
2023/11/7 修正：xs,ys,zs坐标    x = xs*(xhi-xlo)+xlo
@author: Administrator
"""

# -*- coding: utf-8 -*-
"""
表面原子数量
计算反应物与催化剂间距离小于d的原子，输入它们的原子类型标号

表面原子密度=表面原子数/催化剂的上层原子数
"""
 
from collections import Counter
from numba import jit

@jit   # 
# @nb.jit(nopython=True,parallel=True)
def distance (idCHO,tCHO,xCHO,yCHO,zCHO,  idMoS,tMoS,xMoS,yMoS,zMoS) :
    suf = []
    Atomid = []
    for ind in range(len(tCHO)):
        # 反应物原子id 2367原子与MoS中某个原子距离小于cutoff
        # 则不与其他原子比较了，此时记录下id 2367到到Atomid[]，避免下次重复比较
        # 在循环第一步检查id 2367原子是否被录入到Atomid[]中
        # 计算每个反应物原子和催化剂原子之间的距离
        # 添加符合的反应物原子的类型标号到suf[]，方便归类每种类型原子的数量
        if idCHO[ind] not in Atomid :
            x1 = xCHO[ind]
            z1 = zCHO[ind]
            y1 = yCHO[ind]
            
            for jnd in range(len(tMoS)) :
                x2 = xMoS[jnd]
                y2 = yMoS[jnd]
                z2 = zMoS[jnd]
                dist = ((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)  
                
                if dist < target_d**2 :
                    suf.append(tCHO[ind])   # 添加符合的原子的类型标号到suf[]，方便归类每种类型原子的数量
                    # Atomid.append(idCHO[ind])
                    # break          # !!!!!!!! 如果计算全部的连接，注释掉break                    
                                                       
        else :
            print ('distance error atom',idCHO[ind])
            continue
                    
    return suf

@jit  
def distance0 (idCHO,tCHO,xCHO,yCHO,zCHO,  idMoS,tMoS,xMoS,yMoS,zMoS) :
    suf = []
    Atomid = []
    for ind in range(len(tCHO)):
        x1,y1,z1 = xCHO[ind], yCHO[ind],zCHO[ind]      
        for jnd in range(len(tMoS)) :
            x2,y2,z2 = xMoS[jnd], yMoS[jnd],zMoS[jnd]
            dist = ((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)**0.5     
            if abs(dist) < target_d :
                suf.append(tCHO[ind])   # 添加符合的原子的类型标号到suf[]，方便归类每种类型原子的数量
                Atomid.append(idCHO[ind])
    return suf

def ploting(dataplot) :
    import matplotlib.pyplot as plt
    x = dataplot['5']   
    for i in dataplot.columns :
        if i != 'time' :
            y = dataplot[i]
     
        plt.plot(x,y)
    plt.show()

# >>> data paragraph constitute dataframe
def main() :
    import pandas as pd
    
    N_suf = []

    times = []
    try:

        for i in range(0,len(Para)-1,skip) :  # 跳步
            pass
            print ("Total  : ", len(Para), 'now is : ', i, i/len(Para)*100,"%")
            print ("TIMESTEP : ", txt[Para[i]+1])
            
            # 读取bond文件的每个章节
            dat = [] # 章节数据列表
            suf = [] # 催化剂表面的原子列表
            dat = pd.DataFrame (txt[Para[i]+9:Para[i+1]])   # [1:3] is 1,2   representing [1:3) 
            dat = dat[0].str.split(' |\n',expand=True)
            dat.columns=['id', 'type', 'xs', 'ys', 'zs', 'null']
            
            ''' 
            # >>>>>>>  define the atom type >>>>>>>>>>>>>>>>>>>>>
        '''
            # pandas筛选选定的原子类型的数据行
            CHO = dat[dat['type'].str.contains(reactant)]
            MoS = dat[dat['type'].str.contains(catalyst)]
            # print ('iteration ... ')

            # 所有选定的反应物的原子编号, 原子类型标号， 原子坐标

            if ifscale == 'yes' :
                xCHO = CHO['xs'].astype(float)*(scale[0][1] - scale[0][0]) + scale[0][0] 
                yCHO = CHO['ys'].astype(float)*(scale[1][1] - scale[1][0]) + scale[1][0]
                zCHO = CHO['zs'].astype(float)*(scale[2][1] - scale[2][0]) + scale[2][0]
                xMoS = MoS['xs'].astype(float)*(scale[0][1] - scale[0][0]) + scale[0][0]
                yMoS = MoS['ys'].astype(float)*(scale[1][1] - scale[1][0]) + scale[1][0]
                zMoS = MoS['zs'].astype(float)*(scale[2][1] - scale[2][0]) + scale[2][0]
                xCHO, yCHO, zCHO = xCHO.tolist(), yCHO.tolist(), zCHO.tolist()
                xMoS, yMoS, zMoS = xMoS.tolist(), yMoS.tolist(), zMoS.tolist()
                idCHO, tCHO, idMoS, tMoS = CHO['id'].tolist(), CHO['type'].tolist(), MoS['id'].tolist(), MoS['type'].tolist()

            else :
                idCHO, tCHO, xCHO, yCHO, zCHO = CHO['id'].tolist(), CHO['type'].tolist(), CHO['xs'].tolist(), CHO['ys'].tolist(), CHO['zs'].tolist()
                idMoS, tMoS, xMoS, yMoS, zMoS = MoS['id'].tolist(), MoS['type'].tolist(), MoS['xs'].tolist(), MoS['ys'].tolist(), MoS['zs'].tolist()


            # suf = []
            # for ind in range(len(tCHO)):
            #     x1,y1,z1 = float(xCHO[ind]), float(yCHO[ind]),float(zCHO[ind])
                
            #     for jnd in range(len(tMoS)) :
            #         x2,y2,z2 = float(xMoS[jnd]), float(yMoS[jnd]),float(zMoS[jnd])
            #         dist = ((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)**0.5  
                    
            #         if abs(dist) < 0.2 :
            #             suf.append(tCHO[ind])
            ''' 计算distance '''
            suf = distance (idCHO,tCHO, xCHO, yCHO, zCHO, idMoS,tMoS, xMoS, yMoS, zMoS)
            # atid.extend(atomid)
            
            # 统计催化剂表面每种反应物原子的数量
            n_suf = Counter(suf)    # 一个章节的数据
            N_suf.append (n_suf)    # 所有章节的原子类型列表
            times.append(i)
            
        dfsur = pd.DataFrame(N_suf) 

        #########################################
        # 更换列名为元素符号
        listname = list(dfsur.columns)

        for i in listname :
            x = str(i)
            for ii in range(1,len(Atom_type)+1) :
                x = x.replace(str(ii),str(ele_col[ii]))
            ele_colst.append(x)
                
        dfsur.columns = ele_colst 
        ###################################
        dfsur.insert(0,'time',times)
        print(list(dfsur.columns))
        
        dfsur.to_csv(ouptpath+fo+'.csv') 
        
        
        import matplotlib.pyplot as plt
        x = dfsur['5']   
        for i in dfsur.columns :
            if i != 'time' :
                y = dfsur[i]
         
            plt.plot(x,y)
        plt.show()
        
        ploting(dfsur)
        
    except :
        print ('file error',i)





# ----------
# --------------------------------
# ----------------------------------------------------------------------------------------
# setting 
reactant = '4|5'              # reactant：反应物原子类型标号；
catalyst = '1'            # catalyst：反应物原子类型标号;catalyst = '1|3' 是包含1和3号原子催化剂
Atom_type = ['Mo', 'Ni', 'S', 'C','H']       # 元素符号，按孙旭列出
target_d = 2.5              # 原子间距离 A埃
skip = 2



# ----------
# --------------------------------
# ----------------------------------------------------------------------------------------
# 文件路径，读入的文件名，
import os
p = r'D:\ECH4_MoSNi\Data_CH4Elec\PartE_Ni-2\nonE\result\\'
pp = [ppp for ppp in os.listdir(p) if ('dump' in ppp) and ('.sh' not in ppp)]

ouptpath = r'D:\ECH4_MoSNi\数据处理\表面原子数\Ni\\'

ifscale = 'yes'

for i_pp in pp :
    pass
    # fpaths = r'D:\ECH4_MoSNi\Data_CH4Elec\PartE_Ni-2\nonE\result'+str(i_pp) + r'\result\\'
    fpaths = r'D:\ECH4_MoSNi\Data_CH4Elec\PartE_Ni-2\nonE\result\\'
    fnm = r'dump.reax.lammpstraj'
    
    
    # 元素编号字典  ---------
    ele_col = {}
    ele_colst = []
    for i, e in zip(range(1,len(Atom_type)+1), Atom_type ) :
        ele_col[i] = e
    print(ele_col)
    
    
    ''' 输出文件名 ''' # 催化剂元素-反应物元素
    fo = str(i_pp)+"_"+ ("").join([ele_col[int(i)] for i in catalyst.split('|')]) +"-"+ ("").join([ele_col[int(i)] for i in reactant.split('|')])







    
# -------------------------------------------------------------    
    
    # >>> readlines for large dataset                   # define function can't transfer this so big dataset
    
    with open(fpaths+fnm,'r') as f :        # so we have to write in normal var
        txt = f.readlines()
        txt_len = len(txt)
    f.close()
    scale = [ list(map(float,i.strip().split())) for i in txt[5:8] ]
    # >>> Para[] "ITEM: TIMESTEP"+7 is the data area
    Para = []                       # 章节的起始行号列表
    for i in range(txt_len) :
        a = txt[i].strip()
        if "TIMESTEP" in a :
            Para.append(i)  
    Para.append(txt_len+1)
    print ('read over')
            
    main ()


    





# def Distance (a,b) :
#     x1,y1,z1 = float(a[0]), float(a[1]),float(a[2])
#     x2,y2,z2 = float(b[0]), float(b[1]),float(b[2])
#     dist = ((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)**0.5                          
#     return dist

        
    # for ind, xyz1 in CHO.iterrows():

    #     a = list(xyz1[2:5])
    #     print ('iteration ',len(dat), 'now is : ', ind )
        
    #     for jnd, xyz2 in MoS.iterrows() :
    #         b = list(xyz2[2:5])
    #         d = Distance (a,b)
            
    #         if d < 0.3 :
    #             suf.append(xyz1[1])
                
    # n_suf = Counter()
    # print(n_suf)
    # N_suf.append (n_suf)
        

    


                
            








