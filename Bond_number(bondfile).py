
# 原子编号-关联原子列表，判断：列表中的原子是否为H,O，进而判断与原子形成的键。归纳到键列表NBOND中，去除重复的项。


'''
待修改一下时间的输入
'''
# bonding atom list generate

file_names =r'D:\纤维素裂解\2000-2800\result\\' 
outpath = file_names # r'//home/As135/Documents/Mo-ElectricField/10-7/result/'   
ele_name = ['C', 'H', 'O']
ele_n = len(ele_name)



#from numba import jit
import pandas as pd 
from collections import Counter 
import numpy as np
import os
import time
from tqdm import tqdm  

infs = [fs for fs in os.listdir(file_names) if 'bonds' in fs and 'out' in fs]
# infs = ['bonds.2000.out']
''' 
### record chapter lines
'''
def chpter (lns) :
    # 记录段落首行号，numbr of partcl是chpter行数
    chp_st, chp_len, st, lens = [], [], [], []
    for idx, val in enumerate(lns) :
        if '# id type nb id_1...id_nb mol bo_1...bo_nb abo nlp q' in val :
            chp_st.append((idx+1))
        if '# Number of particles' in val :
            aa = val.strip()
            aa = aa.split(' ')
            chp_len.append(int(aa[-1]))
        
    for i, x in enumerate(chp_len) :
    	if i %2 == 0 :
            lens.append(x)
    for j, y in enumerate(chp_st) :
    	if j %2 == 0 :
            st.append(y)
    # print (st, lens)
    # print (chp_st, chp_len)
    return st, lens



''' 
### record chapter lines
### 用于处理每一chpter，嵌套在product_list函数的循环里
'''
def bond_num (chapter) :
    chps = chapter
        # create 10 list
    d = [[]]*10         # 把每一行的原子组合都放进一个列表中，再去除重复的
                        
    # 原子类型列表
    data = pd.DataFrame(chps)
    data = data[0].str.split(' ', expand = True) 
    
    for i in range(10):  # 这里与上面一致
        d[i] = data[data[2].isin([str(i)])]
        d[i] = list(d[i][1])
    
    # 原子连接列表
    b_list = []
    for i_dx, val_i in enumerate(chps) :
        a = val_i.strip()
        a = a.split(' ')
    
        # 判断nb是否为0
        try : 
            if (len(a) > 2) and (a[2] != 0) :
                blst = a[3: (3+int(a[2]))] 
            
            # adding connecting atom (bond)
                for i_b, val_b in enumerate(blst) :
                    b_list.append(tuple(sorted([a[0],val_b])))  
            else :
                continue 
        except :
            continue
    
    # 去重
    b_list = list(set([x for x in b_list]))
    # 替换
    bb_lst = [list(x) for x in b_list ]
    
    for x, v_x in enumerate(b_list) :
        for xx, v_xx in enumerate(v_x) :
            
            for at in range(1,10) :
                #print (at)
                if v_xx in d[at] :
                    bb_lst[x][xx] = at
                    break
                else :
                    continue
    # (2,1) -> (1,2)        
    bb_lst = [tuple(sorted(x)) for x in bb_lst]
        # Counter 数数

    return bb_lst

    
def producelist(elenumber,chp_st, chp_len, outfname) :
    ele = elenumber+1
    # 6（elenumber）个元素 ，建立dataframe表的列名（1,1），（1,2）……
    bdnm = [tuple(sorted([nm,nn])) for nm in range(1,ele) for nn in range(1,ele)]
    bdnm =list(set([x for x in bdnm]))   # 去除重复的         
    
    # 建立0列表    
    pd_bond = pd.DataFrame(np.zeros((len(chp_st),len(bdnm)),dtype=int),columns=bdnm)
    
    # 主循环，传递每个章节到子函数
    # 用字典dic[name]赋值dataframe列表  pd.loc[i][name] = dict[name]  (列名与字典同名)
    # 
    
    n=0 
 
     
    for val_cs, val_cl in tqdm(zip(chp_st, chp_len),leave=True,position = 0) : 
        # print(f"progressing:  {val_cs}/{chp_st[-1]}    {(int(val_cs)/int(chp_st[-1]))*100:.3f}%")
        # can fun() return
        bb_lst = bond_num (lns[val_cs:val_cs+val_cl]) 
        bb = Counter(bb_lst)      # (1,1) : 2, (1,2):1……
        
        for v in bb :   # v 是name,bb是字典
            pd_bond.loc[n][v] = bb[v]   # good!  对第n行的每个单元赋值     
        n += 1
        
    listname = list(pd_bond.columns)
    # print(listname)
    
    # 更换列名为元素符号 ---------
    ele_col = {}
    ele_colst = []
    for i, e in zip(range(1,len(ele_name)+1), ele_name ) :
        ele_col[i] = e
    # print(ele_col)
    
    for i in listname :
        x = str(i)
        for ii in range(1,len(ele_name)+1) :
            x = x.replace(str(ii),str(ele_col[ii]))
        ele_colst.append(x)
            
    pd_bond.columns = ele_colst 
    # ---------
        
    pd.set_option('display.max_rows', 10000,'display.max_columns', 10000,"display.max_colwidth",10000,'display.width',10000)
    pd_bond.to_csv(outfname+'Bond_Number.csv')


    

''' # >>>>>>>>>> Main pro >>>>>>>>>>>>>>>>>>> '''

### 读取目录下的bond文件
for fs in tqdm(infs,leave=True,position = 0) :
### ---------------- ###

    with open (file_names+fs,'r') as f :   # 文件名
        lns = f.readlines()
    print ('read completed')
    
    outfname = outpath+fs.replace('.out','')
    # try :
        
    chp_st, chp_len = chpter (lns)
    producelist(ele_n, chp_st, chp_len,outfname)   # 元素个数
    
    # except :
    #     with open (outfname+'log.txt','a+') as flog :
    #         print (fs,file=flog)








