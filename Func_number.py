# -*- coding: utf-8 -*-
import os
from rdkit import Chem
from rdkit import RDConfig
from rdkit.Chem import Draw
from rdkit.Chem import FragmentCatalog
from rdkit.Chem import rdRGroupDecomposition as rdRGD
# from rdkit.Chem.Pharm2D.SigFactory import SigFactory
# from rdkit.Chem.Pharm2D import Generate, Gobbi_Pharm2D
from rdkit.ML import InfoTheory
from tqdm import tqdm 
import time
import pandas as pd 
from collections import Counter

fpth = r'/mnt/d/纤维素裂解/数据处理/产物/2600.spec.tab_sps-2.csv'
m = 1   # 出现次数少于m次的产物
# # !!!!!!!!! 下一句，速度慢可以分开读

# -----------------------------------------------------------------
# 获取官能团库,     # 分析每个分子的官能团列表，集合到一起，分类汇总
fName = os.path.join(
    RDConfig.RDDataDir,
    'FunctionalGroups.txt' )

# 根据官能团库实例化一个参数器
fparams = FragmentCatalog.FragCatParams(1, 6, fName)
# 查看官能团库中包含的官能团数量
fparams_num = fparams.GetNumFuncGroups()

mols = []
# 查看每个官能团对应的集团
for i in range(fparams_num):
    mols.append(fparams.GetFuncGroup(i))

# 可视化官能团库
img = Draw.MolsToGridImage(mols, molsPerRow=8)
img.save( 'func.jpg' )

###########################   Modify   #################################

import func_timeout
@func_timeout.func_set_timeout(0.5)
def funs(m, fcgen):
    # 计算分子片段
    fcgen.AddFragsFromMol(m, fcat)
    # 查看分子片段数量
    num_entries = fcat.GetNumEntries()
    # print ('fun_end')
    return num_entries, fcgen

func = []         # 官能团种类
sps_lst_New = []  # 按顺序记录有官能团的物种列表
func_num = {}     # 每种产物的官能团字典

pf = pd.read_csv(fpth, header = 1,index_col = 0) 
print ('有',len(pf.columns))

# 删除全为0的列，减少smiles式的检验
non0Lst = []
for col in pf :
    if (pf[col].eq(0).all()) or (pf[col].sum() < m) : 
        continue
    else :
        non0Lst.append(col)

print ('出现m次以上的产物有',len(non0Lst))

# !!!!!!!!!
pf = pf[non0Lst]    # non0Lst[:1000]   '''   '''

sps_lst = list(pf.columns)[1:]

# 创建一个片段生成器，通过该对象生成片段
fcat = FragmentCatalog.FragCatalog(fparams)

for sps in tqdm(sps_lst) :
    try :
        m = Chem.MolFromSmiles(sps)
        fcgen = FragmentCatalog.FragCatGenerator()
        try : # 超时中断处理
            num_entries, fcgen = funs(m, fcgen)
            # 通过存储器查看片段
            # for i in range(num_entries) :
            #     print(fcat.GetEntryDescription(i))
            # 向存储器传入分子片段id ， 获取片段中所包含的官能团标号 ： GetEntryFuncGroupIds
        except func_timeout.exceptions.FunctionTimedOut:  # 超时异常
            print(sps,'异常')
            continue 
    except :
        continue

    if num_entries > 1 :
        sps_lst_New.append(sps)
        entries_group_ids = fcat.GetEntryFuncGroupIds(num_entries-1)
        # print('matched the function group ids is', list(entries_group_ids))
        # 物种官能团数量,创建字典
        func.extend(list(entries_group_ids))
        cnt = dict(Counter(list(entries_group_ids)))
        func_num[sps] = cnt
        
        # 向参数器传入官能团编号，获取官能团对应的mol对象
        # for ii in entries_group_ids :
        #     fg = fparams.GetFuncGroup(ii)
        #     print('name of group', fg.GetProp('_Name'))  # name of group 1 -C(=O)O

            # img = Draw.MolsToGridImage(mols, molsPerRow=2)
            # img.save(
            #     './mol'+str(ii)+'.jpg'
            #     )


# 官能团信息            
print(len(sps_lst),'个产物中有',len(sps_lst_New),'个包含官能团')
  
# 根据产物是数量列表计算时刻官能团的数量
func = list(set(func))
print ('官能团有：',func)

# 查看每个官能团对应的集团
fun_mol = []
fun_moll = {} # 官能团结构式列表
for i in func :
    fun_mol.append(fparams.GetFuncGroup(i))

# 官能团结构
for i in func :
    fg = fparams.GetFuncGroup(i)
    stru = fg.GetProp('_Name')
    print(i,' of group', stru)  # name of group 1 -C(=O)O
    fun_moll[i]=stru
    
img = Draw.MolsToGridImage(fun_mol, molsPerRow=5)
img.save('./mol'+'.jpg')
 
    
# 官能团数据表
pfunc = pd.DataFrame(None, columns = func) 
# df = df.groupby(df.columns, axis=1).sum() 
for i_sps in sps_lst_New :  # 产物遍历

    for i_func in func :    # 官能团遍历

        try :
            x = func_num[i_sps][i_func]   # 该物种是否有这个官能团，有就计数
            if pfunc[i_func].isnull().any() :  # 判断是否是空列
                pfunc [i_func] = pf [i_sps]*x                   # 官能团数值
            else:
                pfunc [i_func] = pfunc [i_func] + pf [i_sps]*x  # 官能团数值
        except :
            pass
# 更换列名

for i in func :
    pfunc = pfunc.rename(columns = {i : fun_moll[i]})
pfunc.insert(0,'Time', pf.iloc[:,0])

print (pfunc)

pfunc.to_csv('2600_function_number-2-3400.csv')

