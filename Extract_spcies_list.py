# -*- coding: utf-8 -*-
"""
products： 输入要提取的产物名称
"""
# import pandas as pd
import os
import time
import pandas as pd

folder = r'D:\ActiveEnergy\产物\\'
fils = os.listdir(folder)
fils = [fff for fff in fils if ('csv' in fff)]

products ='CH4'
products = products.split()
    
st = time.time()

for ns in products :
    pass
    try :
    # lst=pd.DataFrame()
        lst={}
        for i_file in fils :
            pass
            data = pd.read_csv(folder+i_file,engine='c',low_memory=False)  
            lst['t'+ns+i_file]=data['Time']
            lst[ns+i_file]=data[ns]
        Lst = pd.DataFrame(lst)
        Lst.to_csv(folder+ns+'.csv')
    except:
        print(ns+'   No')
et = time.time()
print (et-st)

