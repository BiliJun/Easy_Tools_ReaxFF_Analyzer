# -*- coding: utf-8 -*-
"""
识别段落
始末语句
Per MPI rank memory allocation (min/avg/max)
Loop time of

"""
import pandas as pd 
import os
import re
from cycler import cycler
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
# Weight_Radius
# Weight_Radius
# time

# 输入

fname = '2200'                                 # 文件名中的字符筛选
fth = r'D:\zh\JP-10\结果\运动位移\转动惯量\\'    # 输入路径
op = r'D:\zh\JP-10\结果\运动位移\\'                # 输出路径
txtlabel = '2200 K'                               # 标注
# 图例标签
# lglabel = ['non-E','$10^{-5}$ V/Å','0.5 V/Å']
# lglabel = ['2100K','2200K','2300K','2400K']
lglabel =['non-E','$10^{-3}$ V/Å','$10^{-2}$ V/Å','0.1 V/Å','0.5 V/Å','$10^{-5}$ V/Å', '$10^{-6}$ V/Å']
#   '$10^{-6}$ V/Å'

# 颜色
colorlist = ['black','red','mediumpurple','lightcoral','deeppink','red','deepskyblue','blue']
#['black','deeppink','blue']
# ['black','deeppink','red','lightcoral','mediumpurple','deepskyblue','blue'] #'turquoise','deepskyblue','slategrey','cornflowerblue','blue'             ]
mpl.rcParams['axes.prop_cycle'] = cycler(color=colorlist)



pth = [p for p in os.listdir(fth) if fname in p and 'csv' in p ]
plt.figure(figsize=(5, 4))

    
for p in pth :
    print(p)
    pf = pd.read_csv(fth+p,sep = ',')
    from scipy.signal import savgol_filter
    y = savgol_filter(pf['Weight_Radius'],100, 3, mode= 'nearest')

    plt.plot( pf['time']/2,y,alpha = 1)
    # plt.plot( pf['time']/2,pf['Weight_Radius'],alpha = 0.1)    
    plt.title(txtlabel,x=0.1,y=0.85,color='black',fontsize=18,fontproperties='Times New Roman') #,backgroundcolor='w')    
    plt.xlabel('Time (ps)',fontdict={'family':'Times New Roman', 'size' : 18})
    plt.ylabel('Weight gyration radius (Å)',fontsize=18,fontproperties='Times New Roman',backgroundcolor='w')
    plt.xticks(fontproperties='Times New Roman',fontsize=12)
    plt.yticks(fontproperties='Times New Roman',fontsize=12)
    # plt.title(sps,x=0.9,y=0.85,fontsize=20,fontproperties='Times New Roman')
    plt.grid(True,linestyle = '--', linewidth = 0.5)
    # patches = [ mpatches.Patch(color=colorlist[i], label="{:s}".format(lglabel[i]) ) for i in range(3) ]
    # plt.legend(handles=patches, loc=1, ncol=6) 
    # plt.legend( lglabel, loc=3, ncol=1, bbox_to_anchor=(0.65,0.7) )
    fig = plt.gcf()
   
for p in pth :
    pass
    pf = pd.read_csv(fth+p,sep = ',')
    # plt.plot( pf['0'],pf['0.1'],'-',alpha = 1)
    print(p)
    y = savgol_filter(pf['Weight_Radius'],100, 3, mode= 'nearest')
    plt.plot( pf['time']/2,pf['Weight_Radius'],alpha = 0.1)    
    # plt.plot( pf['time']/2,y,alpha = 1)
    plt.title(txtlabel,x=0.1,y=0.85,color='black',fontsize=18,fontproperties='Times New Roman') #,backgroundcolor='w')    
    plt.xlabel('Time (ps)',fontdict={'family':'Times New Roman', 'size' : 18})
    plt.ylabel('Weight gyration radius (Å)',fontsize=18,fontproperties='Times New Roman',backgroundcolor='w')
    plt.xticks(fontproperties='Times New Roman',fontsize=12)
    plt.yticks(fontproperties='Times New Roman',fontsize=12)
    # plt.title(sps,x=0.9,y=0.85,fontsize=20,fontproperties='Times New Roman')
    plt.grid(True,linestyle = '--', linewidth = 0.5)
    # patches = [ mpatches.Patch(color=colorlist[i], label="{:s}".format(lglabel[i]) ) for i in range(3) ]
    # plt.legend(handles=patches, loc=1, ncol=6) 
    plt.legend( lglabel, loc=3, ncol=3, bbox_to_anchor=(0.05,1)) 
    fig = plt.gcf()
    
plt.show()
fig.savefig(op+fname+'.tif', bbox_inches='tight',dpi=900)  