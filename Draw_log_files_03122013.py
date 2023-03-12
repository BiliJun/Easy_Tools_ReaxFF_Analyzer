# -*- coding: utf-8 -*-
"""
"""
import argparse
import pandas as pd 
import os
import re
from cycler import cycler
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import sys
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
import lmpstrVal
parser = argparse.ArgumentParser(description='Test for argparse')
parser.add_argument('--transparency', '-t', help='Curve transparency,default = 0',default=10)
parser.add_argument('--filefolder', '-i', help='input file path,defaultpath is current folder',default=os.getcwd()+r'\\')
parser.add_argument('--outfolder', '-o', help= '''input file path,the path ending charcter need tobe "\" or "/" in win or linux system 
                    otherwise create in upper folder, defaultpath is current folder''',default=os.getcwd())
parser.add_argument('--figlab', '-l', help='figure label',default='')
args = parser.parse_args()

colorlist = ['black','deeppink','red','lightcoral','mediumpurple','deepskyblue','blue'] #'turquoise','deepskyblue','slategrey','cornflowerblue','blue'             ]
mpl.rcParams['axes.prop_cycle'] = cycler(color=colorlist)
fth = args.filefolder

import platform
 
if platform.system().lower() == 'windows':
    print("windows")
    os.mkdir(args.outfolder+r'/results/')
    outfolder = args.outfolder+r'/results/'
elif platform.system().lower() == 'linux':
    os.mkdir(args.outfolder+r'\\results\\')
    outfolder = args.outfolder+r'\\results\\'
    print("linux")
    
pth = [p for p in os.listdir(fth) if ('csv' not in p) and ('tif' not in p) and ('.py' not in p)]
txtlabel = args.figlab

def fun_para (path) :   
    with open (path,'r') as fun_f :
        fun_ff = fun_f.readlines()
        f_sart_l = [fun_i for fun_i , fun_fl in enumerate(fun_ff) if 'Per MPI rank' in fun_fl ] 
        f_end_l = [fun_i for fun_i , fun_fl in enumerate(fun_ff) if 'Loop time of' in fun_fl ] 
        f_time_l = [fun_i for fun_i , fun_fl in enumerate(fun_ff) if 'Time step' in fun_fl ] 
    if f_end_l[-1] < f_sart_l[-1] :
        f_endl = len(fun_ff)-1
    else :
        f_endl = f_end_l[-1]-1 
    
    fun_pf = pd.DataFrame([re.split(r"[ ]+",fun_fl.strip()) for fun_fl in fun_ff[f_sart_l[-1]+2:f_endl]]).astype(float)
    f_col = re.split(r"[ ]+",fun_ff[f_sart_l[-1]+1].strip())

    lammps_value = lmpstrVal.lmpstrs() 
    fcol = [lammps_value[f_var] if f_var in lammps_value else f_var for f_var in f_col]
    fun_pf.columns = fcol

    if fun_pf.loc[0]['Step'] > 10000 :
        fun_pf['Step'] = fun_pf['Step'] - fun_pf.loc[0]['Step'] 
    else :
        fun_pf['Step'] = fun_pf['Step']*f_time_l[-1]/1000
    
    return fun_pf, fcol[1:]
        
def fun_pands (xcol,pth) :
    pff = pd.DataFrame()
    for p in pth :
        pf, temp = fun_para (fth+p)
        pff[p] = pf[xcol]
    return pff

temp, Xcol = fun_para (fth+pth[0])
Xlabel = Xcol

for xcol,xlabel in zip(Xcol,Xlabel) :
    try :
        xlist = fun_pands (xcol,pth)
        xlist.to_csv(outfolder+xcol[:5]+'.csv')
        from scipy.signal import savgol_filter
        for xl in xlist : 
            y= savgol_filter(xlist[xl], args.transparency, 3, mode= 'nearest')
            plt.plot(y,label=xl)
            print(xl)
            
        plt.title(txtlabel,x=0.85,y=0.85,color='black',fontsize=18,fontproperties='Times New Roman') #,backgroundcolor='w')    
        plt.xlabel('Time (ps)',fontdict={'family':'Times New Roman', 'size' : 18})
        plt.ylabel(xlabel,fontsize=18,fontproperties='Times New Roman',backgroundcolor='w')
        plt.xticks(fontproperties='Times New Roman',fontsize=12)
        plt.yticks(fontproperties='Times New Roman',fontsize=12)
        plt.grid(True,linestyle = '--', linewidth = 0.5)
        plt.legend( loc=3, ncol=3, bbox_to_anchor=(0.1,1)) 
        fig = plt.gcf()
        plt.figure()
        plt.show()
        fig.savefig(outfolder+txtlabel+xcol[:5]+'.tif', bbox_inches='tight',dpi=900) 
    except Exception as err:
        print(err)
        



    

