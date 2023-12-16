# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 13:13:47 2023

@author: zhouw
"""

import requests
import pubchempy
import pandas as pd
from tqdm import tqdm
CACTUS = "https://cactus.nci.nih.gov/chemical/structure/{0}/{1}"


def smiles_to_iupac(smiles):
    rep = "iupac_name"
    url = CACTUS.format(smiles, rep)
    response = requests.get(url)
    response.raise_for_status()
    return response.text

chem, chF = [], []
pf = pd.read_csv(r'D:\zh\1.csv')
for smi in tqdm(pf['0']) :
    pass
    try :
        # Eng name
        chem.append(smiles_to_iupac(smi))
        smi = 'OCC[C]1O[C@H]([C@H]([C@H]1O)O)CO'
        # 系统命名法 
        compounds = pubchempy.get_compounds(smi, namespace='smiles')
        match = compounds[0]
        chF.append(match.iupac_name)
    except:
        chem.append('error')
        chF.append('error')
pf['chem'] = chem
pf['sysChem'] = chF

pf.to_csv(r'D:\zh\2000chem.csv')
