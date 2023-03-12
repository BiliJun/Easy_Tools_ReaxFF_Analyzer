# -*- coding: utf-8 -*-
"""


"""
lmpstring ='''
Temp = Temperature (K)
PotEng = Potentiol energy (kcal/mol)
KinEng = Kinetic energy (kcal/mol)
E_pair = E_pair
TotEng =Total energy (kcal/mol)
Press = Pressure 
Density = Density (g/cm^-3)
Vol = Volume
vol = Volume
dt = Timestep size
temp = Temperature
press = Pressure
pe = Total potential energy
ke = Kinetic energy
eb = Bond energy
ea = Atom energy
elp = Lone-pair energy
emol = Molecule energy
ev = Valence angle energy
epen = double-bond valence angle penalty
ecoa = Valence angle conjugation energy
ehb = Hydrogen bond energy
et = Torsion energy
eco = Conjugation energy
ew = Van der Waals energy
ep = Coulomb energy
efi = Electric field energy
eqeq = Charge equilibration energy
etotal = Total energy
evdwl = van der Waals pairwise energy
ecoul = Coulombic pairwise energy
epair = pairwise energy
ebond = Bond energy
eangle = Angle energy
edihed = Dihedral energy
eimp = Improper energy
elong = Long-range kspace energy
etail = Van der Waals energy long-range tail correction
enthalpy = Enthalpy
ecouple = Cumulative energy change due to thermo/baro statting fixes
'''
def lmpstrs () :
    lmpstr = lmpstring.split('\n')
    lmpstr = [s.split('=') for s in lmpstr if s != '']
    lmpstrs = {}
    for s in lmpstr :
        if s[0][0] == 'e' :
            lmpstrs['v_'+s[0].strip()]=s[1].strip()
        else :
            lmpstrs[s[0].strip()] = s[1].strip()    
    return lmpstrs




