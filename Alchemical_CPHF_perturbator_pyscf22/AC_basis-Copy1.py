import basis_set_exchange as bse
import scipy as sp
from pyscf import gto,scf
import copy
import numpy as np
from scipy.interpolate import CubicSpline,make_interp_spline as MIS
from pyscf.data.elements import _symbol


zvals510=np.arange(5,11)
zvals310=np.arange(3,11)


#pcX-1 uncontracted BS constant contraction B-Ne
coeffs_pcx=[np.array([e[1][0] for e in \
      gto.basis.load(bse.get_basis("pcX-1",fmt="nwchem",elements=[int(z)]),_symbol(int(z)))]) for z in zvals510]
cspl_pcx=CubicSpline(zvals510,coeffs_pcx)
def pcX1(z):
    fb=gto.basis.load(bse.get_basis("pcX-1",fmt="nwchem",elements=[int(z)]),_symbol(int(z)))
    if z<=4 or z>=12:
        return fb
    for i in range(len(fb)):
        fb[i][1][0]=cspl_pcx(z)[i]  
    return (fb)


# 3-21 G - 6-31G pople's double Z constant contraction Li-Ne

coeffs_321=[np.array([g for ao in \
      gto.basis.load(bse.get_basis("3-21g",fmt="nwchem",elements=[int(z)]),_symbol(int(z))) \
          for g in ao[1:]]) for z in zvals310]
cspl_321=CubicSpline(zvals310,coeffs_321)
def p321g(z):
    fb=gto.basis.load(bse.get_basis("3-21g",fmt="nwchem",elements=[int(z)]),_symbol(int(z)))
    if z<=4 or z>=12:
        return fb
    ci=0
    for i in range(len(fb)):
        for j in range(1,len(fb[i])): 
            fb[i][j]=list(cspl_321(z)[ci])
            ci+=1
    return (fb)
coeffs_631=[np.array([g for ao in \
      gto.basis.load(bse.get_basis("6-31g",fmt="nwchem",elements=[int(z)]),_symbol(int(z))) \
          for g in ao[1:]]) for z in zvals310]
cspl_631=CubicSpline(zvals310,coeffs_631)
def p631g(z):
    fb=gto.basis.load(bse.get_basis("6-31g",fmt="nwchem",elements=[int(z)]),_symbol(int(z)))
    if z<=3 or z>=12:
        return fb
    ci=0
    for i in range(len(fb)):
        for j in range(1,len(fb[i])): 
            fb[i][j]=list(cspl_631(z)[ci])
            ci+=1
    return (fb)


# def2-SVP  / def2-TZVP constant contraction Li-Be / B-Ne
coeffs_svp=[np.array([g for ao in \
      gto.basis.load(bse.get_basis("def2-SVP",fmt="nwchem",elements=[int(z)]),_symbol(int(z))) \
          for g in ao[1:]]) for z in zvals510]

cspl_svp=CubicSpline(zvals510,coeffs_svp)
def def2_svp(z):
    fb=gto.basis.load(bse.get_basis("def2-SVP",fmt="nwchem",elements=[int(np.round(z,decimals=0))]),\
                     _symbol(int(np.round(z,decimals=0))))
    if z<=4 or z>=12:
        return fb
    ci=0
    for i in range(len(fb)):
        for j in range(1,len(fb[i])): 
            fb[i][j]=list(cspl_svp(z)[ci])
            ci+=1
    return (fb)

coeffs_tzvp=[np.array([g for ao in \
      gto.basis.load(bse.get_basis("def2-TZVP",fmt="nwchem",elements=[int(z)]),_symbol(int(z))) \
          for g in ao[1:]]) for z in zvals510]
cspl_tzvp=CubicSpline(zvals510,coeffs_tzvp)
def def2_tzvp(z):
    fb=gto.basis.load(bse.get_basis("def2-TZVP",fmt="nwchem",elements=[int(np.round(z,decimals=0))]),\
                     _symbol(int(np.round(z,decimals=0))))
    if z<=4 or z>=12:
        return fb
    ci=0
    for i in range(len(fb)):
        for j in range(1,len(fb[i])): 
            fb[i][j]=list(cspl_tzvp(z)[ci])
            ci+=1
    return (fb)



# cc-pvdz constant contraction Li-Ne
coeffs_pvdz=[[ao[1:] for ao in \
      gto.basis.load(bse.get_basis("cc-pvdz",fmt="nwchem",elements=[int(z)]),_symbol(int(z))) ] for z in zvals310]
coeffs_pvdz_L =[np.array([l[0] for l in coeffs_pvdz]),\
                np.array([l[1] for l in coeffs_pvdz]),\
                np.array([l[2] for l in coeffs_pvdz])]
cspl_pvdz=[CubicSpline(zvals310,c) for c in coeffs_pvdz_L]
def ccpvDZ(z):
    fb=gto.basis.load(bse.get_basis("cc-pvdz",fmt="nwchem",elements=[int(np.round(z,decimals=0))]),\
                 _symbol(int(np.round(z,decimals=0))))
    if z<=2 or z>=11:
        return fb
    l=0
    fb=[]
    for i in range(len(cspl_pvdz)):
        fb.append([i])
        for y in cspl_pvdz[i](z):
            fb[-1].append(list(y) )
    return (fb)


        # cc-pvtz  constant contraction B-Ne
coeffs_pvtz=[[ao[1:] for ao in \
      gto.basis.load(bse.get_basis("cc-pvtz",fmt="nwchem",elements=[int(z)]),_symbol(int(z))) ] for z in zvals510]
coeffs_pvtz_L =[np.array([l[0] for l in coeffs_pvtz]),\
                np.array([l[1] for l in coeffs_pvtz]),\
                np.array([l[2] for l in coeffs_pvtz]),\
                np.array([l[3] for l in coeffs_pvtz])]
cspl_pvtz=[CubicSpline(zvals510,c) for c in coeffs_pvtz_L]
def ccpvTZ(z):
    fb=gto.basis.load(bse.get_basis("cc-pvtz",fmt="nwchem",elements=[int(np.round(z,decimals=0))]),\
                 _symbol(int(np.round(z,decimals=0))))
    if z<=4 or z>=11:
        return fb
    l=0
    fb=[]
    for i in range(len(cspl_pvtz)):
        fb.append([i])
        for y in cspl_pvtz[i](z):
            fb[-1].append(list(y) )
    return (fb)

