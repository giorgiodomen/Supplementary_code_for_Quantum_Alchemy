from pyscf.gto.moleintor import getints
from pyscf.gto.mole import make_env
from pyscf.scf._vhf import incore
import numpy as np

from NucNuc_derivs import first_deriv_nuc_nuc
from alch_deriv import parse_charge,DeltaV
from FcMole import FcM,FracMole
from pyscf.data.elements import charge as Z_charge
import copy


def AC_S_g_h1(mol,pvec,acbsf):
    fc=mol.atom_charges()+np.asarray(pvec)
    at_ident=[A[0] for A in mol._atom]
    if len(set(at_ident))!=len(at_ident):
        for i in range(len(mol._atom)):
            new_name=mol._atom[i][0]+str(i)
            mol._atom[i]=(new_name,mol._atom[i][1])
    at_ident=[A[0] for A in mol._atom]
    new_bas={}
    for i,idf in enumerate(at_ident):
        new_bas[idf]=acbsf(fc[i])
    return S_g_h1(mol,new_bas,pvec)

def S_g_h1(mol,basis,pvec):
    a,b,e = make_env(mol._atom, basis, mol._env, mol.nucmod, mol.nucprop)
    a=np.array(a,dtype=np.float64)
    S=getints("int1e_ovlp_sph",a,b,e,hermi=1,aosym='s1')
    G=getints("int2e_sph",a,b,e,hermi=1,aosym='s8')
    h1_kin=getints("int1e_kin_sph",a,b,e,hermi=1,aosym='s1')    
    h1_nuc=np.zeros(h1_kin.shape)
    for i in range(len(a)):
        a[i][0]=0
    for i in range(len(pvec)):
        a[i][0]=1
        h1_nuc+=getints("int1e_nuc_sph",a,b,e,hermi=1,aosym='s1')*(mol.atom_charges()[i]+ pvec[i])
        a[i][0]=0
    a,b,e = make_env(mol._atom, basis, mol._env, mol.nucmod, mol.nucprop)
    return S,G,h1_kin,h1_nuc

def make_rdm1e(mo_energy, mo_coeff, mo_occ):
    '''Energy weighted density matrix'''
    mo0 = mo_coeff[:,mo_occ>0]
    mo0e = mo0 * (mo_energy[mo_occ>0] * mo_occ[mo_occ>0])
    return np.dot(mo0e, mo0.T.conj())


def AD1(mf,pvec,acbsf,dl=.02):
    pvec=np.array(pvec,dtype=np.float64)
    assert len(pvec)==len(mf.mol._atom)
    Sp,Gp,h1_kin_p,h1_nuc_p=AC_S_g_h1(mf.mol,pvec*dl,acbsf)
    Sm,Gm,h1_kin_m,h1_nuc_m=AC_S_g_h1(mf.mol,-pvec*dl,acbsf)
    dS,dG,dh1_kin,dh1_nuc=(Sp-Sm)/(2*dl),(Gp-Gm)/(2*dl),(h1_kin_p-h1_kin_m)/(2*dl),(h1_nuc_p-h1_nuc_m)/(2*dl) #cfd diff for the AO ints
    dH1_alch=DeltaV(mf.mol,parse_charge(pvec))[0]
    dVnn=first_deriv_nuc_nuc(mf.mol,parse_charge(pvec))
    P=mf.make_rdm1()
    C=mf.mo_coeff
    O=mf.mo_occ
    e=mf.mo_energy
    W=make_rdm1e(e,C,O)
    J,K=incore(dG,P)
    AD=[np.einsum("ij,ij",P,dh1_kin+dh1_nuc),dVnn,+0.5*np.einsum("ij,ij",P,J),-np.einsum("ij,ji",P,K)/4 ,-np.einsum("ij,ij",W,dS)]

    return(np.sum(AD))




def FcMfBS(atom,pvec,acbsf,**kwargs):
    at_ident=[a.split()[0] for a in atom.split(';')]
    at_ident_new=copy.deepcopy(at_ident)
    atom_new=""
    if len(set(at_ident))!=len(at_ident):
        for i in range(len(pvec)):
            at_ident_new[i]=at_ident[i]+str(i)
        atom_split=atom.split(';')
        for i,l in enumerate(atom_split):
            atom_new+=l.replace(at_ident[i],at_ident_new[i])+';'
    else: atom_new=atom
    up_ch=[Z_charge(i) for i in at_ident_new]
    new_bas={}
    fc=np.array(up_ch)+np.array(pvec)
    for i,idf in enumerate(at_ident_new):
        new_bas[idf]=acbsf(fc[i])
    return FcM(fcs=list(pvec),atom=atom_new,basis=new_bas,**kwargs)
    



