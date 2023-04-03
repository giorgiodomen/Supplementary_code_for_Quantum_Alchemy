from pyscf.gto import M
import numpy as np
import basis_set_exchange as bse
import re
import sys

atoms="G,H,He,Li,Be,B,C,N,O,F,Ne,Na,Mg,Al"
atoms=atoms.split(',')

def abse_atom(ref,targ,method,xcf,bs="pcX-2"):
    if ref==targ: return 0
    spin=(atoms.index(targ))%2
    if '# Basis Set Exchange' in bs:
        for lin in bs.split('\n')[:10]:
            if "#   Basis set:" in lin:
                bs=lin.split(':')[-1].strip()
                break
    T=M(atom='{} 0 0 0'.format(targ),spin=spin,\
             basis=bse.get_basis(bs,fmt="nwchem",elements=[atoms.index(targ)]),verbose=0)
    TatR=M( atom='{} 0 0 0'.format(targ),spin=spin,\
             basis=bse.get_basis(bs,fmt="nwchem",elements=[atoms.index(ref)]),verbose=0)
    mf=method(T)
    if xcf: mf.xc=xcf
    eT=mf.scf()
    mf=method(TatR)
    if xcf: mf.xc=xcf
    eTatR=mf.scf(dm0=mf.init_guess_by_1e())
    return eT-eTatR

def absec(ref,targ,method,xcf,bs="pcX-2"):
    reflist= re.sub( r"([A-Z])", r" \1", ref).split()
    targlist= re.sub( r"([A-Z])", r" \1", targ).split()
    if len(reflist) != len(targlist):
        print(reflist,targlist,"reference and target lengths do not match!", sys.exc_info()[0])
        raise 
    bsae=0
    for i in range(len(reflist)):
        bsae+=abse_atom(reflist[i],targlist[i],method,xcf,bs)
    return (bsae)