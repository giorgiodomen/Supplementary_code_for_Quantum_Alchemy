from numpy.linalg import norm

def first_deriv_nuc_nuc(mol,dL):
    """dL=[[i1,i2,i3],[c1,c2,c3]]"""
    dnn=0
    for j in range(len(dL[0])):
        r2 = mol.atom_coord(dL[0][j]) 
        for i in range(mol.natm):
            if i != dL[0][j]:
                q1 = mol.atom_charge(i)
                r1 = mol.atom_coord(i)
                r = norm(r1-r2)
                dnn += (q1 * dL[1][j])/ r
    return dnn

def second_deriv_nuc_nuc(mol,dL):
    """dL=[[i1,i2,i3],[c1,c2,c3]]"""
    dnn=0
    for j in range(len(dL[0])):
        r2 = mol.atom_coord(dL[0][j]) 
        for i in range(len(dL[0])):
            if dL[0][i] > dL[0][j]:
                r1 = mol.atom_coord(dL[0][i])
                r = norm(r1-r2)
                dnn += (dL[1][i] * dL[1][j])/ r
    return 2*dnn