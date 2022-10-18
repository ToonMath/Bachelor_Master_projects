import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
import matplotlib.pyplot as plt
from random import randint
#define pauli matrix for Sz -------------------------------------
sigmaz = np.array([[1,0],[0,-1]])

#define operator matrices -------------------------------------
Sz = sigmaz/2
Splus = np.array([[0.,1.],[0.,0.]])
Smin = np.array([[0,0],[1,0]])
I = sp.identity(2)

#ddefine matrix and reshape ----------------------------------
H2 = sp.kron(Sz,Sz)+(sp.kron(Splus,Smin)+sp.kron(Smin,Splus))/2




# get Hamiltonian -----------------------------

def getsparseham(H2, N, bc = 'periodic'):
    "Construct nearest-neighbor Hamiltonian (sparse) on an N-site lattice"    
    d = int(np.sqrt(H2.shape[0])) # local dimension of a site
    # initialize sparse matrices    
    idn = sp.csr_matrix(1.)
    id1 = sp.csr_matrix(np.eye(d))
    H2 = sp.csr_matrix(H2)            
    Hk = H2    # initial term
    for k in range(N-2):
        # increase size of identity operator
        idn = sp.kron(idn, id1, 'csr')
        # add new terms to Hamiltonian
        Hk = sp.kron(Hk, id1, 'csr') + sp.kron(idn, H2, 'csr')    
    if bc == 'periodic':
        # create H2 x I x I x I x I
        Hl = sp.kron(H2, idn, 'csr')
        # and make a shift of the indices to get periodic term    
        ind = np.arange(d**N)
        iN = ind % d 
        ind = (ind - iN) / d + iN*d**(N-1)
        ind = ind.astype(int)
        Hl = Hl[ind,:][:,ind]
        # add term to Hamiltonian
        Hk = Hk + Hl
    return Hk

H = getsparseham(H2,16,'periodic')
eigval,eigvec = spl.eigsh(H,1,which='SA')
newsize2 = int(np.sqrt(eigvec[:,0].size))


# compute full decomposition ---------------------------

def fulldecomposition(H):   
    minsize = []
    #perform first split
    eigval,eigvec = spl.eigsh(H,1,which='SA')
    newsize2 = int(np.sqrt(eigvec[:,0].size))
    newground = np.reshape(eigvec[:,0],(newsize2,newsize2))
    svdR,svdλ,svdL = np.linalg.svd(newground,full_matrices=False)
    notminsize = [svdR,svdL]
    
    while svdR.shape!=(2,2) and len(notminsize)!=0:
        eigval,eigvec = spl.eigsh(svdR,1,which='SA')
        newsize2 = int(np.sqrt(eigvec[:,0].size))
        newground = np.reshape(eigvec[:,0],(newsize2,newsize2))
        svdR,svdλ,svdL = np.linalg.svd(newground,full_matrices=False)    
        H = svdR
        notminsize.append(svdL)
        if svdR.shape==(2,2):
            minsize.append(svdR)
            print(len(notminsize),'\n')
            if len(notminsize)>=1:
                if notminsize[-1].shape==(2,2):
                    minsize.append(notminsize[-1])
                    del notminsize[-1]
                    svdR = notminsize[-1]
                if notminsize[-1].shape!=(2,2):
                    svdR = notminsize[-1]
                    del notminsize[-1]
    return minsize


def statecomputation(MPSarray):
    state = MPSarray[0]
    size_of_state = 2**len(MPSarray)
    newstate = np.array((size_of_state,size_of_state))
    for i in range(1,len(MPSarray)-1):
        newstate = np.tensordot(state,MPSarray[i],axes=0) 
        state = newstate
    return state


minsize = fulldecomposition(H)
state = statecomputation(minsize)
print(state.reshape(H.shape))























