import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
import matplotlib.pyplot as plt


#define pauli matrix for Sz -------------------------------------
sigmaz = np.array([[1,0],[0,-1]])

#define operator matrices -------------------------------------
Sz = sigmaz/2
Splus = np.array([[0.,1.],[0.,0.]])
Smin = np.array([[0,0],[1,0]])
I = np.array([[1,0],[0,1]])

#ddefine matrix and reshape ----------------------------------
H2 = np.kron(Sz,Sz)+(np.kron(Splus,Smin)+np.kron(Smin,Splus))/2
NewH2 = np.reshape(H2,(2,2,2,2))

print(H2,NewH2,sep='\n\n')

#find eigenvalues and eigenvectors and reshape--------------------------
eigval,eigvec = np.linalg.eigh(H2)

eigvec0 = np.reshape(eigvec[0],(2,2))
print('Find eigenvector and reshape')
print(eigvec[0],eigvec0,sep='\n\n')
#calculate tensorproduct <ground|H|ground> and <ground|ground>
Tprod = np.tensordot(np.tensordot(eigvec0, NewH2), eigvec0)
Tground = np.tensordot(eigvec0,eigvec0)
print(Tprod,Tground)
#---------- get transposed ground state matrice and compute again
newground = np.transpose(eigvec0)
Tprod = np.tensordot(np.tensordot(newground, NewH2), newground)
Tground = np.tensordot(newground,newground)
print('Using transposed ground state matrix')
print(Tprod,Tground)

print('\nCheck if matrix is hermitian, if this is the case than H-H^T = 0')
print(H2,np.transpose(H2),sep='\n\n')
