import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
import matplotlib.pyplot as plt


#define pauli matrices -------------------------------------
sigmaz = np.array([[1,0],[0,-1]])

#define operator matrices for S = 1/2
Sz = sigmaz/2
Splus = [[0.,1.],[0.,0.]]
Smin = [[0.,0.],[1.,0.]]
I = np.identity(2)


#define matrices for 2 and 3 sites-------------------------
H2 = np.kron(Sz,Sz)+(np.kron(Splus,Smin)+np.kron(Smin,Splus))/2
H3 = np.kron(H2,I)+np.kron(I,H2)


#------------------EXERCISE 7.2.1--------------------------
H4 = np.kron(H3,I)+np.kron(np.kron(I,I),H2)
eigval,sweepmatrix = np.linalg.eigh(H4)
groundE_open = eigval[0]
print('Ground state energy open boundaries:\t',groundE_open)

#------------------EXERCISE 7.2.2--------------------------
Hprime = np.kron(H2,I)+np.kron(I,H2)
Ip = np.kron(I,I)
H4 = np.kron(Hprime,I)+np.kron(Ip,H2)+np.kron(np.kron(Sz,Ip),Sz)+0.5*(np.kron(np.kron(Splus,Ip),Smin)+np.kron(np.kron(Smin,Ip),Splus))
eigval,sweepmatrix = np.linalg.eigh(H4)
groundE_periodic = eigval[0]
print('Ground state energy for periodic boundaries:\t',groundE_periodic)

#-----------------------------------------------------------
#-------------------EXERCISE 7.2.3--------------------------
Szs = sp.csr_matrix(Sz)
Spluss = sp.csr_matrix(Splus)
Smins = sp.csr_matrix(Smin)
H2s = sp.kron(Szs,Szs,'csr')+(sp.kron(Spluss,Smins,'csr')+sp.kron(Smins,Spluss,'csr'))/2
H3s = sp.kron(H2,I)+sp.kron(I,H2)

H4s = sp.kron(H3s,I,'csr')+sp.kron(sp.kron(I,I,'csr'),H2s,'csr')

groundE_sparse,sweep = spl.eigsh(H4s, 4, which='SA')
print('Ground state energy open boundaries sparse matrix:\t',groundE_sparse[0])

#----- check memory usage--------
densmem = H4.nbytes
sparsemem = H4s.data.nbytes
print('Memory used by dense matrix: %d bytes'%densmem,'\nMemory used by sparse matrix: %d bytes'%sparsemem)


















  
        