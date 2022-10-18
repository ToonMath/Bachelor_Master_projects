
# Exercise 7.1-----------------------------

import numpy as np
import cmath


# Function for S=1/2 for 2 sites----------------------------
def spin_half_heisenberg():
    # define pauli matrices
    sigmax = np.array([[0,1],[1,0]])
    sigmay = np.array([[0,complex(0,1)],[complex(0, -1),0]])
    sigmaz = np.array([[1,0],[0,-1]])

    #define operator matrices for S = 1/2
    Sz = sigmaz/2
    Splus = [[0.,1.],[0.,0.]]
    Smin = [[0.,0.],[1.,0.]]

    #get hamiltonian and eigenvalue for S = 1/2 heisenberg model
    hz = np.kron(Sz,Sz)+(np.kron(Splus,Smin)+np.kron(Smin,Splus))/2

    eigval,sweepmatrix = np.linalg.eigh(hz)
    print('S=1/2 Heisenberg model\n---------------------\n')
    print('Hamiltonian is :\n',hz)
    print('With eigenvalue: \n',np.round(eigval,4),'and eigenvectors:\n',np.round(sweepmatrix,4),sep='\n')
 
    
# Function for S=1 for 2 sites----------------------------
def spin_one_heisenberg():
    # define operatros for S=1
    Sz = np.array([[1,0,0],[0,1,0],[0,0,1]])
    Splus = np.array([[0,1,0],[0,0,1],[0,0,0]])*2**0.5
    Smin = np.array([[0,0,0],[1,0,0],[0,1,0]])*2**0.5
    
    #get hamiltonian for S=1
    hz = np.kron(Sz,Sz)+(np.kron(Splus,Smin)+np.kron(Smin,Splus))/2
    
    #determine eigenvalues and eigenvectors
    eigval,sweepmatrix = np.linalg.eigh(hz)
    print('\nS=1 Heisenberg model\n-------------------\n')
    print('Hamiltonian is :\n',hz)
    print('With eigenvalue: \n',np.round(eigval,4),'and eigenvectors:\n',np.round(sweepmatrix,4),sep='\n')

spin_half_heisenberg()
spin_one_heisenberg()

