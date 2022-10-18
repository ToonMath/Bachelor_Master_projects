#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 10 11:12:34 2022

@author: toonmuskens

7.2 Exact diagonalization
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
import matplotlib.pyplot as plt


"-----------------------------Define rc variables--------------------"
plt.rcParams['text.usetex'] = True
plt.rcParams.update({
"text.usetex": True,
"font.family": "sans-serif",
"font.serif": ['Computer Modern'],
"font.size":30})
plt.rcParams["figure.figsize"] = (10,6.5)


"------------define filename to store results------------------------"

filename = "7_2_Spin-chain.txt"
f = open(filename, "w")

"-------- define S-operators as matrices for 2 site S=1/2-----------"

Sz      = np.array([[0.5,0],[0,-0.5]])
Splus   = np.array([[0,1],[0,0]])
Smin    = np.array([[0,0],[1,0]])
"-------------------------------------------------------------------"
"define function that returns hamiltonian for n-chain with open boundaries or periodic"

#define lowest number of sites with nonzero hamiltonian ------
H2 = sp.kron(Sz,Sz,'csr')+(sp.kron(Splus,Smin,'csr')+sp.kron(Smin,Splus,'csr'))/2
#create class to reuse calculated variables
class sparse:
    
    
    def __init__(self, Hold,n):
        self.Hold = Hold
        self.Energy = np.zeros((n-1,4))
        self.n = n
    def getEnergy(self):
    #define first nonzero hamiltonian for lowest number of sites
        I = np.identity(2)
        
        for i in range(2,self.n+1):
            if i==2:
                H = H2
                E,sweep = spl.eigsh(H,3,which ='SA')
                self.Energy[i-2] = [i,E[0],E[1],E[2]]
            if i>2:    
                In = np.identity(2**(i-2))
                Hnew = self.Hold
                H = sp.kron(In,H2,'csr')+sp.kron(Hnew,I,'csr')
                Hprime = sp.kron(Sz,sp.kron(In,Sz,'csr'),'csr')+(sp.kron(Splus,sp.kron(In,Smin,'csr'),'csr')+sp.kron(Smin,sp.kron(In,Splus,'csr'),'csr'))/2.
                self.Hold =  H + Hprime
                E,sweep = spl.eigsh(H,3,which ='SA')
                self.Energy[i-2] = [i,E[0],E[1],E[2]]
        return self.Energy
            
    
    def getEgap(self):
        E = self.Energy
        Egap = E[:,2]-E[:,1]
        return E[:,0],Egap
#-------------------EXERCISE 7.2.4--------------------------

key = sparse(H2,14)
E = key.getEnergy()
Egap = key.getEgap()
plt.plot(E[:,0],E[:,1])
# plt.plot(Egap[0],Egap[1])
#get ground energy data and create plot

#-------------------EXERCISE 7.2.5--------------------------   



    
    
    