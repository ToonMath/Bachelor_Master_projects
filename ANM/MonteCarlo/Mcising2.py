import numpy as np
import matplotlib.pyplot as pl
import binning
import math
from numba import jit

def plotlatt(latt):
    """ Plot the lattice """
    pl.figure()
    pl.matshow(latt)
        
    
def getenergy(latt):
    """ Compute the energy of a configuration """
    # efficient code based on numpy-array operations
    energy = - np.sum(np.sum(latt * (np.roll(latt, 1, axis=0) + np.roll(latt, 1, axis=1))))      
    # Less efficient code based on for-loop
    #L = latt.shape[0]
    #energy = 0 
    #for i in range(L):
    #   for j in range(L):
    #        sij = latt[i,j]
    #        val = latt[(i+1)%L, j] + latt[i,(j+1)%L]             
    #        energy += -val*sij
    return energy

# Optional - improves the speed of the loop 
@jit()       
def doMCsweep(latt, energy, mag, beta):
    """ Do a MC sweep and return current energy and magnetization"""
    L = latt.shape[0]
    N = L*L
    # do one sweep (N MC steps)
    for k in range(N):
        i = np.random.randint(L)
        j = np.random.randint(L)
        sij = latt[i,j]
        val = latt[(i+1)%L, j] + latt[i,(j+1)%L] + latt[(i-1)%L,j] + latt[i,(j-1)%L]
        ediff = 2*val*sij
        prob = math.exp(-beta*ediff)
        if (np.random.random() < prob):  
            latt[i,j] = - latt[i,j]
            energy += ediff
            mag -= 2*sij
    return (energy, mag)

def getCluster(latt2,beta): 
        cluster = []
        cluster.append((np.random.randint(L),np.random.randint(L)))
        initSpin = latt2[cluster[0]]
        lat2 = latt2.copy()
        neighbor = [latt2[(cluster[0][0]+1)%L,(cluster[0][1])%L],latt2[(cluster[0][0])%L,(cluster[0][1]+1)%L],latt2[(cluster[0][0]-1)%L,(cluster[0][1])%L],latt2[(cluster[0][0])%L,(cluster[0][1]-1)%L]]
        neighborsites = [[(cluster[0][0]+1)%L,(cluster[0][1])%L],[(cluster[0][0])%L,(cluster[0][1]+1)%L],[(cluster[0][0]-1)%L,(cluster[0][1])%L],[(cluster[0][0])%L,(cluster[0][1]-1)%L]]
        
        for i in range(len(neighbor)):
            if neighbor[i] == initSpin:
                prob = 1-math.exp(-beta*2)
                if (np.random.random() < prob):
                    cluster.append(neighborsites[i])
                    latt2[neighborsites[i]]*-2
        return latt2.copy()

# A class to evaluate a measurement
class Obsresult:
    """Class to evaluate the measurements"""
    def __init__(self, vals, initSpin, name="obs"):
        self.mean = np.mean(vals)
        self.name = name
        # print(vals)
        (self.binerrors, self.tau, self.converged, self.binvals) = binning.binning_analysis(vals)
        self.err = self.binerrors[-1]
        self.initSpin = initSpin
        print(name, ":", self.mean, " +/- ", self.err)
    def plotbinning(self):
        pl.figure()
        pl.plot(self.binerrors,'ko')
        pl.xlabel('binning step')
        pl.ylabel('estimated error')
    def relerrors(self):
        return (self.binerrors[1:] - self.binerrors[:-1]) / self.binerrors[-1]


def runMCsim(L, T, Nsweeps, Ntherm, seed):    
    """run a Monte Carlo simulation"""
    N = L*L    
    beta = 1./T
    np.random.seed(seed)
    
    
    
    latt = np.random.rand(L,L)
    latt = np.sign(latt-0.5)
    latt = getCluster(latt, beta)
    # create lattice
    # get initial values
    energy = getenergy(latt)
    mag = np.sum(latt)  

    # store measurements in:
    venergy = np.zeros((Nsweeps, 1))
    vmag = np.zeros((Nsweeps, 1))

    print(('Do 2D Ising simulation with L=%i, T=%f, Ntherm=%i, Nsweeps=%i, seed=%i' \
    % (L, T, Ntherm, Nsweeps, seed)))
    

    print(('do %i thermalization sweeps...' % Ntherm))
    for k in range(Ntherm):
        latt = getCluster(latt, beta)
        (energy, mag) = doMCsweep(latt, energy, mag, beta)
        if k%1000==0:
            plotlatt(latt)
        
 
    print('done!')



    print(('perform %i MC sweeps and do measurements...' % Nsweeps))
    for k in range(Nsweeps):
        latt = getCluster(latt, beta)
        (energy, mag) = doMCsweep(latt, energy, mag, beta)
        if k%10000==0:
            plotlatt(latt)
        # now store values
        venergy[k] = energy
        vmag[k] = mag
        
    print('done!')
  
    # create results
    print("--------------------")
    print("Evaluate observables")
    rese = Obsresult(venergy/N, "Energy per site")
    resmag = Obsresult(vmag/N, "Magnetization per site")
    resm = Obsresult(np.abs(vmag)/N, "m")
    resE2 = Obsresult(venergy*venergy, "Energy squared")
    resM2 = Obsresult(vmag*vmag, "Magnetization squared")
    resM4 = Obsresult(vmag*vmag*vmag*vmag, "Magnetization fourth")
    
    # put all the results in a dictionary and return all results
    res = {"Es" : rese, "Ms" : resmag, "m" : resm, \
    "E2": resE2, "M2": resM2,"M4": resM2}
    print("--------------------\n")
    return res
    



if __name__ == "__main__":
    # This is a test run
    L = 24               # linear dimension of the lattice, lattice-size= N x N
    T = 2.              # temperature
    Nsweeps = 2**17     # total number of Monte Carlo steps
    Ntherm = 5000      # Number of thermalization sweeps
    seed = 100
    
    res = runMCsim(L, T, Nsweeps, Ntherm, seed)
