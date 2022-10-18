import numpy as np
import matplotlib.pyplot as plt
import math
import random

"""create code for Mcising model 4.1"""
#define start configuration wherep spin up=1 and spin down =-1

#get new configuration

def proposeFlip(array,L):
    newConfiguration = array
    x = np.random.randint(0,L)
    y = np.random.randint(0,L)
    newConfiguration[x][y] = (-1)*newConfiguration[x][y]
    Spin = [x,y]
    return newConfiguration,Spin

#get total energy
J = 1
def getTotalEnergy(c,L):
    energySum = 0.0
    for i in range(0,L):
        for j in range(0,L):
            energySum = energySum-J*c[i][j]*c[(i-1)%L][j]
            energySum = energySum-J*c[i][j]*c[(i+1)%L][j]
            energySum = energySum-J*c[i][j]*c[i][(j-1)%L]
            energySum = energySum-J*c[i][j]*c[i][(j+1)%L]
    return energySum/2
#get energy change
def getEnergyChange(c,randSpin,L):
    deltaE = 0
    deltaE = deltaE-c[randSpin[0]][randSpin[1]]*c[(randSpin[0]-1)%L][randSpin[1]]
    deltaE = deltaE-c[randSpin[0]][randSpin[1]]*c[(randSpin[0]+1)%L][randSpin[1]]
    deltaE = deltaE-c[randSpin[0]][randSpin[1]]*c[randSpin[0]][(randSpin[1]-1)%L]
    deltaE = deltaE-c[randSpin[0]][randSpin[1]]*c[randSpin[0]][(randSpin[1]+1)%L]
    return deltaE*2*J

# define magnetization
def getMagnetization(c,L):
    spinarray = np.reshape(c,(L**2,1))
    M = sum(spinarray)
    return M
def getMagnetizationChange(c,randSpin,L):
    return c[randSpin[0]][randSpin[1]]  

def binning_analysis(samples):
    """Perform a binning analysis over samples and return 
    errors: an array of the error estimate at each binning level, 
    tau: the estimated integrated autocorrelation time, 
    converged: a flag indicating if the binning has converged, and 
    bins: the last bin values"""
    minbins = 2**6 # minimum number of bins     
    maxlevel = int(np.log2(len(samples)/minbins)) # number of binning steps
    maxsamples = minbins * 2**(maxlevel)   # the maximal number of samples considered 
    bins = np.array(samples[-maxsamples:]) 
    errors = np.zeros(maxlevel+1)
    for k in range(maxlevel):
        errors[k] = np.std(bins)/np.sqrt(len(bins)-1.)
        bins = np.array((bins[::2] + bins[1::2])/2.)
        
    errors[maxlevel] = np.std(bins)/np.sqrt(len(bins)-1.)    
    tau = 0.5*((errors[-1]/errors[0])**2 - 1.)
    relchange = (errors[1:] - errors[:-1]) / errors[1:]
    meanlastchanges = np.mean(relchange[-3:])    # get the average over last changes
    converged = 1
    
    return (errors, tau, converged, bins)

   
def MCising(Nmc,Ns,temperatur,initMatrix,L):
    oldEnergy = getTotalEnergy(initMatrix,L)
    initMagnetization = getMagnetization(initMatrix,L)
    oldMatrix = initMatrix
    newMatrix = np.zeros((L,L))
    randSpin = []
    avgE = []
    avgM = []
    avgMsquare = []
    Order = []
    Nt = int(0.1*Nmc)
    Ntotal = int((Nmc+Nt)*L**2)
    T = temperatur
    for i in range(0,Ntotal):
        dE = 0
        dM = 0
        flipresult = proposeFlip(np.copy(oldMatrix),L)
        newMatrix = np.copy(flipresult[0])
        randSpin=(flipresult[1])
        #get difference in energy
        dE = getEnergyChange(newMatrix,randSpin,L) 
        newEnergy = oldEnergy+dE
        newMagnetization = getMagnetization(newMatrix,L)
        #plt.matshow(newMatrix)
        if newEnergy<=oldEnergy:
            oldMatrix = np.copy(newMatrix)
            oldEnergy = newEnergy
            initMagnetization = newMagnetization
            
        elif newEnergy>oldEnergy:
            rnum = random.uniform(0,1)
            if rnum<math.exp(-(1/T)*dE):
                oldMatrix = np.copy(newMatrix)
                oldEnergy = newEnergy 
                initMagnetization = newMagnetization
        if i%L**2 == 0:
          avgE.append(oldEnergy/L**2)
          avgM.append(initMagnetization/L**2)
          avgMsquare.append(initMagnetization**2)
          Order.append(abs(initMagnetization)/L**2)
    # remove thermalization step
    
    avgE = avgE[Nt:]
    avgM = avgM[Nt:]
    avgMsquare = avgMsquare[Nt:]
    Order = Order[Nt:]
    
    avgE = avgE[:Ns]
    avgM = avgM[:Ns]
    avgMsquare = avgMsquare[:Ns]
    Order = Order[:Ns]
    
    
    avg = [np.mean(avgE),np.mean(avgM),np.mean(avgMsquare),np.mean(Order)]
    
    Eerrors, Etau, Econv, Ebins = binning_analysis(avgE)
    Merrors, Mtau, Mconv, Mbins = binning_analysis(avgM)
    M2errors, M2tau, M2conv, M2bins = binning_analysis(avgMsquare)
    Oerrors, Otau, Oconv, Obins = binning_analysis(Order)
    
    Eerror = round(np.mean(Eerrors),6)
    Merror = round(np.mean(Merrors),6)
    M2error = round(np.mean(M2errors),6)
    Oerror = round(np.mean(Oerrors),6)
    
    error = [Eerror,Merror,M2error,Oerror]
    
    for i in range(len(avg)):
        if abs(avg[i])<error[i]:
            avg[i] = 0
            error[i] =0
    return avg[2],error[2],Otau  
#define paramters to use for MC simulation 

 
"""collect data for plots 4.2"""

l = [4,8,12,16,20,24]
T = np.linspace(2,3.4,16)
f = open('MagDataMc.csv','w+')
f.write("L T M error tau")
f.write('\n')
for i in range(len(l)):
    for j in range(len(T)):
        L = l[i]
        datalist = MCising(2**14,100,T[j-1],np.ones((L,L)),L)
        value = datalist[0]
        error = datalist[1]
        data = str(L)+' '+str(T[j])+' '+str(value)+' '+str(error)+' '+str(datalist[2])
        f.write(data)
        f.write('\n')
f.close()