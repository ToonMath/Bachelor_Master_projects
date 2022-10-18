#
#


#-----------------------Exercise 7.2.6------------------------
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
import matplotlib.pyplot as plt


#define operator matrices -------------------------------------
Sz = np.array([[1,0,0],[0,0,0],[0,0,-1]])
Splus = np.array([[0,1,0],[0,0,1],[0,0,0]])*(2**0.5)
Smin = np.array([[0,0,0],[1,0,0],[0,1,0]])*(2**0.5)
I = np.array([[1,0,0],[0,1,0],[0,0,1]])

#define lowest number of sites with nonzero hamiltonian ------
H2 = sp.kron(Sz,Sz)+(sp.kron(Splus,Smin)+sp.kron(Smin,Splus))/2
E0,eigvec = spl.eigsh(H2,3,which ='SA')
#create class to reuse calculated variables
class variables:
    def __init__(self, Hprime, Ip,E0):
        self.Hprime = Hprime
        self.Ip = Ip
        self.Energy = [E0]
    def getEigenval(self):
    #define first nonzero hamiltonian for lowest number of sites
        Hprime = self.Hprime
        Ip = self.Ip
        H = sp.kron(Hprime,I,'csr')+sp.kron(Ip,H2,'csr')+sp.kron(sp.kron(Sz,Ip,'csr'),Sz)+0.5*(sp.kron(sp.kron(Splus,Ip,'csr'),Smin,'csr')+sp.kron(sp.kron(Smin,Ip,'csr'),Splus,'csr'))
        self.Hprime = sp.kron(Hprime,I,'csr')+sp.kron(Ip,H2,'csr')
        self.Ip = sp.kron(Ip,I,'csr')
        Energy,sweep = spl.eigsh(H,3,which ='SA')
        self.Energy.append(Energy)
#-------------------Plot ground state energy--------------------------



def getEnergy(N):
    var = variables(H2,I,E0)
    for j in range(0,len(N)-1):
        var.getEigenval() 
    return var.Energy

#get ground energy data and create plot
def plotgroundEnergy(N):
    sites=np.arange(2,N+1,1)
    results = np.array(getEnergy(sites))
    plt.xlabel('N [number of sites]')
    plt.ylabel('$E_{ground}$')
    plt.title('Ground energy as function of number of sites')
    plt.plot(sites,results[:,0])
    plt.savefig('GroundEnergy_spin_one.pdf')
#-------------------Plot Energy gap--------------------------   
def gapEnergy(N):
    sites = np.arange(2,N+1,1)
    energy = np.array(getEnergy(sites))
    energygap = energy[:,1]-energy[:,0]
    return sites,energygap

def plotEnergyGap(N):
    results = gapEnergy(N)
    devide = np.ones(len(results[0][2::2]))*-1
    plt.xlabel('1/N [$sites^{-1}$]')
    plt.ylabel('$\Delta E$')
    plt.title('Energy gap per sites')
    xvalues,yvalues = results[1][2::2],results[0][2::2]**devide
    plt.plot(xvalues,yvalues,label='data')
    #interpolate to Egap = 0 ---------------------------------
    extrpol = np.polyfit(xvalues[-4:],yvalues[-4:],1)
    x = np.linspace(0,xvalues[0],5)
    y = extrpol[1]+x*extrpol[0]
    label = 'y='+str(round(extrpol[1],3))+'+'+str(round(extrpol[0],3))+'x'
    plt.plot(x,y,'--',label=label)
    plt.legend()
    plt.savefig('Energygap_spin_one.pdf')
  
def createAndSaveplots(N):
    plotgroundEnergy(N)
    plt.clf()
    plotEnergyGap(N)
    print('Finished')
    
createAndSaveplots(10)    
    
    
    