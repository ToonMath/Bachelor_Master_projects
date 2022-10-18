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
    

H16 = getsparseham(H2,16,'periodic')
eigval,eigvec = spl.eigsh(H16,1,which='SA')
newsize2 = int(np.sqrt(eigvec[:,0].size))


#reshape groundstate and split ---------------------------
newground = np.reshape(eigvec[:,0],(newsize2,newsize2))
svd = np.linalg.svd(newground,full_matrices=False)

S = - np.dot(svd[1]**2,np.log(svd[1]**2))
fig,ax = plt.subplots(2,1,figsize=(10,6))
sumeig = np.sum(svd[0]**2)
normconst = 1/sumeig

ax[0].semilogy(svd[1]**2,label='ground state')
ax[0].set_ylabel('$p_k$')
ax[0].set_title('probability values ($p_k$)')

ax[1].semilogy(normconst*svd[1]**2,label='ground state')
ax[1].set_xlabel('n')
ax[1].set_ylabel('$p_k$')
ax[1].set_title('probability values ($p_k$) normalised')


#--------- create random state------------

randomstate = np.random.rand(newsize2,newsize2)
randeigvec = randomstate.reshape((newsize2**2,1))
svdrand = np.linalg.svd(randomstate,full_matrices=False)
Srand = - np.dot(svdrand[1]**2,np.log(svdrand[1]**2))
randeigval = np.sqrt(np.sum(H16.dot(randeigvec)**2)/np.sum(randeigvec**2))

sumeig = np.sum(svdrand[0]**2)
normconst = 1/sumeig
ax[0].semilogy(svdrand[1]**2,label='random matrix')
ax[1].semilogy(normconst*svdrand[1]**2,label='random matrix')

ax[0].legend()
ax[1].legend()

fig.savefig('Eigvalues_densMatrix.pdf')

#------------ energy error -------------

svd_ground_a = np.linalg.eig(svd[0])[0][0]
svd_ground_b = np.linalg.eig(svd[2])[0][0]
Eground = svd[1][2:40]*svd_ground_a*svd_ground_b

svdrand_ground_a = np.linalg.eig(svdrand[0])[0][0]
svdrand_ground_b = np.linalg.eig(svdrand[2])[0][0]
Eground_rand = svd[1][2:40]*svdrand_ground_a*svdrand_ground_b

error = -1*Eground+eigval
error_rand = -1*Eground_rand+randeigval
norm_error = np.sum(np.abs(error)**2)
norm_error_rand = np.sum(np.abs(error_rand)**2) 

fig1,ax1 = plt.subplots(figsize=(10,6))
ax1.semilogy((np.abs(error)**2)/norm_error,label = 'ground')
ax1.semilogy((np.abs(error_rand)**2)/norm_error_rand,label = 'random')
ax1.legend()

fig1.savefig('Energyerror.pdf')



