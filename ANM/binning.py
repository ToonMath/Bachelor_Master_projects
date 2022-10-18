import numpy as np

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
    if meanlastchanges > 0.05:
        print("warning: binning maybe not converged, meanlastchanges:", meanlastchanges)
        converged = 0
    return (errors, tau, converged, bins)

