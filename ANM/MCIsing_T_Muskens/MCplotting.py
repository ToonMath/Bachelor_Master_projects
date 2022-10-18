import numpy as np
import re
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.cm as cm
import itertools

filename = "/Users/toonmuskens/Desktop/ANM/DataMc.csv"
filename1 = "/Users/toonmuskens/Desktop/ANM/MagDataMc.csv"
data = pd.read_csv(filename, sep = ' ')
data1 = pd.read_csv(filename1, sep=' ')
LData = np.array(data['L'].tolist())
LData2 = list( dict.fromkeys(LData) )
TData = np.array(data['T'].tolist())
mData = np.array(data['order'].tolist())
eData = np.array(data['error'].tolist())
tauData = np.array(data['tau'].tolist()) 
MData = np.array(data1['M'].tolist())
MeData = np.array(data1['error'].tolist())
nData =16

fig, ax = plt.subplots(2,1,figsize=(10,6))
fig2,ax2 = plt.subplots()



for i in range(len(LData2)):
    ax[0].plot(TData[nData*i:(i+1)*nData],mData[nData*i:(i+1)*nData],label=LData2[i])
    ax[0].set_xlabel('T')
    ax[0].set_ylabel('m')
    ax[0].errorbar(TData[nData*i:(i+1)*nData],mData[nData*i:(i+1)*nData] , eData[nData*i:(i+1)*nData],fmt='+')
    ax[0].legend()
    
    ax[1].plot(TData[nData*i:(i+1)*nData],tauData[nData*i:(i+1)*nData],label=LData2[i])
    ax[1].set_xlabel('T')
    ax[1].set_ylabel('t')
    ax[1].legend()

    ax2.plot(TData[nData*i:(i+1)*nData],MData[nData*i:(i+1)*nData],label=LData2[i])
    ax2.set_xlabel('T')
    ax2.set_ylabel('M')
    ax2.errorbar(TData[nData*i:(i+1)*nData],MData[nData*i:(i+1)*nData] , MeData[nData*i:(i+1)*nData],fmt='+')
    ax2.legend()
fig.savefig('/Users/toonmuskens/Desktop/ANM/exercise4_2.pdf')
fig2.savefig('/Users/toonmuskens/Desktop/ANM/exercise4_2b.pdf')