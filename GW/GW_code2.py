import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.constants


M_sun = 1.9891*10**30
#import data and devide in corresponding variables
filename = '/Users/toonmuskens/Downloads/strain_data_exec3.dat'
data = pd.read_csv(filename, sep = ' ')
time = np.array(data['time'].tolist())
h_plus = np.array(data['strain_plus'].tolist())
h_cross = np.array(data['strain_cross'].tolist())
#plot polerizations
fig, ax = plt.subplots(figsize=(10,5))
ax.set_xlim([time[0], time[len(time)-1]])
ax.plot(time,h_plus,'r',label ='$h_+$(t)')
ax.plot(time,h_cross,'b',label ='$h_x$ (t)')
ax.set_title('Amplitude as function of time')
leg = ax.legend();

#find the instantaneous frequency
h_plus_square = h_plus**2
h_cross_square = h_cross**2

freq = []
time_2 = []
i=0
freq_i = 0
fig1, ax1 = plt.subplots(2,1,figsize=(10,6))
while time[i]<0:
    h_plus_derivative = (h_plus[i+1]-h_plus[i])/(time[i+1]-time[i])
    h_cross_derivative = (h_cross[i+1]-h_cross[i])/(time[i+1]-time[i])
    freq_i = (2*np.pi*(h_plus_square[i]+h_cross_square[i]))**(-1)*(h_cross_derivative*h_plus[i]-h_plus_derivative*h_cross[i])
    freq.append(float(freq_i))
    time_2.append(float(time[i]))
    i+=1
#plot frequency
ax1[0].set_xlim([-0.5, time[len(time)-1]])
ax1[0].plot(time_2,freq)

ax1[0].set_xlabel('t')
ax1[0].set_ylabel('f(t)[Hz]')
ax1[0].set_title('instantaneous frequency till t=0')

fig.savefig('/Users/toonmuskens/Desktop/Frequency_evo.pdf')

#find chirp mass
temp = 0
inv_freq = []
Mc = [] 
for j in range(i):
    inv_freq.append((freq[j]*8*np.pi)**(-8/3)*5)
    if j>=1:
        Mc.append((-(1)*(inv_freq[j]-inv_freq[j-1])/(time_2[j]-time_2[j-1]))**(3/5)*scipy.constants.c**3*scipy.constants.gravitational_constant**(-1)/M_sun)
time_2 = time[:len(time_2)-1]
ax1[1].set_xlim([-0.5, time[len(time)-1]])
ax1[1].plot(time_2,Mc)

ax1[1].set_xlabel('t')
ax1[1].set_ylabel('Mc[$M_{sol}$]')
ax1[1].set_title('Chirp mass')


fig1.savefig('/Users/toonmuskens/Desktop/chirp_mass_est.pdf')



