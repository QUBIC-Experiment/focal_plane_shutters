import numpy as np
import matplotlib.pyplot as plt
import string 
from scipy import optimize
import scipy.signal as signal
from numpy import pi, r_
import csv




data_pre = np.loadtxt("/Users/albertopolini/Desktop/filescicloswitch1e2/180419082413_CAsw01n00trA_pre.txt")
data_post = np.loadtxt("/Users/albertopolini/Desktop/filescicloswitch1e2/180419082413_CAsw01n00trA_post.txt")

startdata = 0
enddata = 5000
dt = 1e-6

x = np.arange(0, (enddata-startdata)*dt,dt)
y_pre = data_pre [startdata:enddata,0]
z_pre = data_pre [startdata:enddata,1]
t_pre = data_pre [startdata:enddata,2]

y_post = data_post [startdata:enddata,0]
z_post = data_post [startdata:enddata,1]
t_post = data_post [startdata:enddata,2]

#print len(x), len(y_pre)
def analisi (x,y):
    #Creo il filtro FIR
    nyq_rate = len(x)/2
    cut_off = 60.
    N = 100
    taps = signal.firwin (N, cut_off/nyq_rate)

    #Filtro x 
    fy = signal.lfilter (taps, 1.0, y)

    #Trovo la variazione di fase del segnale originale
    shift = (0.5*(N-1))*dt


    
    #Derivo e plotto la funzione filtrata
    der = np.diff(fy[N-1:])/dt
    

    #Fit col metodo dei minimi quadrati
    fitfunc = lambda p, x: p[0]*np.cos(2*np.pi/0.001024*x+p[1])  # Target function
    errfunc = lambda p, x, y: fitfunc(p, x[N:]) - der # Distance to the target function
    p0 = [5000., 0.] # Initial guess for the parameters
    p1, success = optimize.leastsq(errfunc, p0[:], args=(x, der), factor=1, maxfev = 100000)

    return p1[0],p1[1]
    
'''prova1, prova2 = analisi(x,y_pre)
print prova2'''

amp1_pre,Phi1_pre = analisi(x,y_pre)
amp2_pre, Phi2_pre = analisi(x,z_pre)
amp3_pre, Phi3_pre = analisi(x,t_pre)

amp1_post, Phi1_post =analisi (x,y_post)
amp2_post, Phi2_post =analisi (x,z_post)
amp3_post, Phi3_post =analisi (x,t_post)

dPhi12_20 = np.abs(Phi1_pre-Phi2_pre)
dPhi13_20 = np.abs(Phi1_pre-Phi3_pre)

dPhi12_90 = np.abs(Phi1_post-Phi2_post)
dPhi13_90 = np.abs(Phi1_post-Phi3_post)

variazione1 = np.abs(dPhi12_20-dPhi12_90)
variazione2 = np.abs(dPhi13_20-dPhi13_90)

print "differenza di fase colonna 1-2 20mA, 90mA = %.4f rad" %variazione1
print "differenza di fase colonna 1-3 20mA, 90mA = %.4f rad" %variazione2

f1_pre = amp1_pre*np.cos(2*np.pi/0.001024*x+Phi1_pre) 
f1_post = amp1_post*np.cos(2*np.pi/0.001024*x+Phi1_post)


plt.figure(1)
plt.plot(x,f1_pre,"b")
plt.plot(x,f1_post*10,"r")
plt.grid(True)
plt.title("colonna 1, ampiezza 250mA aumentata di un fattore 10")

f2_pre = amp2_pre*np.cos(2*np.pi/0.001024*x+Phi2_pre) 
f2_post = amp2_post*np.cos(2*np.pi/0.001024*x+Phi2_post)


plt.figure(2)
plt.plot(x,f2_pre,"b")
plt.plot(x,f2_post*10,"r")
plt.grid(True)
plt.title("colonna 2, ampiezza 250mA aumentata di un fattore 10")

f3_pre = amp3_pre*np.cos(2*np.pi/0.001024*x+Phi3_pre) 
f3_post = amp3_post*np.cos(2*np.pi/0.001024*x+Phi3_post)


plt.figure(3)
plt.plot(x,f1_pre,"b")
plt.plot(x,f1_post*10,"r")
plt.grid(True)
plt.title("colonna 3, ampiezza 250mA aumentata di un fattore 10")



plt.show()













   

