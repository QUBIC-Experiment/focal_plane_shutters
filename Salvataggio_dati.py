import numpy as np
import matplotlib.pyplot as plt
import string 
from scipy import optimize
import scipy.signal as signal
from numpy import pi, r_
import csv


dPhi = np.arange(10.)
TdPhi = np.arange(64.)
Tdev = np.arange(64.)
n_iterazioni = np.arange(10)
n_switch = np.arange(64)
dPhi1 = np.arange(10.)
TdPhi1 = np.arange(64.)
Tdev1 = np.arange(64.)
amp = np.arange(10.)
amp1 = np.arange(10.)
Tdevamp1 = np.arange(64.)
amp2 = np.arange(10.)
Tdevamp2 = np.arange(64.)
Tdevamp = np.arange(64.)
Tamp = np.arange(64.)
Tamp1 = np.arange(64.)
Tamp2 = np.arange(64.)
Phi = np.arange(10.)
Phi1 = np.arange(10.)
Phi2 = np.arange(10.)
TPhi = np.arange(64.)
TPhi1 = np.arange(64.)
TPhi2= np.arange(64.)
Tdevp  = np.arange(64.)
Tdevp1 = np.arange(64.)
Tdevp2 =np.arange(64.)
matrice = np.zeros((5000.,64.))
lunghezza = np.arange(0,5000)

for l in n_switch:
    for i in n_iterazioni:
        if (l<9):
            path = "MisureQubic/TcryoVAsw63ma0020n09trA/TcryoVAsw0" + np.str(l+1) + "ma0020n0" + np.str(i) + "trA.txt"
        else:
            path ="MisureQubic/TcryoVAsw63ma0020n09trA/TcryoVAsw" + np.str(l+1) + "ma0020n0" + np.str(i) + "trA.txt"
        #print l, i, path
        data=np.loadtxt(path)
    
        startdata = 0
        enddata = 5000
        dt = 1e-6
        x = np.arange (0, dt*(enddata-startdata), dt)
        y = data [startdata:enddata,0]
        z = data [startdata:enddata,1]
        p = data [startdata:enddata,2]
        #Creo il filtro FIR
        nyq_rate = len(x)/2
        cut_off = 60.
        N = 100
        taps = signal.firwin (N, cut_off/nyq_rate)

        #Filtro x 
        fy = signal.lfilter (taps, 1.0, y)

    
        #Segnale originale
        #plt.figure(2)
        #plt.plot(x,y,"orange")

        #Trovo la variazione di fase del segnale originale
        shift = (0.5*(N-1))*dt


        #Plot del segnale filtrato spostato (I primi N-1 posti sono influenzati dalle condizioni iniziali)
        #plt.figure(2)
        #plt.plot(x[N-1:]-shift,fy[N-1:], "teal",linewidth=2)
        #plt.grid(True)
        #plt.xlabel("T (microsecondi)")
        #plt.title("segnale originale V con segnale filtrato")


        #Derivo la funzione filtrata
        der = np.diff(fy[N-1:])/dt
    

        #Fit col metodo dei minimi quadrati
        fitfunc = lambda p, x: p[0]*np.cos(2*np.pi/0.001024*x+p[1])  # Target function
        errfunc = lambda p, x, y: fitfunc(p, x[N:]) - der # Distance to the target function
        p0 = [5000., 0.] # Initial guess for the parameters
        p1, success = optimize.leastsq(errfunc, p0[:], args=(x, der), factor=1, maxfev = 100000)


        #print("Frequenza = %.2f Hz" %(1/0.001024))
        #print ("Phi = %.2f rad" %p1[1])


        #Filtro z 
        fz = signal.lfilter (taps, 1.0, z)

        #Segnale originale
        #plt.figure(5)
        #plt.plot(x,z,"orange")

        #Trovo la variazione di fase del segnale originale
        shift = (0.5*(N-1))*dt


        #Plot del segnale filtrato spostato (I primi N-1 posti sono influenzati dalle condizioni iniziali)
        #plt.figure(5)
        #plt.plot(x[N-1:]-shift,fz[N-1:], "teal",linewidth=2)
        #plt.grid(True)
        #plt.xlabel("T (microsecondi)")
        #plt.title("segnale originale I con segnale filtrato")


        #Derivo la funzione filtrata
        der = np.diff(fz[N-1:])/dt
    

        #Fit col metodo dei minimi quadrati
        fitfunc = lambda p, x: p[0]*np.cos(2*np.pi/0.001024*x+p[1])  # Target function
        errfunc = lambda p, x, y: fitfunc(p, x[N:]) - der # Distance to the target function
        p0 = [5000., 0.] # Initial guess for the parameters
        pz, success = optimize.leastsq(errfunc, p0[:], args=(x, der), factor=1, maxfev = 100000)


  
        #print("Frequenza = %.2f Hz" %(1/0.001024))
        #print ("Phi = %.2f rad" %pz[1])
        
        #Filtro p 
        fp = signal.lfilter (taps, 1.0, p)

       
        #Segnale originale
        plt.figure(8)
        plt.plot(x,p,"orange")

        #Trovo la variazione di fase del segnale originale
        shift = (0.5*(N-1))*dt


        #Plot del segnale filtrato spostato (I primi N-1 posti sono influenzati dalle condizioni iniziali)
        #plt.figure(9)
        #plt.plot(x[N-1:]-shift,fp[N-1:], "teal",linewidth=2)
        #plt.grid(True)
        #plt.xlabel("T (microsecondi)")
        #plt.title("segnale originale V con segnale filtrato")


        #Derivo e plotto la funzione filtrata
        der = np.diff(fp[N-1:])/dt
        #plt.figure(10)
        #plt.plot(x[N:],der,"")

        #Fit col metodo dei minimi quadrati
        fitfunc = lambda p, x: p[0]*np.cos(2*np.pi/0.001024*x+p[1])  # Target function
        errfunc = lambda p, x, y: fitfunc(p, x[N:]) - der # Distance to the target function
        p0 = [5000., 0.] # Initial guess for the parameters
        pp, success = optimize.leastsq(errfunc, p0[:], args=(x, der), factor=1, maxfev = 100000)


        #plt.plot(x, fitfunc(pp, x), "r-",linewidth=3) # Plot of the data and the fit
        #plt.xlabel("T (microsecondi)")
        #plt.title ("Derivata della funzione V con fit")
        #plt.grid(True)
        #print("Frequenza = %.2f Hz" %(1/0.001024))
        #print ("Phi = %.2f rad" %pp[1])

        #print  "differenza di fase tra la prima e la seconda colonna e': %.2f rad" %(p1[1]-pz[1])
        #print  "differenza di fase tra la prima e la terza colonna e': %.2f rad" %(p1[1]-pp[1])

        
        
        
        
        #print  "differenza di fase = %.2f rad" %(p1[1]-pz[1])
        dPhi[i]  = p1[1]-pz[1]
        dPhi1[i] = p1[1]-pp[1]
        amp[i]   = p1[0]
        amp1[i]  = pz[0]
        amp2[i]  = pp[0]
        Phi[i]   = p1[1]
        Phi1[i]  = pz[1]
        Phi2[i]  = pp[1]
        
        #plt.show()
    #print""
    #print "Vettore variazione di fase " 
    #print ""   
    #print Phi
    mean  =  np.mean(dPhi)
    dev   = np.std(Phi)
    mean1 = np.mean(dPhi1)
    dev1  = np.std(dPhi1)
    ph    = np.mean(Phi)
    ph1   = np.mean(Phi1) 
    ph2   = np.mean(Phi2)
    dph   = np.std(Phi)
    dph1  = np.std(Phi1) 
    dph2  = np.std(Phi2)
    a     = np.mean(amp)
    a1    = np.mean(amp1)
    a2    = (np.mean(amp2))
    da     = np.std(amp)
    da1    = np.std(amp1)
    da2    = np.std(amp2)
    
    
    
    #media e deviazione stansard della differenza tra prima terza colonna e tra prima seconda colonna
    Tdev1[l] = dev1
    TdPhi1[l] = mean1
    TdPhi[l]  = mean
    Tdev[l]   = dev
    #media delle fasi (tutte e tre le colonne)
    TPhi[l]   = ph
    TPhi1[l]  = ph1
    TPhi2[l]  = ph2
    #deviazioni standard delle fasi
    Tdevp[l]  = dph
    Tdevp1[l] = dph1
    Tdevp2[l] = dph2
    #ampiezze medie delle tre colonne
    Tamp[l]    =  a
    Tamp1[l]   = a1
    Tamp2[l]   = a2
    #deviazione standard ampiezze
    Tdevamp[l]= da
    Tdevamp1[l]=da1
    Tdevamp2[l]=da2
    
print "variazione di fase prima seconda colonna "
print TdPhi
print "variazione di fase prima terza colonna"
print TdPhi1
print "fase prima colonna"
print TPhi
print "fase seconda colonna"
print TPhi1
print "fase terza colonna"    
print TPhi2
print "ampiezza prima colonna"
print Tamp
print "ampiezza seconda colonna"
print Tamp1
print "ampiezza terza colonna"
print Tamp2

for i in n_switch:
    plt.plot(x,Tamp1[i]*np.cos(2*np.pi/0.001024*x+TPhi1[i]))
    
for i in n_switch:
    for o in lunghezza:
        matrice[o][i] =Tamp2[i]*np.cos(2*np.pi/0.001024*o+TPhi2[i])
    
#B = zip(h)
np.savetxt("MisureQubic/Plot3.txt",matrice, fmt = "%.5f", delimiter = " " )

A = zip(TdPhi,TdPhi1,TPhi,TPhi1,TPhi2,Tamp,Tamp1,Tamp2)
#print A
with open("MisureQubic/Tcryodata.txt","w") as f:
    writer = csv.writer(f, delimiter= " ")
    writer.writerows(A)
plt.show()   

quit()
    
