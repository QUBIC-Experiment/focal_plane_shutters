import numpy as np
from numpy import pi, r_
import matplotlib.pyplot as plt
from scipy import optimize
import scipy.signal as signal


n_iterazioni = np.arange(10)



for i in n_iterazioni :
    data=np.loadtxt( "dati shutters 18042018" + "/test_18042018_01VAsw13" + "ma0250n0" + np.str(i) + "trA.txt")

    startdata = 00
    enddata = 8000
    dt = 1e-6
    x = np.arange (0, dt*(enddata-startdata), dt)
    y = data [startdata:enddata+1,0]
    z = data [startdata:enddata+1,1]
    p = data [startdata:enddata+1,2]

    #print len(x), len(y)
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

    plt.figure(1)
    plt.plot(x[N:],der)
    plt.plot(x, fitfunc(p1, x),linewidth=3) # Plot of the data and the fit
    
    print("Frequenza = %.2f Hz" %(1/0.001024))
    print ("Phi = %.2f rad" %p1[1])


    #Filtro z 
    fz = signal.lfilter (taps, 1.0, z)

       #Trovo la variazione di fase del segnale originale
    shift = (0.5*(N-1))*dt


    #Derivo e plotto la funzione filtrata
    der1 = np.diff(fz[N-1:])/dt
    
    #Fit col metodo dei minimi quadrati
    fitfunc = lambda p, x: p[0]*np.cos(2*np.pi/0.001024*x+p[1])  # Target function
    errfunc = lambda p, x, y: fitfunc(p, x[N:]) - der1 # Distance to the target function
    p0 = [5000., 0.] # Initial guess for the parameters
    pz, success = optimize.leastsq(errfunc, p0[:], args=(x, der1), factor=1, maxfev = 100000)
    
    plt.figure(2)
    plt.plot(x[N:],der1)
    plt.plot(x, fitfunc(pz, x),linewidth=3) # Plot of the data and the fit
    plt.xlabel("T (microsecondi)")
    plt.title ("Derivata della funzione I con fit")
    plt.grid(True)
    print("Frequenza = %.2f Hz" %(1/0.001024))
    print ("Phi = %.2f rad" %pz[1])


    #Filtro p 
    fp = signal.lfilter (taps, 1.0, p)

   
    #Trovo la variazione di fase del segnale originale
    shift = (0.5*(N-1))*dt

    #Derivo e plotto la funzione filtrata
    der2 = np.diff(fp[N-1:])/dt

    

    #Fit col metodo dei minimi quadrati
    fitfunc = lambda p, x: p[0]*np.cos(2*np.pi/0.001024*x+p[1])  # Target function
    errfunc = lambda p, x, y: fitfunc(p, x[N:]) - der2 # Distance to the target function
    p0 = [5000., 0.] # Initial guess for the parameters
    pp, success = optimize.leastsq(errfunc, p0[:], args=(x, der2), factor=1, maxfev = 100000)

    plt.figure(3)
    plt.plot(x, fitfunc(pp, x),"y",linewidth=3) # Plot of the data and the fit
    plt.xlabel("T (microsecondi)")
    plt.title ("Derivata della funzione V con fit")
    plt.grid(True)
    #print len(x[N:]),len(der)
    plt.plot(x[N:],der2)
    print("Frequenza = %.2f Hz" %(1/0.001024))
    print ("Phi = %.2f rad" %pp[1])





    print  "differenza di fase tra la prima e la seconda colonna e': %.2f rad" %(p1[1]-pz[1])
    print  "differenza di fase tra la prima e la terza colonna e': %.2f rad" %(p1[1]-pp[1])




plt.show()

