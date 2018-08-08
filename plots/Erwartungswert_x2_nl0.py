import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import scipy as sp
import csv
import matplotlib.ticker as plticker
from uncertainties import ufloat

x=np.linspace(-6,6,2000) #wie arange, nur anzahl der intervalle wird festgelegt (1001), nicht intervallgroesse

omega=2 #frequenz der treibekraft S,S+,S-, T=2pi/omega
F=1 #amplitude der treibkraft vom originalen x1
S_p=F/np.sqrt(2)
S_m=-F/np.sqrt(2)
def treibkraft(S,omega,t):  #treibkraft, ohne phase
	return S*np.sin(omega*t)
k=1
k_2=0   #3 x_1 ist 0, 0 x_2 ist 0
m=1
omega_p=np.sqrt(k/m) #omega+
omega_m=np.sqrt((k+2*k_2)/m) #omega-

#klassische Lsg und zeitliche Ableitung
def klass_lsg(omega_0,m,S,omega,t):
	return (treibkraft(S,omega,t))/(m*(omega_0**2 - omega**2)) 
#omega0 ist omega+- dann fuer zeta +-, S ist fuer S+-

def abl_klass_lsg(omega_0,m,S,omega,t):
	return (S*omega*np.cos(omega*t))/(m*(omega_0**2 - omega**2))
#ableitung treibkraft so hingeschrieben



#zeta+-
def zeta_p(t):
	return klass_lsg(omega_p,m,S_p,omega,t)
def zeta_m(t):
	return klass_lsg(omega_m,m,S_m,omega,t)

#zeta+-_punkt (zeitliche ableitung)
def abl_zeta_p(t):
	return abl_klass_lsg(omega_p,m,S_p,omega,t)
def abl_zeta_m(t):
	return abl_klass_lsg(omega_m,m,S_m,omega,t)




#ungetriebener oszi energien
h_quer=1 #natuerliche einheiten
E_p_null=(1/2)*omega_p*h_quer
E_m_null=(1/2)*omega_m*h_quer

#lagrange funktionen
def L_p(t):
	return (1/2)*m*(abl_zeta_p(t))**2 - (1/2)*m*omega_p**2*(
zeta_p(t))**2 + zeta_p(t)*treibkraft(S_p,omega,t)

def L_m(t):
	return (1/2)*m*(abl_zeta_m(t))**2 - (1/2)*m*omega_m**2*(
zeta_m(t))**2 + zeta_m(t)*treibkraft(S_m,omega,t)




plt.rcParams.update({'font.size': 20}) #allgemeine font size, hier nur achsen zahlen, rest nochmal geandert
#plt.rc('legend', fontsize=100) soll legend aendert klappt nicht
#plt.legend(loc=2, prop={'size': 100})


#plot <x12>00
plt.plot(x,(1/np.sqrt(2))*(zeta_p(x) - zeta_m(x)), "tomato", label="$<x_1>$", linewidth=5)
plt.plot(x,(1/np.sqrt(2))*(zeta_p(x) + zeta_m(x)), "maroon", label="$<x_2>$",linewidth=5)
plt.xlabel("$t/s$", fontsize=26)
plt.ylabel("$Auslenkung/m$", fontsize=26)
#plt.title("Erwartungswert $<x_1>_{00}(t)$,$<x_2>_{00}(t)$ \n fuer den Grundzustand $\psi_{00}$")
plt.grid()
plt.legend(loc="best", fontsize=23)
plt.tight_layout()
plt.savefig("<x2>nl0")
