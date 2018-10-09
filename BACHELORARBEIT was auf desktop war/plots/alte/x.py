import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import scipy as sp
import csv
import matplotlib.ticker as plticker
from uncertainties import ufloat

x=np.linspace(-5,5,2000) #wie arange, nur anzahl der intervalle wird festgelegt (1001), nicht intervallgroesse

omega=2 #frequenz der treibekraft S,S+,S-, T=2pi/omega
F=1 #amplitude der treibkraft vom originalen x1
S_p=F/np.sqrt(2)
S_m=-F/np.sqrt(2)
def treibkraft(S,omega,t):  #treibkraft, ohne phase
	return S*np.sin(omega*t)
k=1
k_2=2
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
#plt.plot(x, L_P(x), "brown", label=("hi")


#plot zeta+
plt.plot(x, zeta_p(x), "royalblue", label="$zeta_+$")
plt.xlabel("$t/s$")
plt.ylabel("Auslenkung $\zeta /m$")
plt.title("Klassische Loesung zeta plus")
plt.grid()
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("zeta_+")
#plt.show()
#file.close()

#plot zeta- 
plt.plot(x, zeta_m(x), "blue", label="$zeta_-$")
plt.xlabel("$t/s$")
plt.ylabel("Auslenkung $\zeta /m$")
plt.title("Klassische Loesung zeta minus")
plt.grid()
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("zeta_-")
#plt.show()
#file.close()




#plot <x1>00
plt.plot(x,(1/np.sqrt(2))*(zeta_p(x) - zeta_m(x)), "red", label="<x1>")
plt.xlabel("$t/s$")
plt.ylabel("Auslenkung $x1 /m$")
plt.title("Erwartungswert <x1> fuer grundzustand $\psi_{00}$")
plt.grid()
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("<x1>00")

#plot <x2>00
plt.plot(x,(1/np.sqrt(2))*(zeta_p(x) + zeta_m(x)), "maroon", label="<x2>")
plt.xlabel("$t/s$")
plt.ylabel("Auslenkung $x2 /m$")
plt.title("Erwartungswert <x2> fuer grundzustand $\psi_{00}$")
plt.grid()
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("<x2>00")


#plot <H>00
#plt.plot(x, E_p_null+E_m_null -L_p(x)-L_m(x) +m*(abl_zeta_p(x))**2+m*(abl_zeta_m(x))**2, "black", label="<H>")
#plt.xlabel("$t/s$")
#plt.ylabel("<H>")
#plt.title("Erwartungswert <H> fuer grundzustand $\psi_{00}$")
#plt.grid()
#plt.legend(loc="best")
#plt.tight_layout()
#plt.savefig("<H>00")









	
