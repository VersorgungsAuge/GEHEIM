import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import scipy as sp
import csv
import matplotlib.ticker as plticker
from uncertainties import ufloat

x=np.linspace(-5,5,2000) #wie arange, nur anzahl der intervalle wird festgelegt (1001), nicht intervallgroesse

omega=4 #frequenz der treibekraft S,S+,S-, T=2pi/omega
F=12 #amplitude der treibkraft vom originalen x1
S_p=F/np.sqrt(2)
S_m=-F/np.sqrt(2)
def treibkraft(S,omega,t):  #treibkraft, ohne phase
	return S*np.sin(omega*t)
k=13
k_2=5
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



#<H>00 durchschnitt
H_dach_nullnull = E_p_null+E_m_null - (1/8*m)*F**2*(1/(omega_p**2-omega**2) + 1/(omega_m**2-omega**2) - 2*omega**2/(omega_p**2-omega**2)**2 - 2*omega**2/(omega_m**2-omega**2)**2)

print(H_dach_nullnull)

#plot <H>00
plt.plot(x, E_p_null+E_m_null -L_p(x)-L_m(x) +m*(abl_zeta_p(x))**2+m*(abl_zeta_m(x))**2, "lime", label="$<H>_{00}$")
plt.axhline(y=H_dach_nullnull, color='m', linestyle='-', label="$H_{00,mittel}$")
plt.xlabel("$t/s$")
plt.ylabel("$Energie/J$")
plt.title("Erwartungswert $<H>_{00}$ fuer Grundzustand $\psi_{00}$  \n und ueber T gemittelter Erwartungswert $H_{00,mittel}$")
plt.grid()
plt.legend(loc="best")
plt.tight_layout()
#plt.show()
plt.savefig("<H>00")