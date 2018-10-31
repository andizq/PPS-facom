import sys
import numpy as np
import matplotlib.pyplot as plt
import math
import time
from scipy import optimize, integrate, interpolate
from scipy.interpolate import UnivariateSpline
from scipy.misc import derivative

###############################################################################################
#Constants
G = 6.674e-11 		#Cavendish ct
au = 1.496e11		#Astronomical Unit
mH = 1.672e-27		#Proton mass
kb = 1.38e-23 		#Boltzmann constant
yr = 86400*365		#yr in seconds
Mearth = 5.98e24        #kg Earth mass

###############################################################################################
#Sun properties
Msun = 1.989e30		#Solar mass
Lsun = 3.98e26		#Solar luminosity
Tsun = 5780.		#Solar Eff. Temperature
Rsun = 6.957e8		#Solar radius in m.
Zsun = 0.0134		#Solar Metallicity 
Xsun = 0.7381		#Solar Hydrogen mass fraction
Ysun = 0.2484		#Solar Helium Mass Fraction	
f_dgsun = Zsun/Xsun	#Dust-Gas ratio for Sun

###############################################################################################
#Disc properties and parameters
alpha=0.001 			#Viscosity coefficient
mp= 2.0*mH 			#Mass of the gas molecules
mu= 2.3 #mp/mH			#Mean particles mass
E_10 = 750. 			#Surface gas density of the disc at 10au in Kg m^-2
Es_10 = 5 * 3.2 / 1.51414 	#Surface dust density of the disc at 10au in Kg m^-2 (Only silicates!!!)
tau_dep= 2e6 * yr   		#Disc mean life-time: Replace for eq. (6) Ida&Lin-2008
tau_acc= 1.* tau_dep      	#Dust mean life-time
q=1.5 				#Power-law exponent in the surface density
qs=1.5                  	#Power-law exponent in the dust-surface density
T_ionH=1.e4			#Ionization temperature of Hydrogen in K
Beta_ionH=1.0           	#Loss-mass factor inside Rg 
PHI = 1.e41             	#photons s^-1 EUV Flux (Alexander-Armitage 2.1.3)
alpha_B = 2.6e-19       	#m3 s^-1 Recombination coefficient for atomic hydrogen (Alexander-Armitage Appendix) 
fg0 = 1.                	#Initial gas fraction compared with MMSN (Ida VI)
#fd = fg = 1.			#Values for MMNS from Ida&Lin 2005
C1_mig = 1.              	#Delay constant for migration type I
m_mean = 1e19           	#Typical Mass of planetecimals 
delta_ac = 10           	#In Hill radius

class dict2obj(object):
    def __init__(self,dic={}):self.__dict__.update(dic)
    def __add__(self,other):
        for attr in other.__dict__.keys():exec("self.%s=other.%s"%(attr,attr))
        return self
