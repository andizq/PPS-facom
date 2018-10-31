# This script produces synthetic planetary masses at different semimajor axis as a function of Mstar.
# Based mainly on Benz et.al (2014) and Ida & Lin (2008).
# Developed by FACom - SEAP group 
# Pablo Cuartas Restrepo (pablo.cuartas@udea.edu.co) - Andres Izquierdo (andres.izquierdo@udea.edu.co) - Jaime Alvarado (jaimea.alvarado@udea.edu.co)
# V1: 09-2016 ->
# Universidad de Antioquia - FCEN - Instituto de Fisica

from util import *

###############################################################
# 1. Vertical and Radial Structure
###############################################################

#Disc properties and parameters:

###################################################################
# General conditions for the disk
###################################################################

#Keplerian frecuency
def Omega(r,Mstar):
  Omega = np.sqrt(G*Mstar/r**3)
  return Omega

#Disk temperature profile (Armitage, AoPF - Hayashi)
def Td(r,Tstar,Rstar,Lstar):
  #Td = Tstar*(np.pi**-1*(np.arcsin(Rstar/r)-(Rstar/r)*(1-(Rstar/r)**2)**0.5))**0.
  Td = 280.*(r/au)**-0.5*(Lstar/Lsun)**0.25
  return Td

#Speed of sound (Armitage, AoPF)
def cs(T):
  return (kb * T / (mu*mH))**0.5

#Vertical disk scale-height (Armitage, AoPF)   
def h(R,omega,Tstar,Rstar,Lstar):
  Tdisc = Td(R,Tstar,Rstar,Lstar)
  cs_disc = cs(Tdisc)
  h = cs_disc / omega
  return h

#Ices line (Benz et.al, 2014)
def a_ice(Lstar):
  a_ice = 2.7*(Lstar/Lsun)**0.5*au
  return a_ice


###################################################################
# Disc Properties and scales
###################################################################

#Gravitational radius
def Rg(T_ionH,Mstar):
  Rg = 8.9*(Mstar/Msun)* (cs(T_ionH) / 1e4)**-2*au
  return Rg

#Initial gas fraction compared to MMSN (Ida 1)
def fg0(Mstar): 
  return Mstar/Msun

#Disk gas fraction (Benz et.al. 2014)
def fg(t,tau_dep,fg0):
  return fg0 * np.exp(-t / tau_dep)

#Dust mean life time: tau_acc minimization function for eqs 9 and 10 from Ida&Lin IV 
def tau_acc(r,t_obs, Mstar, fg0, fg, fd0, a_ice):

  eta_ice = 1
  expon0 = 27/10. + (qs-3/2.) + 2/5.*(q-3/2.)
  expon1 = 3 + (qs-3/2.)

  def fun0(x,r,t_obs): 
    f = (x - 2.2e5 * eta_ice**-1 * fd0**-1 * np.exp(t_obs / x) * 
         (fg0**-1 * np.exp(t_obs / tau_dep))** 0.4 * 10**(0.4 * (1.5 - q)) * 
         (r/au)**expon0 * (Mstar/Msun)**(-1/6.) * yr)
    return f 
  def fun1(x,r,t_obs):
    f = (x - 2e7 * eta_ice**-1 * fd0**-1 * np.exp(t_obs / x) * 
         (r/au)**expon1 * (Mstar/Msun)**(-1/2.) * yr)
    return f

  accrete, residuo = [], []
  i, tmp = 0, 0
  #print (t_obs)
  for t in t_obs:
    result=[]
    if fg[i] < 1e-3:
      fun = fun0
      x0 = 1e7
    else:
      fun = fun1
      x0 = 1e5
    for rval in r:
      #while(tmp<0):
      tmp = optimize.root(fun, x0 = x0 * yr, args = (rval,t) ).x 
      result.append(tmp)
      #residuo.append(fun(tmp,rval))
      accrete.append(result)
    i+= 1

  return np.array(result)

#Initial dust fraction (Ida 4, Benz 2014)
def fd0(Zstar,fg0):
  return fg0 * Zstar / Zsun 

#Disk dust fraction
def fd(t,tau_acc,fd0,flag = False):
  if flag:
    fd = [fd0 * np.exp(-t/tau_list) for tau_list in tau_acc]
    return np.array(fd)
  else: return fd0 * np.exp(-1*t/tau_acc) 

#Disk surface gas density
def Eg(r,Rg,fg,flag = False):
  if flag:
    Eg,CorteExp = [],[]
    for r_list in r:
      Eg.append( E_10 * fg * (r_list/(10 * au))**-q)  	#Ida&Lin 2008, eq. (1)
      CorteExp.append(np.exp(-r_list/Rg)) 
    Eg,CortExp = np.array(Eg),np.array(CorteExp)
    return Eg*CorteExp
  else:
    Eg = E_10 * fg * (r/(10 * au))**-q  
    CorteExp = np.exp(-r/Rg) 
    return Eg*CorteExp

#Disk surface solids density
def Es(r,a_ice,fd, flag = False):
  
  if flag:
    eta_ice = np.array([np.where( r_list < a_ice, 1, 4.2) for r_list in r])
    Es = []
    for r_list, eta_list, fd_list in list( (zip(r,eta_ice,fd)) ) :
      Es.append(Es_10 * eta_list * fd_list * (r_list/(10*au))**-qs)
    return np.array(Es) 

  else:
    if r < a_ice: eta_ice = 1.0
    else: eta_ice=4.2
    return Es_10 * eta_ice * fd * (r/(10*au))**-qs

#Embryo and core mass acretion
def Miso(r, Es, Mstar):  
  #term = 0.16 * (Es/100.) ** 1.5 * (r/au)**3 * (delta_ac/10)**1.5 * (Mstar/Msun)**-0.5 * Mearth 
  term = (8*np.pi*3**0.5 * Es * r**2)**1.5 / (3*Mstar)**0.5 
  #term_scat = 1.4e3 * (r/au)**-1.5 * (Mstar/Msun) * Mearth
  return term #term_scat

#Mass growth rate
def Mcore_dot(Mcore, tau_acc):
  return Mcore / tau_acc

#Hydro-mass function
def Mhydro(Mcore,tau_acc):
  return ( (Mcore / tau_acc) / ( 1e-6*Mearth/yr ))**0.25 * Mearth

#Core mass function
def Mcore(r, t, Es, Eg, Miso, Mstar):
  term1 = (t/(4.8e5*yr))**3 * (Es/100.)**3 * (Eg/(2.4e4))**(6/5.) 
  term2 = (m_mean/1e19)**-0.4 * (r/au)**(-9/5.) * (Mstar/Msun)**0.5 * Mearth 
  M = term1 * term2
  return M #np.where(M > Miso, Miso, M)

#----------
#MIGRATION
#----------
def tau_mig1(r, fg, Mcore, Mstar):
  tmig1 = 1.6e5 * C1_mig**-1 * fg**-1 * (Mcore / Mearth)**-1 * (r / au) * (Mstar / Msun)**1.5 * yr
  return tmig1

def da_dt(r, tau_mig1):
  return r/tau_mig1
  
"""
#Disk surface solids density
def EsH(r,t,tau_dep,qs,Lstar,Zstar):
  tauacc = tau_acc(r,t,tau_dep,qs,q,Mstar,Zstar)
  eta_ice=0
  if r<a_ice(Lstar): 
    eta_ice=1.
  else: 
    eta_ice=4.2
  #EsH = 100.*eta_ice*fd*(r/au)**-qs	#As Hayashi 1981, Ida&Lin 2004b eq (1)
  EsH = 100.*eta_ice*fd(r,t,tau_dep,qs,q,Mstar,Zstar)*(r/au)**-1.5
  return EsH
"""
"""

#Disk mean density
def rho(z,r,t,tau_dep,q,Tstar,Rstar,Mstar):
  h_val = h(r,Tstar,Rstar,Mstar)
  E_val = Eg(r,t,tau_dep,q,Mstar)
  rho = np.exp(-z*z/(2*h_val**2))*E_val/(h_val*(2*np.pi)**0.5)
  return rho

#Disk molecular viscosity
"""
"""
#Viscosity function only with parameters in the disk mean plane
def viscosity(dz,r,t,tau_dep,q,cs,omega,Tstar,Rstar,Mstar,flag=0):

  h_val = h(r,Tstar,Rstar,Mstar)
  E_val = Eg(r,t,tau_dep,q,Mstar)  
  nu0=alpha*cs**2/omega
  
  #Array to save the densities along z-axis (height-scale)
  rho_array=[]  
  for z in np.arange(-h_val,h_val+dz,dz):    
    rho_array.append(rho(z,r,t,tau_dep,q,Tstar,Rstar,Mstar))    
  
  Visco = integrate.trapz(nu0*np.array(rho_array),dx=dz) #Integrating nu0*rho
    
  if flag: return Visco
  else: return Visco/E_val
"""
"""

#Viscosity function integrated along the disk vertical scale
def viscosity(dz,r,t,tau_dep,q,Tstar,Rstar,Mstar,flag=0):
  
  #Array to save the densities along z-axis (height-scale)
  rho_array=[]
  nu0=[]
  E_val = Eg(r,t,tau_dep,q,Mstar)
  h_val = h(r,Tstar,Rstar,Mstar)
  
  for z in np.arange(-h_val,h_val+dz,dz):
    r=(r*r+z*z)**0.5
    cs_val = cs(r,Tstar,Rstar,Lstar)
    Omega_val = Omega(r,Mstar)  
    nu0.append(alpha*cs_val**2/Omega_val)
    rho_array.append(rho(z,r,t,tau_dep,q,Tstar,Rstar,Mstar))    
  
  Visco = integrate.trapz(np.array(nu0)*np.array(rho_array),dx=dz) #Integrating nu0*rho
  if flag: return Visco
  else: return Visco/E_val

def Mol_viscosity(r,t,tau_dep,q,cs,omega,Mstar):  
  #Array to save the densities along z-axis (height-scale)
  E_val = Eg(r,t,tau_dep,q,Mstar)
  nu = alpha * cs / omega
  return nu #/E_val

###############################################################
# 3. Other Disc properties, on Photo-evaporation 
###############################################################

#Stellar wind (energetic protons) density at Rg
def n0g(Rg): #Alexander-Armitage 2007 (n0(R14)-Benz))
  n0g1 = 0.14*(3*PHI/(4*np.pi*alpha_B*Rg**3))**0.5
  return n0g1

#Full density profile
def n0(r, Rg): #Benz
  r_Rg = r/Rg
  #n0 = n0g1*(r/(Beta_ionH*Rg))**-2.5 #Benz
  n0 = n0g(Rg) * (2./( r_Rg**7.5 + r_Rg**12.5))**0.2 #Alexander-Armitage Appendix
  return n0

#######################################################################
#Internal Photoevaporation due to EUV radiation
#######################################################################

def Edot_photoINT(r,PHI,alpha_B,Tstar,T_ionH,Rstar,Mstar, Rg, cs):   
  A = n0(r, Rg) * mH
  B = 2 * cs                 #Eq. 16 in Benz (2014) corrects Eq. (7) in Clark (2001)
  Edotphotoint = np.where(r < Rg*Beta_ionH, 0., A*B) 
  return Edotphotoint    

########################################################################
#Gas surface density evolution
########################################################################

def ViscoE(dz,r,t,tau_dep,q,Tstar,Rstar,Mstar):
  V = viscosity(dz,r,t,tau_dep,q,Tstar,Rstar,Mstar,flag=1)#*Eg(r,t,tau_dep,q)
  return V

########################################################################
# 4. Accretion rate
########################################################################

def tau_c_acc(r,t,tau_dep,qs,q,Mstar,Lstar,Zstar):
  aa_ice = a_ice(Lstar)
  if r < a_ice: eta_ice = 1.
  else: eta_ice = 4.2
  fdd = fd(r,t,tau_dep,qs,q,Mstar,Zstar)
  fgg = fg(t,tau_dep,Mstar)
  expon = 0
  if fgg >= 1e-3: 
    expon = 27/10. + (qs-1.5) + 0.4 * (q-1.5) 
    temp = (2.2e5 * eta_ice**-1 * fdd**-1 * fgg**-0.4 
            * 10**(0.4*(1.5-q)) * (r/au)**expon * (Mstar/Msun)**-(1/6.)) * yr 
  else: temp = 2e7 * eta_ice**-1 * fdd**-1 * (r/au)**(3+(qs-1.5)) * (Mstar/Msun)**-0.5 * yr
  return temp


def Mcore(r,t,tau_dep,qs,q,m_mean,Mstar,Lstar,Zstar):
  Ess = Es(r,t,tau_dep,qs,q,Mstar,Lstar,Zstar)
  Egg = Eg(r,t,tau_dep,q,Mstar)
  fdd = fd(r,t,tau_dep,qs,q,Mstar,Zstar)
  term1 = (t/(4.8e5*yr))**3 * (Ess/(10 * 10))**3 * (Egg/(10 * 2.4e3))**(6/5.) 
  term2 = (m_mean / 1e19)**-0.4 * (r/au)**(-9/5.) * (Mstar/Msun)**0.5 * Mearth 

def Mcore_GasDepletion(r,t,tau_dep,qs,q,m_mean,Mstar,Lstar,Zstar,fg0): #Ida y Lin II eq (3). Mcore when gas depletion 
  fdd = fd(r,t,tau_dep,qs,q,Mstar,Zstar)
  if r < a_ice: eta_ice = 1.
  else: eta_ice = 4.2
  term1 = (tau_dep/(4.8e5 * yr))**3 * eta_ice**3 * fdd**3 * fg0**(6/5.) * (m_mean/1e19)**-0.4
  term2 = (r/au)**(-8.1) * (Mstar/Msun)**0.5 * (1 - np.exp(-2*t/(5*tau_dep))**3) * Mearth
  
  return term1 * term2

def Miso(r,t,tau_dep,qs,q,m_mean,Mstar,Lstar,Zstar):
  Ess = Es(r,t,tau_dep,qs,q,Mstar,Lstar,Zstar)
  fdd = fd(r,t,tau_dep,qs,q,Mstar,Zstar)
  if r < a_ice: 
    eta_ice = 1.
  else: eta_ice = 4.2
  term = 0.16 * eta_ice ** 1.5 * fdd**1.5 * (r/au)**(3/4. - 1.5*(qs-1.5)) * (Mstar/Msun)**-0.5 * Mearth
  term_iso = 0.16 * (Ess/(10 * 10))**1.5 * (r/au)**3 * (Mstar/Msun)**-0.5 * Mearth 
  term_scat = 1.4e3 * (r/au)**-1.5 * (Mstar/Msun) * Mearth
  return term_iso # term_scat

"""
"""

###############################################################
# 5. Envelope
###############################################################



###############################################################
# 6. Atmosphere
###############################################################



###############################################################
# 7. Infalling
###############################################################



###############################################################
# 8. Core structure
###############################################################



###############################################################
# 9. Migration
###############################################################

#def tau_mig(r,t,tau_dep,C1_mig,Mstar,Mcore):

###############################################################
# 10. Planet-Planet interaction
###############################################################

"""
