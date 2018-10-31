# This script produces synthetic planetary masses at different semimajor axis as a function of Mstar.
# Based mainly on Benz et.al (2014) and Ida & Lin (2008).
# Developed by FACom - SEAP group 
# Pablo Cuartas Restrepo (pablo.cuartas@udea.edu.co) - Andres Izquierdo (andres.izquierdo@udea.edu.co) - Jaime Alvarado (jaimea.alvarado@udea.edu.co)
# V1: 09-2016 ->
# Universidad de Antioquia - FCEN - Instituto de Fisica

from util import *
import discGen as dG
import sys

t0 = time.time()
#########################################################################
#Central star properties

StarDict = {'M': 1.*Msun, 'L': 1.*Lsun, 'T': 1.*Tsun , 'R': 1.*Rsun,
            'Z': 1.*Zsun, 'Y': 1.*Xsun, 'X': 1.*Xsun}

Star = dict2obj(StarDict)

#########################################################################
#Integration parameters

#DISC
amin = 0.1 * au        #Minimun disc radius
amax = 100 * au	       #Maximun disc radius
R_step = 0.1 * au / 10     #Integration step for "r" in Disc
R = np.arange(amin,amax,R_step)

#########################################################################
#Orbital frequency
omega = dG.Omega(R,Star.M)

#########################################################################
#Scaleheight
dz = 1e-2 * au
h_val = dG.h(R,omega, Star.T, Star.R, Star.L) 
zList = np.array([np.linspace(-h, h , 2*h / dz) 
                  if int(2*h / dz) % 2 == 1 else np.linspace(-h, h , 2*h / dz + 1) 
                  for h in h_val]) #Forcing the number of points to be odd (then the midplane will always exist)
#zList = np.where( np.array( map( lambda x: len(x), zList)) == 1, [np.array([0.])], zList)[0]
args_z = np.where( np.array(list(map( lambda x: len(x), zList))) == 1)
zList[args_z] = [np.array([0.])]
eps = 1e-3
args_z0 = [np.where(abs(zList[i]) < eps)[0] for i in range(len(zList))] #Indices for plane z=0

#########################################################################
#r coordinate for each point
it = iter(R**2)
r = np.array([ np.sqrt( next(it) + z**2 ) for z in zList ])

#########################################################################
#Temp. profile
T = dG.Td(r,Star.T,Star.R,Star.L)

#########################################################################
#Speed of sound 
cs = dG.cs(T)

#########################################################################
#TIME
t_step = 1e4 * yr             #Integration step for "t" in Time
t_max = 1 * tau_dep / 2 
t_min = 0.01 * yr
Time = np.arange(t_min,t_max,t_step)

#########################################################################
#Gravitational Radius
Rg = dG.Rg(T_ionH,Star.M)
fd0 = dG.fd0(Star.Z, fg0)
a_ice = dG.a_ice(Star.L)

#########################################################################
#Accretion mean time
#t_obs = 1e7 * yr
t_obs = Time[-1] / 10 #tau_dep
r_obs = 1 * au
#arg = np.where(Time<=t_obs)[0][-1]

fix_time = True
fix_pos = False

#in R
if fix_time: 
    print ("Calculating tau_acc for a fixed time")
    tt = np.array([t_obs])
    fg = dG.fg(tt, tau_dep, fg0)[0]
    tau_acc = dG.tau_acc(R, tt, Star.M, fg0, np.array([fg]), fd0, a_ice)
    fd = dG.fd(t_obs,tau_acc,fd0,flag=True).T[0]
    Eg = dG.Eg(R,Rg,fg,flag=True) 
    Es = dG.Es(R,a_ice,fd,flag=True)
    Miso0 = dG.Miso(R,Es,Star.M)
    Mcore0 = dG.Mcore(R,t_obs,Es,Eg,Miso0,Star.M) + m_mean
    Mhydro0 = dG.Mhydro(Mcore0,tau_acc.T)[0]
elif fix_pos: 
    print ("Calculating tau_acc for a fixed position")
    fg = dG.fg(Time, tau_dep, fg0)
    tau_acc = dG.tau_acc([r_obs], Time, Star.M, fg0, fg, fd0, a_ice)
    fd = dG.fd(Time,tau_acc,fd0)[0]
    Eg = dG.Eg(r_obs,Rg,fg)
    Es = dG.Es(r_obs,a_ice,fd)
    Miso0 = dG.Miso(r_obs,Es,Star.M)
    Mcore0 = dG.Mcore(r_obs,Time,Es,Eg,Miso0,Star.M) + m_mean
    Mhydro0 = dG.Mhydro(Mcore0,tau_acc.T)[0]
else: sys.exit("Please set either a fixed time or a fixed position to compute the model")

"""
#in r
if fix_time: 
    print ("Calculating tau_acc for a fixed time")
    tau_acc = dG.tau_acc(r, t_obs, Star.M, fg0, fg[arg], fd0, a_ice)
    fd = dG.fd(t_obs,tau_acc,fd0,flag=True)
    Eg = dG.Eg(r,Rg,fg[arg],flag=True)
    Es = dG.Es(r,a_ice,fd,flag=True)
    Eg_z0 = np.array([Eg[i] [args_z0[i]] [0] for i in range(len(args_z0))])
    fd_z0 = np.array([fd[i] [args_z0[i]] [0] for i in range(len(args_z0))])
    Es_z0 = np.array([Es[i] [args_z0[i]] [0] for i in range(len(args_z0))])
elif fix_pos: 
    print ("Calculating tau_acc for a fixed position")
    tau_acc = dG.tau_acc(r_obs, Time, Star.M, fg0, fg, fd0, a_ice)
    fd = dG.fd(Time,tau_acc,fd0)
    Eg = dG.Eg(r_obs,Rg,fg)
else: sys.exit("Please set either a fixed time or a fixed position to compute the model")
"""

plt.title('$\Sigma_s$')
plt.plot(R/au,Es/10, label = "t = %.0e yr"%(t_obs/yr)) 
plt.xscale('log')
plt.yscale('log')
plt.xlim(R[0]/au,100)
plt.ylim(1e-1,1e3)
plt.legend(loc='best', fontsize = 18)
plt.savefig("Sigma_s.png")
plt.close()

plt.title('$\Sigma_g$')
plt.plot(R/au,Eg/10, label = "t = %.0e yr"%(t_obs/yr)) 
plt.xscale('log')
plt.yscale('log')
plt.xlim(R[0]/au,100)
plt.ylim(1e-1,1e5)
plt.legend(loc='best', fontsize = 18)
plt.savefig("Sigma_g.png")
plt.close()

#########################################################################
#FULL CALCULATION ALONG R and t
#########################################################################

#########################################################################
#SOME Initializations
#########################################################################

num_planets = 40 
a_lims = np.log10(5*amin/au), np.log10(0.3*amax/au) #Min and Max possible a's for random planets
r_obs = np.sort( 10**np.random.uniform(a_lims[0], a_lims[1], size = num_planets) ) * au #Random positions
r_obs_n = np.zeros( (num_planets, len(Time)) ) #Matrix where the time evolution of 'a' will be stored

fd, Eg, Es, tau_acc, tau_mig1, da_mig = [[0 for _ in xrange(num_planets)] for __ in xrange(6)] 
r_list, Mgas_list, Mhydro_list, Miso, Mplanet_final,  Mg0 = [[0 for _ in xrange(num_planets)] for __ in xrange(6)]  
flag_iso, flag_hydro = [[0 for _ in xrange(num_planets)] for __ in xrange(2)]  

Mcore, Mplanet, Mplanet_list, Mgas, Mhydro = [np.zeros( (num_planets, len(Time)) ) for _ in xrange(5)]
#########################################################################
#########################################################################

i_ref, i, j, r_now = [0, 1, 0, 0]
K = 1e9 * Mearth**3 * yr

fg = dG.fg(Time, tau_dep, fg0)
for a_obs in r_obs:
    r_now = a_obs
    tau_acc = dG.tau_acc([r_now], [Time[0]], Star.M, fg0, [fg[i]], fd0, a_ice)
    fd = dG.fd(Time[0],tau_acc,fd0)
    Eg = dG.Eg(r_now,Rg,fg[0])
    Es = dG.Es(r_now,a_ice,fd)
    Miso = dG.Miso(r_now,Es,Star.M)
    Mcore0 = dG.Mcore(r_now,Time[0],Es,Eg,Miso,Star.M) + m_mean
    tau_mig1 = dG.tau_mig1(r_now,fg[0],Mcore0,Star.M)
    da_mig = dG.da_dt(r_now, tau_mig1)
    Mhydro = dG.Mhydro(Mcore0,tau_acc)
    Mcore_dot = dG.Mcore_dot(Mcore0,tau_acc)
    Mplanet = ( 3 * K * Mcore_dot )**0.25
    
    r_obs_n[j][0] = r_now
    Mcore[j][0] = Mcore0
    Mplanet_list[j][0] = Mcore[j][0]
    
    for t in Time[1:]:        
        r_obs_n[j][i] = r_now - da_mig*t_step
        #r_now = r_obs_n[j][i]
        #Mcore_j -> now
        if r_obs_n[j][i] < 0:
            r_obs_n[j][i] = 0.1             
        r_now = r_obs_n[j][i]
        tau_acc = dG.tau_acc([r_now], [t], Star.M, fg0, [fg[i]], fd0, a_ice)
        fd = dG.fd(t,tau_acc,fd0)
        Eg = dG.Eg(r_now,Rg,fg[i])
        Es = dG.Es(r_now,a_ice,fd)
        Miso = dG.Miso(r_now,Es,Star.M)

        if flag_iso[j]: Mcore[j][i] = Mcore[j][i_ref]
        else: Mcore[j][i] = dG.Mcore(r_now,t,Es,Eg,Miso,Star.M) + m_mean
        #Mcore[j][i] = dG.Mcore(r_now,t,Es,Eg,Miso,Star.M) + m_mean
        tau_mig1 = dG.tau_mig1(r_now,fg[i],Mcore[j][i],Star.M)
        da_mig = dG.da_dt(r_now , tau_mig1)
        
        Mhydro = dG.Mhydro(Mcore[j][i],tau_acc)
        Mcore_dot = dG.Mcore_dot(Mcore[j][i],tau_acc)
        
        try:
            Mplanet = ( 3 * K * Mcore_dot )**0.25
        except RuntimeWarning:
            print(tau_acc)
        
        if ( Mcore[j][i] >= Miso or ( r_now < a_ice and a_obs > a_ice ) 
             and flag_iso[j] == False ):

            Mcore[j][i] = Miso
            i_ref = i
            flag_iso[j] = True
        
        if Mcore[j][i] >= Mhydro and flag_hydro[j] == False:
            Mg0[j] = 2 * Mcore_dot * tau_dep
            flag_hydro[j] = True
        
        if flag_hydro[j]:
            Mgas = ( 2 * Mcore_dot * tau_dep) - Mg0[j]
            if Mgas < 0: Mgas = 0
        else: Mgas = 0

        Mgas_list[j] = Mgas
        Mplanet_list[j][i] = Mcore[j][i] + Mgas
        Mplanet_final[j] = Mplanet_list[j][i] 
        r_list[j] = r_obs_n[j][i-1]
            
        i += 1
    i = 0
    j += 1

#010618: Arreglar el ultimo paso en t --> ver plot migracion    
    
Mgas_list = np.array(Mgas_list).T
Mplanet_final = np.array(Mplanet_final).T

plt.figure(figsize=(10,5))
Tf_str = '%1.e'%(Time[-1]/yr)
plt.title(r"Type $\mathrm{I}$ Migration, 10$^{%d}$ $\rightarrow$ %s$\times$10$^{%s}$ yr"%tuple([np.log10(Time[0]/yr), Tf_str[0], Tf_str[-1]]))
plt.plot(R/au, Miso0/Mearth, 'r--', label = r'M$_{iso},$ $10^%d$ $yr$'%round(np.log10(t_obs/yr)))
plt.plot(R/au, Mhydro0/Mearth, 'b--', label = r'M$_{hydro},$ $10^%d$ $yr$'%round(np.log10(t_obs/yr)))
plt.yscale('log')
plt.xlabel('R (au)', fontsize=16)
plt.ylabel('Embryo Mass ($M_{\oplus}$)', fontsize=18)
plt.xlim(0.5,12)
plt.ylim(0.001,20)

for j in range(num_planets): plt.plot(r_obs_n[j] / au, Mplanet_list[j] / Mearth, 'o-', color='g', ms = 2)

rplot = np.array(r_list)/au
#plt.plot(rplot, np.array(Miso_list) / Mearth,'r*', ms=5.0, label='M$_{iso}$')
#plt.plot(rplot, np.array(Mhydro_list) / Mearth,'g*', ms=5.0, label='M$_{hydro}$')
#plt.plot(rplot, np.array(Mgas_list) / Mearth,'b*', ms=5.0, label='M$_{gas}$')
plt.plot(rplot, np.array(Mplanet_final) / Mearth,'y*', ms=5.0, label='M$_{planet}$')
plt.legend(fontsize = 10, loc = 'upper right')#, bbox_to_anchor=(0.95, 1.05))
plt.savefig("Migration_type_I.png")
plt.close()

#------------------------------
#PLOTS FOR POSTER, ORACLE 2018
#------------------------------
t_list = np.logspace(5, 7, 3) * yr

Eg, Es = [[0 for _ in xrange(len(t_list))] for __ in xrange(2)]

x_lab = 'R (au)'
y_lab = [r'$\Sigma_g$ (g cm$^{-2}$)', r'$\Sigma_d$ (g cm$^{-2}$)'] 

i = 0
fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(8, 8.))
for t_obs in t_list:
    fg = dG.fg(t_obs, tau_dep, fg0)
    tau_acc = dG.tau_acc(R, np.array([t_obs]), Star.M, fg0, np.array([fg]), fd0, a_ice)
    fd = dG.fd(t_obs,tau_acc,fd0,flag=True).T[0]
    Eg[i] = dG.Eg(R,Rg,fg,flag=True) 
    Es[i] = dG.Es(R,a_ice,fd,flag=True)    
    ax[0].plot(R/au, Eg[i]/10, label = r'10$^{%d}$ yr'%np.log10(t_obs/yr))    
    ax[1].plot(R/au, Es[i]/10, label = r'10$^{%d}$ yr'%np.log10(t_obs/yr))    
    i+=1

ax[0].set_xscale('log')
ax[0].set_yscale('log')
ax[0].set_ylabel(y_lab[0])
#ax[0].set_xlim(1,)

ax[1].set_xscale('log')
ax[1].set_yscale('log')
ax[1].set_xlabel(x_lab)
ax[1].set_ylabel(y_lab[1])

ax[0].legend(loc = "lower left")
ax[1].legend()
plt.savefig("Surface_densities.png")

print ("Time ellapsed",time.time() - t0)        


#PENDIENTES PARA 01-12-17:
# 1. Iniciar acrecion de gas... Formacion de gigantes.
# 2. Incluir migracion tipo II.
