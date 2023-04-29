from matplotlib import pyplot as plt
from matplotlib import rcParams
import numpy as np

rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
rcParams['errorbar.capsize'] = 3 
rcParams['lines.markersize'] = 10
rcParams['figure.figsize'] = 8, 5

kb=1.381e-23
e=1.602e-19

###########################################################################
#PLOTS OF GAS WITH VARYING TEMPERATURE
print("VARYING TEMPERATURE")

temperature_list=[500,525,550,575,600,625,650,675]
lattice_size_list=[50,30,30,30,50,30,30,30]
starting_pressure=1.0
starting_CO_fraction=0.5
 
plt.figure(figsize=(12,5))
plt.subplots_adjust(wspace =0.3)
print("PLOTS PER SITE")

#H2
plt.subplot(1,2,1)
for index in range(len(temperature_list)):
    temperature=temperature_list[index]
    lattice_sites=lattice_size_list[index]**2
    time, time_uncertainty, H2_amount, H2_amount_uncertainty = np.loadtxt(f"H2 Amounts (T={temperature}K,P={np.round(starting_pressure,2)}bar,x(CO)={starting_CO_fraction}).txt",unpack=True)
    H2_per_site=np.array(H2_amount)/lattice_sites
    H2_per_site_uncertainty=np.array(H2_amount_uncertainty)/lattice_sites
    plt.tick_params(axis='both', right=True, top=True, labelsize=15)
    plt.xlabel('Time Elapsed / s', fontsize=15)
    plt.ylabel('Number of H₂ Gas Molecules Per Site', fontsize=15)      
    plt.errorbar(time,H2_per_site,xerr=time_uncertainty,yerr=H2_per_site_uncertainty,fmt='-',label=f"{temperature} K")
plt.axvline(x=0,color='black',linestyle='--')
plt.axhline(y=0,color='black',linestyle='--')
plt.xlim(-0.5,10.5)
#plt.legend(title="Temperature:",fontsize=15,title_fontsize=15)
#plt.show()
  
#CO2
plt.subplot(1,2,2)
for index in range(len(temperature_list)):
    temperature=temperature_list[index]
    lattice_sites=lattice_size_list[index]**2
    time, time_uncertainty, CO2_amount, CO2_amount_uncertainty = np.loadtxt(f"CO2 Amounts (T={temperature}K,P={np.round(starting_pressure,2)}bar,x(CO)={starting_CO_fraction}).txt",unpack=True)
    CO2_per_site=np.array(CO2_amount)/lattice_sites
    CO2_per_site_uncertainty=np.array(CO2_amount_uncertainty)/lattice_sites
    plt.tick_params(axis='both', right=True, top=True, labelsize=15)
    plt.xlabel('Time Elapsed / s', fontsize=15)
    plt.ylabel('Number of CO₂ Gas Molecules Per Site', fontsize=15)      
    plt.errorbar(time,CO2_per_site,xerr=time_uncertainty,yerr=CO2_per_site_uncertainty,fmt='-',label=f"{temperature} K")
plt.axvline(x=0,color='black',linestyle='--')
plt.axhline(y=0,color='black',linestyle='--')
plt.xlim(-0.5,10.5)
plt.legend(title="Temperature:",loc="right",bbox_to_anchor =(1.45, 0.5), ncol = 1, fontsize=15,title_fontsize=15)
plt.show()


print("\nPLOTS PER SITE PER SECOND")
plt.figure(figsize=(12,5))
plt.subplots_adjust(wspace =0.25)

#H2
plt.subplot(1,2,1)
for index in range(len(temperature_list)):
    temperature=temperature_list[index]
    lattice_sites=lattice_size_list[index]**2
    time, time_uncertainty, H2_amount, H2_amount_uncertainty = np.loadtxt(f"H2 Amounts (T={temperature}K,P={np.round(starting_pressure,2)}bar,x(CO)={starting_CO_fraction}).txt",unpack=True)
    H2_per_site_per_time=np.array(H2_amount[1:])/(lattice_sites*np.array(time[1:]))
    H2_per_site_per_time_uncertainty=(1/lattice_sites) * ((np.array(H2_amount_uncertainty[1:]) / np.array(time[1:])) + ( (np.array(H2_amount_uncertainty[1:])*np.array(time_uncertainty[1:])) / (np.array(time[1:]))**2 ))
    plt.tick_params(axis='both', right=True, top=True, labelsize=15)
    plt.xlabel('Time Elapsed', fontsize=15)
    plt.ylabel('H₂ Gas Molecules Per Site Per Second', fontsize=15)      
    plt.errorbar(time[1:],H2_per_site_per_time,xerr=time_uncertainty[1:],yerr=H2_per_site_per_time_uncertainty,fmt='-',label=f"{temperature} K")
#plt.legend(title="Temperature:",fontsize=15,title_fontsize=15)
plt.axvline(x=0,color='black',linestyle='--')
plt.axhline(y=0,color='black',linestyle='--')
plt.xlim(-0.5,10.5)
#plt.ylim(-0.01,0.04)
#plt.show()
  
plt.subplot(1,2,2)
for index in range(len(temperature_list)):
    temperature=temperature_list[index]
    lattice_sites=lattice_size_list[index]**2
    time, time_uncertainty, CO2_amount, CO2_amount_uncertainty = np.loadtxt(f"CO2 Amounts (T={temperature}K,P={np.round(starting_pressure,2)}bar,x(CO)={starting_CO_fraction}).txt",unpack=True)
    CO2_per_site_per_time=np.array(CO2_amount[1:])/(lattice_sites*np.array(time[1:]))
    CO2_per_site_per_time_uncertainty=(1/lattice_sites) * ((np.array(CO2_amount_uncertainty[1:]) / np.array(time[1:])) + ( (np.array(CO2_amount_uncertainty[1:])*np.array(time_uncertainty[1:])) / (np.array(time[1:]))**2 ))
    plt.tick_params(axis='both', right=True, top=True, labelsize=15)
    plt.xlabel('Time Elapsed', fontsize=15)
    plt.ylabel('CO₂ Gas Molecules Per Site Per Second', fontsize=15)      
    plt.errorbar(time[1:],CO2_per_site_per_time,xerr=time_uncertainty[1:],yerr=CO2_per_site_per_time_uncertainty,fmt='-',label=f"{temperature} K")
plt.legend(title="Temperature:",loc="right",bbox_to_anchor =(1.45, 0.5), ncol = 1, fontsize=15,title_fontsize=15)
plt.axvline(x=0,color='black',linestyle='--')
plt.axhline(y=0,color='black',linestyle='--')
plt.xlim(-0.5,10.5)
plt.show()


print("\nPLOTS PER SITE PER SECOND FOR LOWER TEMPERATURES")
plt.figure(figsize=(12,5))
plt.subplots_adjust(wspace =0.25)

#H2
plt.subplot(1,2,1)
for index in range(len(temperature_list[:5])):
    temperature=temperature_list[index]
    lattice_sites=lattice_size_list[index]**2
    time, time_uncertainty, H2_amount, H2_amount_uncertainty = np.loadtxt(f"H2 Amounts (T={temperature}K,P={np.round(starting_pressure,2)}bar,x(CO)={starting_CO_fraction}).txt",unpack=True)
    H2_per_site_per_time=np.array(H2_amount[1:])/(lattice_sites*np.array(time[1:]))
    H2_per_site_per_time_uncertainty=(1/lattice_sites) * ((np.array(H2_amount_uncertainty[1:]) / np.array(time[1:])) + ( (np.array(H2_amount_uncertainty[1:])*np.array(time_uncertainty[1:])) / (np.array(time[1:]))**2 ))
    plt.tick_params(axis='both', right=True, top=True, labelsize=15)
    plt.xlabel('Time Elapsed', fontsize=15)
    plt.ylabel('H₂ Gas Molecules Per Site Per Second', fontsize=15)      
    plt.errorbar(time[1:],H2_per_site_per_time,xerr=time_uncertainty[1:],yerr=H2_per_site_per_time_uncertainty,fmt='-',label=f"{temperature} K")
#plt.legend(title="Temperature:",fontsize=15,title_fontsize=15)
plt.axvline(x=0,color='black',linestyle='--')
plt.axhline(y=0,color='black',linestyle='--')
plt.xlim(-0.5,10.5)
#plt.ylim(-0.01,0.04)
#plt.show()
  
plt.subplot(1,2,2)
for index in range(len(temperature_list[:5])):
    temperature=temperature_list[index]
    lattice_sites=lattice_size_list[index]**2
    time, time_uncertainty, CO2_amount, CO2_amount_uncertainty = np.loadtxt(f"CO2 Amounts (T={temperature}K,P={np.round(starting_pressure,2)}bar,x(CO)={starting_CO_fraction}).txt",unpack=True)
    CO2_per_site_per_time=np.array(CO2_amount[1:])/(lattice_sites*np.array(time[1:]))
    CO2_per_site_per_time_uncertainty=(1/lattice_sites) * ((np.array(CO2_amount_uncertainty[1:]) / np.array(time[1:])) + ( (np.array(CO2_amount_uncertainty[1:])*np.array(time_uncertainty[1:])) / (np.array(time[1:]))**2 ))
    plt.tick_params(axis='both', right=True, top=True, labelsize=15)
    plt.xlabel('Time Elapsed', fontsize=15)
    plt.ylabel('CO₂ Gas Molecules Per Site Per Second', fontsize=15)      
    plt.errorbar(time[1:],CO2_per_site_per_time,xerr=time_uncertainty[1:],yerr=CO2_per_site_per_time_uncertainty,fmt='-',label=f"{temperature} K")
plt.legend(title="Temperature:",loc="right",bbox_to_anchor =(1.45, 0.5), ncol = 1, fontsize=15,title_fontsize=15)
plt.axvline(x=0,color='black',linestyle='--')
plt.axhline(y=0,color='black',linestyle='--')
plt.xlim(-0.5,10.5)
plt.show()

########################################################################

########################################################################
print("\nARRHENIUS PLOTS")
plt.figure(figsize=(14,5))
plt.subplots_adjust(wspace =0.15)

H2_gradients=[]
H2_gradients_uncertainty=[]
for index in range(len(temperature_list)):
    temperature=temperature_list[index]
    lattice_sites=lattice_size_list[index]**2
    time, time_uncertainty, H2_amount, H2_amount_uncertainty = np.loadtxt(f"H2 Amounts (T={temperature}K,P={np.round(starting_pressure,2)}bar,x(CO)={starting_CO_fraction}).txt",unpack=True)
    H2_per_site=np.array(H2_amount)/lattice_sites
    
    fit_par,cov = np.polyfit(time[6:], H2_per_site[6:], 1, cov=True) #, w=fit_weighting
    grad=fit_par[0]
    H2_gradients.append(grad)
    fit_SE = np.sqrt(np.diag(cov))
    grad_SE = fit_SE[0]
    H2_gradients_uncertainty.append(grad_SE)

plt.subplot(1,2,1)    
plt.tick_params(axis='both', right=True, top=True, labelsize=15)
plt.xlabel('Temperature / K', fontsize=15)
plt.ylabel('∂(H₂ per site)/∂t / s$^{-1}$', fontsize=15)     
plt.errorbar(temperature_list, H2_gradients, yerr=H2_gradients_uncertainty, fmt='kx')
plt.plot(temperature_list, H2_gradients, 'k-')
#plt.show()


ln_H2_gradients=np.log(np.array(H2_gradients))
ln_H2_gradients_uncertainty=np.array(H2_gradients_uncertainty)/np.array(H2_gradients)
inverse_temperature = 1/np.array(temperature_list)


plt.subplot(1,2,2)
plt.errorbar(inverse_temperature*1000, -ln_H2_gradients, yerr=ln_H2_gradients_uncertainty, fmt='kx')

inverse_temperature = np.delete(inverse_temperature, [4,0])
ln_H2_gradients = np.delete(ln_H2_gradients, [4,0])

fit_par,cov = np.polyfit(inverse_temperature, ln_H2_gradients, 1, cov=True)
grad=fit_par[0]
icp=fit_par[1]
fitting_inv_temperature=np.linspace(1.4e-3,2.1e-3,20)
ln_H2_gradient_fitted_line=grad*fitting_inv_temperature+icp
fit_SE = np.sqrt(np.diag(cov))
grad_SE = fit_SE[0]
icp_SE = fit_SE[1]

plt.tick_params(axis='both', right=True, top=True, labelsize=15)
plt.xlabel('1/Temperature / kK$^{-1}$', fontsize=15)
plt.ylabel('-ln(∂(H₂ per site)/∂t)', fontsize=15)
plt.plot(fitting_inv_temperature*1000,-ln_H2_gradient_fitted_line,'k-')
plt.show()

print("This linear plot omits the 500 K and 600 K terms in its fitting, as these involve lower H2 levels than would be expected from the otherwise good linear fitting, likely because they are done on 50x50 rather than 30x30 lattices.")
print("\nThe linear fitting shows that the overall reaction CO + H2O -> H2 + CO2 follows Arrhenius behaviour as each of the individual steps do, with the gradient giving an overall activation energy and the intercept the fabled pre-exp factor (although this includes other factors such as amounts of CO and H2O and likely some constant parameters so is not as interesting).")
print(f"\nThe gradient of this plot is {-grad} K with uncertainty {grad_SE} K. This gives the overall activation energy of the WGSR on Pd as {-grad*kb/e} eV with uncertainty {grad_SE*kb/e} eV.")
print(f"\nFor what it's worth, the intercept of this plot is {-icp} with uncertainty {icp_SE}. This gives effective pre-exponential factor (which contains dependence on reactant amounts and lattice size by the looks of it) of {np.exp(icp)} s^-1 with uncertainty {np.exp(icp)*icp_SE} s^-1.")


CO2_gradients=[]
CO2_gradients_uncertainty=[]
for index in range(len(temperature_list)):
    temperature=temperature_list[index]
    lattice_sites=lattice_size_list[index]**2
    time, time_uncertainty, CO2_amount, CO2_amount_uncertainty = np.loadtxt(f"CO2 Amounts (T={temperature}K,P={np.round(starting_pressure,2)}bar,x(CO)={starting_CO_fraction}).txt",unpack=True)
    CO2_per_site=np.array(CO2_amount)/lattice_sites
    
    fit_par,cov = np.polyfit(time[6:], CO2_per_site[6:], 1, cov=True) #, w=fit_weighting
    grad=fit_par[0]
    CO2_gradients.append(grad)
    fit_SE = np.sqrt(np.diag(cov))
    grad_SE = fit_SE[0]
    CO2_gradients_uncertainty.append(grad_SE)

ln_CO2_gradients=np.log(np.array(CO2_gradients))
ln_CO2_gradients_uncertainty=np.array(CO2_gradients_uncertainty)/np.array(CO2_gradients)
inverse_temperature = 1/np.array(temperature_list)

inverse_temperature = np.delete(inverse_temperature, [4,0])
ln_CO2_gradients = np.delete(ln_CO2_gradients, [4,0])

fit_par,cov = np.polyfit(inverse_temperature, ln_CO2_gradients, 1, cov=True)
grad=fit_par[0]
icp=fit_par[1]
fitting_inv_temperature=np.linspace(1.4e-3,2.1e-3,20)
ln_CO2_gradient_fitted_line=grad*fitting_inv_temperature+icp
fit_SE = np.sqrt(np.diag(cov))
grad_SE = fit_SE[0]
icp_SE = fit_SE[1]

print(f"\nSimilar analysis for CO2 gives an activation barrier of {-grad*kb/e} eV with uncertainty {grad_SE*kb/e} eV and an effective pre-exp factor of {np.exp(icp)} s^-1 with uncertainty {np.exp(icp)*icp_SE} s^-1.")

########################################################

########################################################
#PLOTS OF GAS WITH VARYING CO MOLE FRACTION

print("\n\nVARYING REACTANT GAS COMPOSITION")
CO_fraction_list=[0.01,0.1,0.5,0.9]
lattice_size_list=[30,50,50,50]
starting_pressure=1.0
temperature=600

print("PLOTS PER SITE")
plt.figure(figsize=(12,5))
plt.subplots_adjust(wspace =0.25)

#H2
plt.subplot(1,2,1)
for index in range(len(CO_fraction_list)):
    starting_CO_fraction=CO_fraction_list[index]
    lattice_sites=lattice_size_list[index]**2
    time, time_uncertainty, H2_amount, H2_amount_uncertainty = np.loadtxt(f"H2 Amounts (T={temperature}K,P={np.round(starting_pressure,2)}bar,x(CO)={starting_CO_fraction}).txt",unpack=True)
    H2_per_site=np.array(H2_amount)/lattice_sites
    H2_per_site_uncertainty=np.array(H2_amount_uncertainty)/lattice_sites
    plt.tick_params(axis='both', right=True, top=True, labelsize=15)
    plt.xlabel('Time Elapsed / s', fontsize=15)
    plt.ylabel('Number of H₂ Gas Molecules Per Site', fontsize=15)      
    plt.errorbar(time,H2_per_site,xerr=time_uncertainty,yerr=H2_per_site_uncertainty,fmt='-',label=f"{starting_CO_fraction}")
    #plt.errorbar(time,H2_per_site,fmt='-',label=f"{starting_CO_fraction}")
plt.axvline(x=0,color='black',linestyle='--')
plt.axhline(y=0,color='black',linestyle='--')
plt.xlim(-0.5,10.5)
plt.ylim(-0.025,1)
#plt.legend(title="CO fraction in reactant gas:",fontsize=15,title_fontsize=15)
#plt.show()

#CO2
plt.subplot(1,2,2)
for index in range(len(CO_fraction_list)):
    starting_CO_fraction=CO_fraction_list[index]
    lattice_sites=lattice_size_list[index]**2
    time, time_uncertainty, CO2_amount, CO2_amount_uncertainty = np.loadtxt(f"CO2 Amounts (T={temperature}K,P={np.round(starting_pressure,2)}bar,x(CO)={starting_CO_fraction}).txt",unpack=True)
    CO2_per_site=np.array(CO2_amount)/lattice_sites
    CO2_per_site_uncertainty=np.array(CO2_amount_uncertainty)/lattice_sites
    plt.tick_params(axis='both', right=True, top=True, labelsize=15)
    plt.xlabel('Time Elapsed / s', fontsize=15)
    plt.ylabel('Number of CO₂ Gas Molecules Per Site', fontsize=15)      
    plt.errorbar(time,CO2_per_site,xerr=time_uncertainty,yerr=CO2_per_site_uncertainty,fmt='-',label=f"{starting_CO_fraction}")
plt.legend(title="CO fraction in reactant gas:",loc="upper left",bbox_to_anchor =(-0.7, 1.3), ncol = 4, fontsize=15,title_fontsize=15)
plt.axvline(x=0,color='black',linestyle='--')
plt.axhline(y=0,color='black',linestyle='--')
plt.xlim(-0.5,10.5)
plt.ylim(-0.025,1)
plt.show()



print("\nPLOTS PER SITE PER SECOND")
plt.figure(figsize=(12,5))
plt.subplots_adjust(wspace =0.23)

#H2
plt.subplot(1,2,1)
for index in range(len(CO_fraction_list)):
    starting_CO_fraction=CO_fraction_list[index]
    lattice_sites=lattice_size_list[index]**2
    time, time_uncertainty, H2_amount, H2_amount_uncertainty = np.loadtxt(f"H2 Amounts (T={temperature}K,P={np.round(starting_pressure,2)}bar,x(CO)={starting_CO_fraction}).txt",unpack=True)
    H2_per_site_per_time=np.array(H2_amount[1:])/(lattice_sites*np.array(time[1:]))
    H2_per_site_per_time_uncertainty=(1/lattice_sites) * ((np.array(H2_amount_uncertainty[1:]) / np.array(time[1:])) + ( (np.array(H2_amount_uncertainty[1:])*np.array(time_uncertainty[1:])) / (np.array(time[1:]))**2 ))
    plt.tick_params(axis='both', right=True, top=True, labelsize=15)
    plt.xlabel('Time Elapsed', fontsize=15)
    plt.ylabel('H₂ Gas Molecules Per Site Per Second', fontsize=15)      
    plt.errorbar(time[1:],H2_per_site_per_time,xerr=time_uncertainty[1:],yerr=H2_per_site_per_time_uncertainty,fmt='-',label=f"{starting_CO_fraction}")
    #plt.errorbar(time[1:],H2_per_site_per_time,fmt='-',label=f"{starting_CO_fraction}")
plt.axvline(x=0,color='black',linestyle='--')
#plt.axhline(y=0,color='black',linestyle='--')
plt.xlim(-0.5,10.5)
#plt.legend(title="CO fraction in reactant gas:",fontsize=15,title_fontsize=15)
#plt.show()

#CO2 
plt.subplot(1,2,2)
for index in range(len(CO_fraction_list)):
    starting_CO_fraction=CO_fraction_list[index]
    lattice_sites=lattice_size_list[index]**2
    time, time_uncertainty, CO2_amount, CO2_amount_uncertainty = np.loadtxt(f"CO2 Amounts (T={temperature}K,P={np.round(starting_pressure,2)}bar,x(CO)={starting_CO_fraction}).txt",unpack=True)
    CO2_per_site_per_time=np.array(CO2_amount[1:])/(lattice_sites*np.array(time[1:]))
    CO2_per_site_per_time_uncertainty=(1/lattice_sites) * ((np.array(CO2_amount_uncertainty[1:]) / np.array(time[1:])) + ( (np.array(CO2_amount_uncertainty[1:])*np.array(time_uncertainty[1:])) / (np.array(time[1:]))**2 ))
    plt.tick_params(axis='both', right=True, top=True, labelsize=15)
    plt.xlabel('Time Elapsed', fontsize=15)
    plt.ylabel('CO₂ Gas Molecules Per Site Per Second', fontsize=15)      
    plt.errorbar(time[1:],CO2_per_site_per_time,xerr=time_uncertainty[1:],yerr=CO2_per_site_per_time_uncertainty,fmt='-',label=f"{starting_CO_fraction}")
    #plt.errorbar(time[1:],CO2_per_site_per_time,fmt='-',label=f"{starting_CO_fraction}")
plt.legend(title="CO fraction in reactant gas:",loc="upper left",bbox_to_anchor =(-0.7, 1.3), ncol = 4, fontsize=15,title_fontsize=15)
plt.axvline(x=0,color='black',linestyle='--')
#plt.axhline(y=0,color='black',linestyle='--')
plt.xlim(-0.5,10.5)
plt.show()
##################################################################

########################################################
#PLOTS OF GAS WITH VARYING TOTAL REACTANT GAS PRESSURES

print("\n\nVARYING REACTANT GAS TOTAL PRESSURE")
starting_pressure_list=[0.1,1.0,10.0,100.0]
lattice_size_list=[50,50,50,50]
starting_CO_fraction=0.5
temperature=600

print("PLOTS PER SITE")
plt.figure(figsize=(12,5))
plt.subplots_adjust(wspace =0.25)

#H2
plt.subplot(1,2,1)
for index in range(len(starting_pressure_list)):
    starting_pressure=starting_pressure_list[index]
    lattice_sites=lattice_size_list[index]**2
    time, time_uncertainty, H2_amount, H2_amount_uncertainty = np.loadtxt(f"H2 Amounts (T={temperature}K,P={np.round(starting_pressure,2)}bar,x(CO)={starting_CO_fraction}).txt",unpack=True)
    H2_per_site=np.array(H2_amount)/lattice_sites
    H2_per_site_uncertainty=np.array(H2_amount_uncertainty)/lattice_sites
    plt.tick_params(axis='both', right=True, top=True, labelsize=15)
    plt.xlabel('Time Elapsed / s', fontsize=15)
    plt.ylabel('Number of H₂ Gas Molecules Per Site', fontsize=15)      
    plt.errorbar(time,H2_per_site,xerr=time_uncertainty,yerr=H2_per_site_uncertainty,fmt='-',label=f"{starting_pressure} bar")
    #plt.errorbar(time,H2_per_site,fmt='-',label=f"{starting_CO_fraction}")
plt.axvline(x=0,color='black',linestyle='--')
plt.axhline(y=0,color='black',linestyle='--')
plt.xlim(-0.5,10.5)
plt.ylim(-0.025,0.35)
#plt.legend(title="CO fraction in reactant gas:",fontsize=15,title_fontsize=15)
#plt.show()

#CO2
plt.subplot(1,2,2)
for index in range(len(starting_pressure_list)):
    starting_pressure=starting_pressure_list[index]
    lattice_sites=lattice_size_list[index]**2
    time, time_uncertainty, CO2_amount, CO2_amount_uncertainty = np.loadtxt(f"CO2 Amounts (T={temperature}K,P={np.round(starting_pressure,2)}bar,x(CO)={starting_CO_fraction}).txt",unpack=True)
    CO2_per_site=np.array(CO2_amount)/lattice_sites
    CO2_per_site_uncertainty=np.array(CO2_amount_uncertainty)/lattice_sites
    plt.tick_params(axis='both', right=True, top=True, labelsize=15)
    plt.xlabel('Time Elapsed / s', fontsize=15)
    plt.ylabel('Number of CO₂ Gas Molecules Per Site', fontsize=15)      
    plt.errorbar(time,CO2_per_site,xerr=time_uncertainty,yerr=CO2_per_site_uncertainty,fmt='-',label=f"{starting_pressure} bar")
plt.legend(title="Total reactant gas pressure:",loc="upper left",bbox_to_anchor =(-1, 1.3), ncol = 4, fontsize=15,title_fontsize=15)
plt.axvline(x=0,color='black',linestyle='--')
plt.axhline(y=0,color='black',linestyle='--')
plt.xlim(-0.5,10.5)
plt.ylim(-0.025,0.35)
plt.show()



print("\nPLOTS PER SITE PER SECOND")
plt.figure(figsize=(12,5))
plt.subplots_adjust(wspace =0.3)

#H2
plt.subplot(1,2,1)
for index in range(len(starting_pressure_list)):
    starting_pressure=starting_pressure_list[index]
    lattice_sites=lattice_size_list[index]**2
    time, time_uncertainty, H2_amount, H2_amount_uncertainty = np.loadtxt(f"H2 Amounts (T={temperature}K,P={np.round(starting_pressure,2)}bar,x(CO)={starting_CO_fraction}).txt",unpack=True)
    H2_per_site_per_time=np.array(H2_amount[1:])/(lattice_sites*np.array(time[1:]))
    H2_per_site_per_time_uncertainty=(1/lattice_sites) * ((np.array(H2_amount_uncertainty[1:]) / np.array(time[1:])) + ( (np.array(H2_amount_uncertainty[1:])*np.array(time_uncertainty[1:])) / (np.array(time[1:]))**2 ))
    plt.tick_params(axis='both', right=True, top=True, labelsize=15)
    plt.xlabel('Time Elapsed', fontsize=15)
    plt.ylabel('H₂ Gas Molecules Per Site Per Second', fontsize=15)      
    plt.errorbar(time[1:],H2_per_site_per_time,xerr=time_uncertainty[1:],yerr=H2_per_site_per_time_uncertainty,fmt='-',label=f"{starting_pressure} bar")
    #plt.errorbar(time[1:],H2_per_site_per_time,fmt='-',label=f"{starting_CO_fraction}")
plt.axvline(x=0,color='black',linestyle='--')
#plt.axhline(y=0,color='black',linestyle='--')
plt.xlim(-0.5,10.5)
#plt.legend(title="CO fraction in reactant gas:",fontsize=15,title_fontsize=15)
#plt.show()

#CO2 
plt.subplot(1,2,2)
for index in range(len(starting_pressure_list)):
    starting_pressure=starting_pressure_list[index]
    lattice_sites=lattice_size_list[index]**2
    time, time_uncertainty, CO2_amount, CO2_amount_uncertainty = np.loadtxt(f"CO2 Amounts (T={temperature}K,P={np.round(starting_pressure,2)}bar,x(CO)={starting_CO_fraction}).txt",unpack=True)
    CO2_per_site_per_time=np.array(CO2_amount[1:])/(lattice_sites*np.array(time[1:]))
    CO2_per_site_per_time_uncertainty=(1/lattice_sites) * ((np.array(CO2_amount_uncertainty[1:]) / np.array(time[1:])) + ( (np.array(CO2_amount_uncertainty[1:])*np.array(time_uncertainty[1:])) / (np.array(time[1:]))**2 ))
    plt.tick_params(axis='both', right=True, top=True, labelsize=15)
    plt.xlabel('Time Elapsed', fontsize=15)
    plt.ylabel('CO₂ Gas Molecules Per Site Per Second', fontsize=15)      
    plt.errorbar(time[1:],CO2_per_site_per_time,xerr=time_uncertainty[1:],yerr=CO2_per_site_per_time_uncertainty,fmt='-',label=f"{starting_pressure} bar")
    #plt.errorbar(time[1:],CO2_per_site_per_time,fmt='-',label=f"{starting_CO_fraction}")
plt.legend(title="Total reactant gas pressure:",loc="upper left",bbox_to_anchor =(-1, 1.3), ncol = 4, fontsize=15,title_fontsize=15)
plt.axvline(x=0,color='black',linestyle='--')
#plt.axhline(y=0,color='black',linestyle='--')
plt.xlim(-0.5,10.5)
plt.show()
##################################################################