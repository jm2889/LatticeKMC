from matplotlib import pyplot as plt
from matplotlib import rcParams
import numpy as np

rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
rcParams['errorbar.capsize'] = 3 
rcParams['lines.markersize'] = 10
rcParams['figure.figsize'] = 8, 5

#data and data_uncertainties should be a list of lists of the relevant values
def plots_over_time(time, data, time_uncertainties, data_uncertainties, xtitle, ytitle, data_labels, colour_list, legend_toggle):
    plt.axhline(y=0,color='black',linestyle='--')
    plt.axvline(x=0,color='black',linestyle='--')
    plt.tick_params(axis='both', right=True, top=True, labelsize=15)
    plt.xlabel(xtitle, fontsize=15)
    plt.ylabel(ytitle, fontsize=15)
    for data_list_index in range(len(data)):
        plt.errorbar(time,data[data_list_index],xerr=time_uncertainties,yerr=data_uncertainties[data_list_index],fmt=f'{colour_list[data_list_index]}-',label=data_labels[data_list_index])
        #plt.plot(time,data[data_list_index],f'{colour_list[data_list_index]}-')
    if legend_toggle==True:
        plt.legend(loc="upper center",bbox_to_anchor =(0.5, 1.15), ncol = min(len(data),3), fontsize=15)
    #plt.show()


#User inputs are broken and crash Spyder for some reason on the version Campus uses
temperature=int(input("What temperature is this data taken at (i.e. what temperature is on the file names)?"))
starting_pressure=float(input("What is the total starting pressure of reactant gas?"))
starting_CO_fraction=float(input("What is the starting mole fraction of CO in the reactant gas mixture?"))


#temperature=600
#starting_pressure=1.0
#starting_CO_fraction=0.5
repeat_run=1

######################################################
#INDEXING, AVERAGING AND PLOTTING OF DATA THAT DEPENDS ON TIME

#Values of time (after periodic averaging)
time=np.loadtxt(f"T={temperature}K,P={np.round(starting_pressure,2)}bar,x(CO)={starting_CO_fraction},run={repeat_run} time_values.txt",unpack=True)

index=[np.ceil(time[data_index_num]/0.5) for data_index_num in range(len(time))]

index_set=sorted(list(set(index)))

#Coverage data of each surface species against time
coverage_CO,coverage_H2O,coverage_OH,coverage_H,coverage_O,coverage_COOH,coverage_CO2=np.loadtxt(f"T={temperature}K,P={np.round(starting_pressure,2)}bar,x(CO)={starting_CO_fraction},run={repeat_run} coverage.txt",unpack=True) 

#Amount of product gas present in the simulation over time.
CO_amount,H2O_amount,H2_produced,CO2_produced=np.loadtxt(f"T={temperature}K,P={np.round(starting_pressure,2)}bar,x(CO)={starting_CO_fraction},run={repeat_run} gases.txt",unpack=True)
#Note: could later see about simulating preventing re-reaction of these gases 

mean_time=[None for item in range(len(index_set))]
stdev_time=[None for item in range(len(index_set))]

mean_coverage_CO=[None for item in range(len(index_set))]
stdev_coverage_CO=[None for item in range(len(index_set))]

mean_coverage_H=[None for item in range(len(index_set))]
stdev_coverage_H=[None for item in range(len(index_set))]

mean_coverage_OH=[None for item in range(len(index_set))]
stdev_coverage_OH=[None for item in range(len(index_set))]

mean_H2_produced=[None for item in range(len(index_set))]
stdev_H2_produced=[None for item in range(len(index_set))]

mean_CO2_produced=[None for item in range(len(index_set))]
stdev_CO2_produced=[None for item in range(len(index_set))]

for item in range(len(index_set)):
    current_time=[]
    current_coverage_CO=[]
    current_coverage_H=[]
    current_coverage_OH=[]
    current_H2_produced=[]
    current_CO2_produced=[]
    for meta_index in range(len(index)):
        if index[meta_index]==index_set[item]:
            current_time.append(time[meta_index])
            current_coverage_CO.append(coverage_CO[meta_index])
            current_coverage_H.append(coverage_H[meta_index])
            current_coverage_OH.append(coverage_OH[meta_index])
            current_H2_produced.append(H2_produced[meta_index])
            current_CO2_produced.append(CO2_produced[meta_index])
            
    mean_time[item]=np.mean(current_time)
    stdev_time[item]=np.std(current_time)            
            
    mean_coverage_CO[item]=np.mean(current_coverage_CO)
    stdev_coverage_CO[item]=np.std(current_coverage_CO)

    mean_coverage_H[item]=np.mean(current_coverage_H)
    stdev_coverage_H[item]=np.std(current_coverage_H)
    
    mean_coverage_OH[item]=np.mean(current_coverage_OH)
    stdev_coverage_OH[item]=np.std(current_coverage_OH) 
    
    mean_H2_produced[item]=np.mean(current_H2_produced)
    stdev_H2_produced[item]=np.std(current_H2_produced)    
    
    mean_CO2_produced[item]=np.mean(current_CO2_produced)
    stdev_CO2_produced[item]=np.std(current_CO2_produced)       


plt.figure(figsize=(8,10))
plt.subplots_adjust(hspace =0.2)

plt.subplot(2,1,1)
plt.plot([],[],'r-',label='CO')
plt.plot([],[],'c-',label='H')
plt.plot([],[],'g-',label='OH')
plt.legend(loc="upper center",bbox_to_anchor =(0.5, 1.15), ncol = 3, fontsize=15)
coverage_data=[mean_coverage_CO]
coverage_uncertainties=[stdev_coverage_CO]
colour_list=['r']
surface_species_labels=["CO"]
plots_over_time(mean_time, coverage_data, stdev_time, coverage_uncertainties, "Time Elapsed During Simulation / s", "Coverage Fraction of Species", surface_species_labels, colour_list,False)

plt.subplot(2,1,2)
coverage_data=[mean_coverage_OH,mean_coverage_H]
coverage_uncertainties=[stdev_coverage_OH,stdev_coverage_H]
colour_list=['g','c']
surface_species_labels=["OH", "H"]
plots_over_time(mean_time, coverage_data, stdev_time, coverage_uncertainties, "Time Elapsed During Simulation / s", "Coverage Fraction of Species", surface_species_labels, colour_list,False)

plt.savefig(f"Coverage Plot(T={temperature}K,P={np.round(starting_pressure,2)}bar,x(CO)={starting_CO_fraction}).png",bbox_inches='tight')
plt.show()


plt.figure(figsize=(12,5))
plt.subplots_adjust(wspace =0.3)

plt.subplot(1,2,1)
plt.plot([],[],'b-',label="H₂")
plt.plot([],[],'r-',label="CO₂")
plt.legend(loc="upper center",bbox_to_anchor =(1.1, 1.15), ncol = 2, fontsize=15)
colour_list=['b']
product_gas_data=[mean_H2_produced]
product_gas_uncertainties=[stdev_H2_produced]
product_gases_labels=["H₂"]
plots_over_time(mean_time, product_gas_data, stdev_time, product_gas_uncertainties, "Time Elapsed During Simulation / s", "Total Number of Gas Molecules", product_gases_labels,colour_list,False)

plt.subplot(1,2,2)
colour_list=['r']
product_gas_data=[mean_CO2_produced]
product_gas_uncertainties=[stdev_CO2_produced]
product_gases_labels=["CO₂"]
plots_over_time(mean_time, product_gas_data, stdev_time, product_gas_uncertainties, "Time Elapsed During Simulation / s", "Total Number of Gas Molecules", product_gases_labels,colour_list,False)

plt.savefig(f"Product Gas Plot(T={temperature}K,P={np.round(starting_pressure,2)}bar,x(CO)={starting_CO_fraction}).png",bbox_inches='tight')
plt.show()
############################################################

############################################################
#PROCESSING AND PLOTTING OF REACTION FREQUENCY DATA (NOT DEPENDENT ON TIME)

#Number of times each reaction/process occurred during the simulation (should be a relatively small file with one number for each reaction and its reverse)
reaction_frequencies=np.loadtxt(f"T={temperature}K,P={np.round(starting_pressure,2)}bar,x(CO)={starting_CO_fraction},run={repeat_run} processes.txt",unpack=True)
#Assumes that it is saved in the order: reaction 1, -1, 2 , -2 etc. These are also in the same order as in the Chutia paper

mean_reaction_frequencies=[None for reaction_index in range(len(reaction_frequencies))]
stdev_reaction_frequencies=[None for reaction_index in range(len(reaction_frequencies))]
    
for reaction_index in range(len(reaction_frequencies)):
    mean_reaction_frequencies[reaction_index]=np.mean(reaction_frequencies[reaction_index])
    stdev_reaction_frequencies[reaction_index]=np.std(reaction_frequencies[reaction_index])

#Convoluted way of splitting the reaction frequencies into forwards and reverse reactions
forward_index_list=[0,2,4,6,8,10,12,14,16,18,19,20,21,22] #List of indices of forward reactions (inc. diffusion)
reverse_index_list=[1,3,5,7,9,11,13,15,17,None,None,None,None,23]
forward_frequencies=[]
reverse_frequencies=[]
forward_uncertainties=[]
reverse_uncertainties=[]
for forward_index in forward_index_list:
    forward_frequencies.append(mean_reaction_frequencies[forward_index])
    forward_uncertainties.append(stdev_reaction_frequencies[forward_index])
for reverse_index in reverse_index_list:
    if reverse_index==None:
        reverse_frequencies.append(0) #Diffusion steps have no reverse, so just represent these reverses as happening 0 times
        reverse_uncertainties.append(0)
    else:
        reverse_frequencies.append(mean_reaction_frequencies[reverse_index])
        reverse_uncertainties.append(stdev_reaction_frequencies[reverse_index])
    
forward_index=np.arange(0,len(forward_frequencies))
reverse_index=np.arange(0,len(reverse_frequencies))

process_name=['H₂O ⇌ OH + H','OH + OH ⇌\nH₂O + O','OH ⇌ O + H','H + H ⇌ H₂','CO + O ⇌ CO₂','CO + OH ⇌\nCOOH','COOH ⇌\nCO₂ + H','CO Adsorption','H₂O Adsorption','O Diffusion','H Diffusion','OH Diffusion','CO Diffusion','CO₂ Desorption']
rcParams['figure.figsize'] = 12, 5
plt.bar(forward_index, forward_frequencies, 0.2, yerr=forward_uncertainties, color='cyan', edgecolor='blue', label='Forward reactions')
plt.bar(reverse_index+0.2, reverse_frequencies, 0.2, yerr=reverse_uncertainties, color='orange', edgecolor='red', label='Reverse reactions')
if max(mean_reaction_frequencies) > 1000:
    plt.yscale('log')
plt.xticks(forward_index+0.1, process_name, fontsize=15, rotation=70)
plt.ylabel('Number of Occurrences of Each Process', fontsize=15)
plt.yticks(fontsize=15)
plt.legend(loc="upper center",bbox_to_anchor =(0.5, 1.15), ncol = 3, fontsize=15)
plt.savefig(f"All Process Plot(T={temperature}K,P={np.round(starting_pressure,2)}bar,x(CO)={starting_CO_fraction}).png",bbox_inches='tight')
plt.show()

process_name=['H₂O ⇌ OH + H','OH + OH ⇌\nH₂O + O','OH ⇌ O + H','H + H ⇌ H₂','CO + O ⇌ CO₂','CO + OH ⇌\nCOOH','COOH ⇌\nCO₂ + H']
rcParams['figure.figsize'] = 8, 5
plt.bar(forward_index[:7], forward_frequencies[:7], 0.2, yerr=forward_uncertainties[:7], color='cyan', edgecolor='blue', label='Forward reactions')
plt.bar(reverse_index[:7]+0.2, reverse_frequencies[:7], 0.2, yerr=reverse_uncertainties[:7], color='orange', edgecolor='red', label='Reverse reactions')
if max(mean_reaction_frequencies[:13]) > 1000:
    plt.yscale('log')
plt.xticks(forward_index[:7]+0.1, process_name, fontsize=15, rotation=55)
plt.ylabel('Number of Occurrences of Each Process', fontsize=15)
plt.yticks(fontsize=15)
plt.legend(loc="upper center",bbox_to_anchor =(0.5, 1.15), ncol = 3, fontsize=15)
plt.savefig(f"Reaction Process Plot(T={temperature}K,P={np.round(starting_pressure,2)}bar,x(CO)={starting_CO_fraction}).png",bbox_inches='tight')
plt.show()
###############################################################

#Save a file containing the H2 amounts over time at this temperature for use in more data analysis code that will come later (that will look at H2 production over all temperatures)
f=open(f"H2 Amounts (T={temperature}K,P={np.round(starting_pressure,2)}bar,x(CO)={starting_CO_fraction}).txt",'w')
for time_index in range(len(mean_time)):
    f.writelines(str(mean_time[time_index]) + " " + str(stdev_time[time_index])+ " " + str(mean_H2_produced[time_index])+ " " + str(stdev_H2_produced[time_index]))
    f.writelines("\n")
f.close()

#Save a file containing the H2 amounts over time at this temperature for use in more data analysis code that will come later (that will look at CO2 production over all temperatures)
f=open(f"CO2 Amounts (T={temperature}K,P={np.round(starting_pressure,2)}bar,x(CO)={starting_CO_fraction}).txt",'w')
for time_index in range(len(mean_time)):
    f.writelines(str(mean_time[time_index]) + " " + str(stdev_time[time_index])+ " " + str(mean_CO2_produced[time_index])+ " " + str(stdev_CO2_produced[time_index]))
    f.writelines("\n")
f.close()
