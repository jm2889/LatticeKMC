# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 10:54:55 2023

@author: aljom
"""

from matplotlib import pyplot as plt
from matplotlib import rcParams
import numpy as np

rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
rcParams['errorbar.capsize'] = 3 
rcParams['lines.markersize'] = 10
rcParams['figure.figsize'] = 8, 5


temperature = 600
starting_pressure = 1.0
starting_CO_fraction = 0.5
repeat_run = 1


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
process_name=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14']
plt.bar(forward_index, forward_frequencies, 0.2, yerr=forward_uncertainties, color='cyan', edgecolor='blue', label='Forward reactions')
plt.bar(reverse_index+0.2, reverse_frequencies, 0.2, yerr=reverse_uncertainties, color='orange', edgecolor='red', label='Reverse reactions')
if max(mean_reaction_frequencies) > 1000:
    plt.yscale('log')
plt.xticks(forward_index+0.1, process_name, fontsize=15)# rotation=70)
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