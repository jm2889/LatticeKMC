#O (0 in all_lattice_populations)
#2 = H2O
#3 = OH
#4 = H
#5 = O
#6 = H2 (not used-gas)
#7 = CO2 (not used-gas)
#8 = COOH (5 in all_lattice_populations)
#9 = CO2 srf (6 in all_lattice_populations)

import random
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rcParams

########################################
#THESE ARE THE PARAMETERS TO ADJUST AND TOGGLE

##Simulation parameters
d = 50 #Length of the lattice (in both directions) in terms of number of sites
max_simulation_time=10 #The length of simulated time for which it will go on for (in s)
number_of_runs = 3 # how many runs at a given set of parameters
number_of_temperatures = 3 # how many different temperatures you want to run in one go
temperature_increment = 25 


##Condition parameters
temperature = 600  #In K
starting_CO_fraction = 0.5 #What fraction of the starting reactant gas mixture (of CO and H2O) is CO?
starting_pressure = 1e5  #Total pressure of both gases at time, t=0 (1e5 Pa = 1 bar)    

##Qualitative toggles
replenish_reactant_gas=True #if True, keep the number of reactant gas molecules the same regardless of if they react at all. This also means not increasing the number id they desorb before reacting.
supress_reverse_product_reaction=True #if True, stop H2 and CO2 from being able to undergo their reverse reactions (always set entry 7 and 9 in rate_constant_matrix to 0). Also stop CO2 readsorption.

##Parameters and toggles for periodic averaging and saving of data
repeat_run=1 #This is just for file names so we can store runs of data taken at the same conditions
save_image=True #Will images be saved at the times given in the list below?
image_saving_times=[0.1, 1, 10] #Times (in s) at which a snapshot of the lattice will be saved
period=0.5 #Determines the lengths of time over which the data is averaged (larger period means fewer data points generated from averaging more raw data points over longer time ranges)

##Stiffness scaling parameters
scaling_factor = 10
delta = 0.5

########################################

number_of_sites=d**2

rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
rcParams['errorbar.capsize'] = 3 
rcParams['lines.markersize'] = 10
rcParams['figure.figsize'] = 10, 5

lattice = [[0 for i in range(d)] for j in range(d)] #Starts off completely empty

#Put CO at one corner and H2O at the furthest point away to get things going but get as close to starting empty as possible (without actually doing so, as this would mean nothing could ever happen)
lattice[0][0]=1
lattice[round(d/2)][round(d/2)]=2

#some constants, or at least terms that should probably not be changed
kb = 1.381e-23  #Boltzmann constant (for convenience)
e = 1.602e-19  #Conversion between eV and J (for convenience)
pre_exp_factor = 5e12  #Pre-exponential factor needed for Arrhenius equation for rate constants of desorption steps (adsorption uses more complicated, gas-dependent equation)
Pd_atom_spacing = 2.751e-10  #Distance between adjacent Pd atoms (so adjacent sites)
CO_mass = 28*1.661e-27  #Mass of CO molecule (in kg)
H2O_mass = 18*1.661e-27
CO2_mass = 44*1.661e-27
H2_mass = 2*1.661e-27
height_of_gas = 5e-6  #Height above the Pd surface which contains the gas that can adsorb (in m) (somewhat arbitrary)    

starting_H2O_fraction = 1-starting_CO_fraction  #Should not be changed


number_of_CO = round((starting_CO_fraction*starting_pressure*(d**2)*(Pd_atom_spacing**2)*height_of_gas)/(kb*temperature))
number_of_H2O = round((starting_H2O_fraction*starting_pressure*(d**2)*(Pd_atom_spacing**2)*height_of_gas)/(kb*temperature))

images_saved=0 #A counter to cycle through the times in the list above

# keep track of the amount of gas released
number_of_CO2 = 0
number_of_H2 = 0

#List of activation barriers for all reaction and diffusion steps (in the same order as before)
activation_barrier_matrix=[1.320,0.240,0.880,1.410,0.320,0.350,0.810,0.018,0.690,0.439,0.720,1.471,0.540,0.580,0.09,0.4,0.017,0]
#For entry 7 (reverse reaction of H2) and 9 (reverse of CO2) which require these species straight from the gas phase, the rate constants should be treated like those for adsorption in terms of dependence on number of gas molecules

#Rate constant list now varies with temperature
rate_constant_matrix = [pre_exp_factor * np.exp(-(activation_barrier_matrix[index]*e)/(kb*temperature)) for index in range(len(activation_barrier_matrix))]
if supress_reverse_product_reaction==True:
    rate_constant_matrix[7]=0
    rate_constant_matrix[9]=0
else:
    #Reverse H2 reaction
    rate_constant_matrix[7] = (number_of_H2/(number_of_sites*height_of_gas)) * np.sqrt((kb*temperature)/(2*np.pi*H2_mass)) * np.exp(-activation_barrier_matrix[7]*e/(kb*temperature))
    #Reverse CO2 reaction
    rate_constant_matrix[9] = (number_of_CO2/(number_of_sites*height_of_gas)) * np.sqrt((kb*temperature)/(2*np.pi*CO2_mass)) * np.exp(-activation_barrier_matrix[9]*e/(kb*temperature))

#Counts the number of times each process occurs
#Key: regular reactions, CO and H2O ads/des, diffusion, CO2 des/ads (same order as in Chutia paper)
process_count=[0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,  0,0, 0,0,  0,0,0,0,  0,0] 


f = open(f"T={temperature}K,P={np.round(starting_pressure/1e5,2)}bar,x(CO)={starting_CO_fraction},run={repeat_run} time_values.txt", "w")
f.close()


f = open(f"T={temperature}K,P={np.round(starting_pressure/1e5,2)}bar,x(CO)={starting_CO_fraction},run={repeat_run} coverage.txt", "w")
f.close()


f = open(f"T={temperature}K,P={np.round(starting_pressure/1e5,2)}bar,x(CO)={starting_CO_fraction},run={repeat_run} gases.txt", "w")
f.close()

print("")
f = open(f"T={temperature}K,P={np.round(starting_pressure/1e5,2)}bar,x(CO)={starting_CO_fraction},run={repeat_run} processes.txt", "w")
f.close()


#####################################################
#Simple visualisation of the lattice at that current point in time
def latticePlot(lattice,all_changed_i,all_changed_j,show_changed,save_image):
    #Lists of x(j) and y(i) positions of each species
    all_empties_x=[]
    all_empties_y=[]
    all_CO_x=[]
    all_CO_y=[]
    all_H2O_x=[]
    all_H2O_y=[]
    all_OH_x = []
    all_OH_y = []
    all_H_x, all_H_y = [], []
    all_O_x, all_O_y = [], []
    all_COOH_x, all_COOH_y = [],[]
    all_CO2_x, all_CO2_y = [],[]

    #Fills the above lists based on the value at the position j,i
    for i in range(d):
        for j in range(d):
            if lattice[i][j]==0:
                #Fill empty site coordinate lists
                all_empties_x.append(j)
                all_empties_y.append(i)
            elif lattice[i][j]==1:
                #Fill CO coordinate lists
                all_CO_x.append(j)
                all_CO_y.append(i)
            elif lattice[i][j]==2:
                #Fill H2O coordinate lists
                all_H2O_x.append(j)
                all_H2O_y.append(i)
            elif lattice[i][j] == 3:
                all_OH_x.append(j)
                all_OH_y.append(i)
            elif lattice[i][j] == 4:
                all_H_x.append(j)
                all_H_y.append(i)
            elif lattice[i][j] == 5:
                #print("plot O")
                all_O_x.append(j)
                all_O_y.append(i)
            elif lattice[i][j] == 8:
                #print("plot COOH")
                all_COOH_x.append(j)
                all_COOH_y.append(i)
            elif lattice[i][j] == 9:
                #print("plot COOH")
                all_CO2_x.append(j)
                all_CO2_y.append(i)

    #Plots the final lattice using different colours for each species
    rcParams['figure.figsize'] = 5, 5
    plt.tick_params(axis='both', right=True, top=True, labelsize=11)
    plt.axis([-0.5,d-0.5,d-0.5,-0.5])
    plt.xticks([])
    plt.yticks([])
    for site_axes in range(d): #Gives the plots the 'grid' effect
        plt.axhline(y=site_axes+0.5, color='black', linewidth=10/d)
        plt.axvline(x=site_axes+0.5, color='black', linewidth=10/d)
    if show_changed==True:
        plt.plot(all_changed_j, all_changed_i, 'kx', markersize=200/d, label='Changed')
    plt.plot(all_CO_x, all_CO_y, 'r.', markersize=200/d, label='Adsorbed CO')
    plt.plot(all_H2O_x, all_H2O_y, 'b.', markersize=200/d, label='Adsorbed H2O')
    plt.plot(all_OH_x, all_OH_y, 'g.', markersize=200/d, label='Adsorbed OH')
    plt.plot(all_H_x, all_H_y, 'c.', markersize=200/d, label='Adsorbed H')
    plt.plot(all_O_x, all_O_y, 'm.', markersize=200/d, label='Adsorbed O')
    plt.plot(all_COOH_x, all_COOH_y, 'k.', markersize=200/d, label='Adsorbed COOH')
    plt.plot(all_CO2_x, all_CO2_y, 'y.', markersize=200/d, label='Adsorbed CO2')
    plt.legend(fontsize=11,loc="upper center",bbox_to_anchor =(0.5, 1.21), ncol = 3, markerscale=d/20)
    
    #if save_image==True:
     #   plt.savefig(f"T={temperature}K,P={np.round(starting_pressure/1e5,2)}bar,x(CO)={starting_CO_fraction},run={repeat_run} lattice_at_{image_saving_times[images_saved]}s.png",bbox_inches='tight')
    
    plt.show()
#####################################################

#####################################################
#Determines which of the four processes can occur based on the identity of the species (or lack thereof) at that lattice site
#These are the so called selection rules that would need to be changed if simulating a different reaction
def whichProcess(species, particles_around, number_of_CO2, number_of_H2, i, j, position_of_processes_that_can_happen):
    reaction_possible = 0
    if ((species == 2) and (0 in particles_around)) or ((species == 0) and (2 in particles_around)): # if there is a water and empty site
        processes_that_can_happen.append(0)
        reaction_possible = 1
    if ((species == 3) and (4 in particles_around)) or ((species == 4) and (3 in particles_around)): # if there is an OH and H next to each each other
        processes_that_can_happen.append(1)
        reaction_possible = 1
    if ((species == 3) and (3 in particles_around)): # if two OH next to one another
        processes_that_can_happen.append(2)
        reaction_possible = 1
    if ((species == 2) and (5 in particles_around)) or ((species == 5) and (2 in particles_around)): #H20 and O
        processes_that_can_happen.append(3)
        reaction_possible = 1
    if ((species == 3) and (0 in particles_around)) or ((species == 0) and (3 in particles_around)): #Oh and empty
        processes_that_can_happen.append(4)
        reaction_possible = 1
    if ((species == 5) and (4 in particles_around)) or ((species == 4) and (5 in particles_around)): # O and H
        processes_that_can_happen.append(5)
        reaction_possible = 1
    if ((species == 4) and (4 in particles_around)): # two H next to each other
        processes_that_can_happen.append(6)
        reaction_possible = 1
# this next step deals with H2 reverse reaction, and checks if any H2 has been produced AND there is an adjacent empty site
    if supress_reverse_product_reaction==False:
        if species == 0 and number_of_H2 != 0 and (0 in surround(lattice, i, j, d)):
            processes_that_can_happen.append(7) 
            reaction_possible = 1
    if ((species == 1) and (5 in particles_around)) or ((species == 5) and (1 in particles_around)): # CO and O
        processes_that_can_happen.append(8)
        reaction_possible = 1
    # CO2 reverse reaction - checks if Co2 is available and an extra empty site
    if supress_reverse_product_reaction==False:
        if ((species == 0) and number_of_CO2 != 0 and (0 in particles_around)):
            processes_that_can_happen.append(9)
            reaction_possible = 1
    if ((species == 1) and (3 in particles_around)) or ((species == 3) and (1 in particles_around)): #CO and OH
        processes_that_can_happen.append(10)
        
        reaction_possible = 1
    if ((species == 8) and (0 in particles_around)) or ((species == 0) and (8 in particles_around)): # COOH and empty
        processes_that_can_happen.append(11)
        processes_that_can_happen.append(12)
        reaction_possible = 1
    if ((species == 9) and (4 in particles_around)) or ((species == 4) and (9 in particles_around)): #surface CO2 and H
        processes_that_can_happen.append(13)
        reaction_possible = 1
    
    if (species==1 and 0 in particles_around):
        no_lat_interactions=[]
        for suitable_site in suitableSurroundings(lattice,i,j,0):
            if superSurround(lattice,suitable_site[0],suitable_site[1],d).count(1) < 2:
                no_lat_interactions.append(suitable_site)
        if len(no_lat_interactions) != 0:                        
            processes_that_can_happen.append(14)
            reaction_possible = 1
    if (species==0 and 1 in particles_around):
        if superSurround(lattice,i,j,d).count(1) < 2:
            processes_that_can_happen.append(14)
            reaction_possible = 1            
      #CO diff
        
    if (species==3 and 0 in particles_around) or (species==0 and 3 in particles_around):
       processes_that_can_happen.append(15)
       reaction_possible = 1
      #OH diff 
       
    if (species==4 and 0 in particles_around) or (species==0 and 4 in particles_around):
        processes_that_can_happen.append(16)
        reaction_possible = 1
      #H diff        
        
    if (species==5 and 0 in particles_around) or (species==0 and 5 in particles_around):
        processes_that_can_happen.append(17)
        reaction_possible = 1
      #O diff
    
    if reaction_possible == 1:
        position_of_processes_that_can_happen.append([i, j])
     
    return processes_that_can_happen, position_of_processes_that_can_happen #This is a list of indices for some of the four reactions (rate constants), not of constants themselves
#################################################################  

#####################################################
#Standard KMC code. The heart and soul of the programme if you want to be poetic about it
def KMC_Sweep(rate_constants, position_of_selected_process):
    partial_sum = 0
    selected_process = 0
    partial_sum += rate_constants[selected_process]
    decision_rng = random.uniform(0, 1)
    total_rate_constant = sum(rate_constants)
    rng_ktot = decision_rng*total_rate_constant
    
    while rng_ktot > partial_sum: #Standard KMC procedure. Stops counting through processes when it reaches the one with probability encompassing the random number generated
        selected_process += 1
        if selected_process == len(rate_constants) - 1:
            #selected_process -= 1
            return selected_process
        else:
            partial_sum += rate_constants[selected_process]
    return selected_process # returns which path we go down, corresponding to which rate const was picked
##################################################### 

#####################################################
#look for a given particle
def look_for_particle(surround, a):
    desired_particle_around = []
    if a == surround[0]: # particle looking for is above
        desired_particle_around.append(0)
    if a == surround[1]: # particle is below
        desired_particle_around.append(1)
    if a == surround[2]: # left
        desired_particle_around.append(2)
    if a == surround[3]: # right
        desired_particle_around.append(3)
    elif a not in surround:
        return 0
    random_particle_around = random.randint(0, len(desired_particle_around))
    if random_particle_around == 0:      
        return desired_particle_around[0]
    else:
        pos_to_return = desired_particle_around[random_particle_around-1]
        return pos_to_return
#####################################################

#####################################################
#look for the locations of a chosen particle around the current (or selected) site. Currently only used for the diffusion biasing of CO
def suitableSurroundings(lattice,i,j,chosen_species):
    locations_of_desired_species=[]
    
    #Above
    if i == 0:
        if lattice[d-1][j]==chosen_species:
            locations_of_desired_species.append([d - 1,j])
    else:
        if lattice[i - 1][j]==chosen_species:
            locations_of_desired_species.append([i - 1,j])
            
    #Below        
    if i == d - 1:
        if lattice[0][j]==chosen_species:
            locations_of_desired_species.append([0,j])        
    else:
        if lattice[i + 1][j]==chosen_species:
            locations_of_desired_species.append([i + 1,j])            

    #Left
    if j == 0:
        if lattice[i][d - 1]==chosen_species:
            locations_of_desired_species.append([i,d - 1])   
    else:
        if lattice[i][j - 1]==chosen_species:
            locations_of_desired_species.append([i,j - 1])           
    
    #Right    
    if j == d - 1:
        if lattice[i][0]==chosen_species:
            locations_of_desired_species.append([i,0])           
    else:
        if lattice[i][j + 1] ==chosen_species:
            locations_of_desired_species.append([i,j + 1] )           

    return locations_of_desired_species
#####################################################

#####################################################
# updates the lattice position for relevant neighbouring sites
def update_above(i, j, particle_to_update_to):
    if i == 0:
        lattice[d-1][j] = particle_to_update_to
        i_changed = d-1
    else:
        lattice[i-1][j] = particle_to_update_to
        i_changed = i-1
    j_changed = (j)
    return i_changed, j_changed


def update_below(i, j, particle_to_update_to):
    if i == d-1:
        lattice[0][j] = particle_to_update_to
        i_changed = (0)
    else:
        lattice[i+1][j] = particle_to_update_to
        i_changed = (i+1)
    j_changed = (j)
    return i_changed, j_changed

def update_left(i, j, particle_to_update_to):
    if j == 0:
        lattice[i][d-1] = particle_to_update_to
        j_changed = (d-1)
    else:
        lattice[i][j-1] = particle_to_update_to
        j_changed = (j-1)
    i_changed = (i)
    return i_changed, j_changed


def update_right(i, j, particle_to_update_to):
    if j == d-1:
        lattice[i][0] = particle_to_update_to
        j_changed = (0)
    else:
        lattice[i][j+1] = particle_to_update_to
        j_changed = (j+1)
    i_changed = (i)
    return i_changed, j_changed


def pos_update(pos, i, j, update_to):
    if pos == 0:
        i_changed, j_changed = update_above(i, j, update_to)
    if pos == 1:
        i_changed, j_changed = update_below(i, j, update_to)
    if pos == 2:
        i_changed, j_changed = update_left(i, j, update_to)
    if pos == 3:
        i_changed, j_changed = update_right(i, j, update_to)
    return i_changed, j_changed
#####################################################

#####################################################
#Performs one of the more than four associated changes depending on the selected process and possible proces(ses) at that site
#Would need to be changed when simulating another reaction
def siteChange(all_lattice_populations, selected_process ,processes_that_can_happen,i,j, species, particles_around, number_of_CO2, number_of_H2, position_of_selected_process):
    changed_i=[]
    changed_j=[]
    
    i = position_of_selected_process[0]
    j = position_of_selected_process[1]
    species = lattice[i][j]
    
    process_number = selected_process
    
    
    if process_number == 0: #Reaction 1   
        process_count[0]+=1
        if species == 2:
            lattice[i][j] = 3 # change the H2O to OH
            position = look_for_particle(surround(lattice, i, j, d), 0)  # look for empty site
            i_changed, j_changed = pos_update(position, i, j, 4) # update empty to H
        elif species == 0:
            lattice[i][j] = 3 # change the empty to OH
            position = look_for_particle(surround(lattice, i, j, d), 2)  # look for H2O
            i_changed, j_changed = pos_update(position, i, j, 4) # update H2O to H
        all_lattice_populations[1]-=1
        all_lattice_populations[2]+=1
        all_lattice_populations[3]+=1
            
    elif process_number == 1: # Reaction -1
        process_count[1]+=1
        if species == 3:
            lattice[i][j] = 2 # change the OH to H2O
            position = look_for_particle(surround(lattice, i, j, d), 4)  # look for H
            i_changed, j_changed = pos_update(position, i, j, 0) # update H to empty
        elif species == 4:
            lattice[i][j] = 2 # change the H to H2O
            position = look_for_particle(surround(lattice, i, j, d), 3)  # look for OH
            i_changed, j_changed = pos_update(position, i, j, 0) # update OH to empty
        all_lattice_populations[1]+=1
        all_lattice_populations[2]-=1
        all_lattice_populations[3]-=1
    
    elif process_number == 2: # Reaction 2
        process_count[2]+=1
        lattice[i][j] = 2 # set one of the OH to H2O
        position = look_for_particle(surround(lattice, i, j, d), 3) # look for the other OH
        i_changed, j_changed = pos_update(position, i, j, 5) # update to O
        all_lattice_populations[2]-=2
        all_lattice_populations[1]+=1
        all_lattice_populations[4]+=1
    
    elif process_number == 3: # Reaction -2
        process_count[3]+=1
        if species == 2:
            lattice[i][j] = 3 # update H2O to OH
            position = look_for_particle(surround(lattice, i, j, d), 5)  # look for O
            i_changed, j_changed = pos_update(position, i, j, 3) # update O to OH
        elif species == 5:
            lattice[i][j] = 3 # update O to OH
            position = look_for_particle(surround(lattice, i, j, d), 2)  # look for H2O
            i_changed, j_changed = pos_update(position, i, j, 3) # update H2O to OH
        all_lattice_populations[2]+=2
        all_lattice_populations[1]-=1
        all_lattice_populations[4]-=1            
        
    elif process_number == 4: # Reaction 3
        process_count[4]+=1
        if species == 3:
            lattice[i][j] = 4 # OH -> H
            position = look_for_particle(surround(lattice, i, j, d), 0) # look for empty
            i_changed, j_changed = pos_update(position, i, j, 5) # empty -> O
        elif species == 0:
            lattice[i][j] = 4 # Empty -> H
            position = look_for_particle(surround(lattice, i, j, d), 3) # look for OH
            i_changed, j_changed = pos_update(position, i, j, 5) # OH -> O
        all_lattice_populations[2]-=1
        all_lattice_populations[4]+=1
        all_lattice_populations[3]+=1            

    elif process_number == 5: # Reaction -3
        process_count[5]+=1
        if species == 4:
            lattice[i][j] = 3 # update H to OH
            position = look_for_particle(surround(lattice, i, j, d), 5) # Look for O
            i_changed, j_changed = pos_update(position, i, j, 0) # update O to empty
        elif species == 5:
            lattice[i][j] = 3 # update O to OH
            position = look_for_particle(surround(lattice, i, j, d), 4) # Look for H
            i_changed, j_changed = pos_update(position, i, j, 0) # update H to empty
        all_lattice_populations[2]+=1
        all_lattice_populations[4]-=1
        all_lattice_populations[3]-=1             
    
    elif process_number == 6: # Reaction 4
        process_count[6]+=1
        lattice[i][j] = 0 # empty the H site
        position = look_for_particle(surround(lattice, i, j, d), 4)
        i_changed, j_changed = pos_update(position, i, j, 0) # empty the other site as well
        number_of_H2 += 1
        all_lattice_populations[3]-=2       
    
    #This is not possible (rate constant set to 0) when supress_reverse_product_reaction is True
    elif process_number == 7: # Reaction -4, H2 re-reaction
        process_count[7]+=1
        lattice[i][j] = 4 # update empty to H
        position = look_for_particle(surround(lattice, i, j, d), 0) # look for the other empty site
        i_changed, j_changed = pos_update(position, i, j, 4) # update the other empty site to H
        number_of_H2 -= 1
        all_lattice_populations[3]+=2

        
    elif process_number == 8: # Reaction 5, 
        process_count[8]+=1
        if species == 1:
            lattice[i][j] = 0 # empty the CO site
            position = look_for_particle(surround(lattice, i, j, d), 5)
            i_changed, j_changed = pos_update(position, i, j, 0) # empty the O site
        elif species == 5:
            lattice[i][j] = 0 # empty the O site
            position = look_for_particle(surround(lattice, i, j, d), 1) # look for the CO site
            i_changed, j_changed = pos_update(position, i, j, 0) # empty the Co site
        number_of_CO2 += 1
        all_lattice_populations[0]-=1
        all_lattice_populations[4]-=1         
    
    #This is not possible (rate constant set to 0) when supress_reverse_product_reaction is True
    elif process_number == 9: # Reaction -5, CO2 re-reaction (not via adsorption)
        process_count[9]+=1
        lattice[i][j] = 1 # empty -> co2
        position = look_for_particle(surround(lattice, i, j, d), 0) # look for the other empty site
        i_changed, j_changed = pos_update(position, i, j, 5)   # fill this empty site with O
        number_of_CO2 -= 1
        all_lattice_populations[0]+=1
        all_lattice_populations[4]+=1 
        
    elif process_number == 10: # Reaction 6
        
        process_count[10]+=1
        if species == 1:
            lattice[i][j] = 8 # update CO to COOH
            position = look_for_particle(surround(lattice, i, j, d), 3) # look for OH
            i_changed, j_changed = pos_update(position, i, j, 0)   # empty this site
        elif species == 3:
            lattice[i][j] = 8 # update OH to COOH
            position = look_for_particle(surround(lattice, i, j, d), 1) # look for CO
            i_changed, j_changed = pos_update(position, i, j, 0)   # empty this site
        all_lattice_populations[0]-=1
        all_lattice_populations[2]-=1
        all_lattice_populations[5]+=1            
    
    elif process_number == 11: # Reaction -6
        process_count[11]+=1
        if species == 8:
            lattice[i][j] = 1 # update COOH to Co
            position = look_for_particle(surround(lattice, i, j, d), 0) # look for empty
            i_changed, j_changed = pos_update(position, i, j, 3)   # update to OH
        elif species == 0:
            lattice[i][j] = 1 # update empty to CO
            position = look_for_particle(surround(lattice, i, j, d), 8) # look for COOH
            i_changed, j_changed = pos_update(position, i, j, 3)   # update to OH
        all_lattice_populations[0]+=1
        all_lattice_populations[2]+=1
        all_lattice_populations[5]-=1             
            
    elif process_number == 12: # Reaction 7
        process_count[12]+=1
        if species == 8:
            lattice[i][j] = 9 # update COOH to CO2(srf)
            position = look_for_particle(surround(lattice, i, j, d), 0) # look for empty
            i_changed, j_changed = pos_update(position, i, j, 4)   # update to H
        elif species == 0:
            lattice[i][j] = 9 # update empty to CO2(srf)
            position = look_for_particle(surround(lattice, i, j, d), 8) # look for COOH
            i_changed, j_changed = pos_update(position, i, j, 4)   # update to H
        all_lattice_populations[5]-=1
        all_lattice_populations[6]+=1
        all_lattice_populations[3]+=1         
    
    #Even though CO2 is undergoing a reverse reaction here, it is stll allowed even when supress_reverse_product_reaction is True, as this CO2 has not been evolved as a gas yet.
    elif process_number == 13: # Reaction -7
        process_count[13]+=1
        if species == 9:
            lattice[i][j] = 8 # update CO2 srf to COOH
            position = look_for_particle(surround(lattice, i, j, d), 4) # look for H
            i_changed, j_changed = pos_update(position, i, j, 0)   # update to empty
        elif species == 4:
            lattice[i][j] = 8 # update H to COOH
            position = look_for_particle(surround(lattice, i, j, d), 9) # look for CO2 srf
            i_changed, j_changed = pos_update(position, i, j, 0)   # update to empty
        all_lattice_populations[5]+=1
        all_lattice_populations[6]-=1
        all_lattice_populations[3]-=1             
    
    elif process_number == 14:
        process_count[21]+=1
        #CO Diff
        if species==1:
            lattice[i][j]=0
            no_lat_interactions=[]
            for suitable_site in suitableSurroundings(lattice,i,j,0):
                if 1 not in superSurround(lattice,suitable_site[0],suitable_site[1],d):
                    no_lat_interactions.append(suitable_site)            
            site_selection = no_lat_interactions[random.randint(0,len(no_lat_interactions)-1)]        
            lattice[site_selection[0]][site_selection[1]] = 1
            i_changed, j_changed = site_selection[0], site_selection[1]
        elif species==0:
            lattice[i][j]=1
            position = look_for_particle(surround(lattice, i, j, d), 1)
            i_changed, j_changed = pos_update(position, i, j, 0)                

    elif process_number == 15:
        process_count[20]+=1
        #OH Diff
        if species==3:
            lattice[i][j]=0
            position = look_for_particle(surround(lattice, i, j, d), 0)
            i_changed, j_changed = pos_update(position, i, j, 3)
        elif species==0:
            lattice[i][j]=3
            position = look_for_particle(surround(lattice, i, j, d), 3)
            i_changed, j_changed = pos_update(position, i, j, 0)            
            

    elif process_number == 16:  
        process_count[19]+=1
        #H Diff
        if species==4:
            lattice[i][j]=0
            position = look_for_particle(surround(lattice, i, j, d), 0)
            i_changed, j_changed = pos_update(position, i, j, 4)
        elif species==0:
            lattice[i][j]=4
            position = look_for_particle(surround(lattice, i, j, d), 4)
            i_changed, j_changed = pos_update(position, i, j, 0)              

    elif process_number == 17:   
        process_count[18]+=1
        #O Diff
        if species==5:
            lattice[i][j]=0
            position = look_for_particle(surround(lattice, i, j, d), 0)
            i_changed, j_changed = pos_update(position, i, j, 5)
        elif species==0:
            lattice[i][j]=5
            position = look_for_particle(surround(lattice, i, j, d), 5)
            i_changed, j_changed = pos_update(position, i, j, 0)              
    
    changed_i.append(i)
    changed_i.append(i_changed)
    changed_j.append(j)
    changed_j.append(j_changed)

    return all_lattice_populations,lattice,changed_i,changed_j,number_of_H2,number_of_CO2 #If nothing happened, lattice does not change
#####################################################     

#####################################################
#Look at the particles around and return them in a list
def surround(lattice, i, j, d):
    
    if i == 0:
        above = lattice[d - 1][j]
        below = lattice[i + 1][j]
    elif i == d - 1:
        above = lattice[i - 1][j]
        below = lattice[0][j]
    else:
        above = lattice[i - 1][j]
        below = lattice[i + 1][j]
    if j == 0:
        left = lattice[i][d - 1]
        right = lattice[i][j + 1]
    elif j == d - 1:
        right = lattice[i][0]
        left = lattice[i][j - 1]
    else:
        left = lattice[i][j - 1]
        right = lattice[i][j + 1]

    surroundlist = [above, below, left, right]
    return surroundlist

#####################################################
#Look at the particles around, including diagonal
def superSurround(lattice, i, j, d):
    
    if i == 0:
        above_pos=d-1
        below_pos=i+1
    elif i == d-1:
        above_pos =i-1
        below_pos=0
    else:
        above_pos =i-1
        below_pos=i+1 
        
    if j == 0:
        left_pos=d-1
        right_pos=j+1
    elif j == d-1:
        left_pos =j-1
        right_pos=0
    else:
        left_pos =j-1
        right_pos=j+1  
    
    above_left=lattice[above_pos][left_pos]
    above=lattice[above_pos][j]
    above_right=lattice[above_pos][right_pos]
    
    left=lattice[i][left_pos]
    right=lattice[i][right_pos]
    
    below_left=lattice[below_pos][left_pos]
    below=lattice[below_pos][j]
    below_right=lattice[below_pos][right_pos]

    super_surroundlist = [above_left, above, above_right, left, right, below_left, below, below_right]
    return super_surroundlist
######################################################

#####################################################
def particle_position_list(species):
    position_list = []
    for i in range(d):
        for j in range(d):
            if lattice[i][j] == species:
                position_list.append([i,j])
    return position_list
########################################################

#####################################################
#Takes care of the adosrption
def adsorption(time_elapsed, CO_adsorb_rate_constant, H2O_adsorb_rate_constant, CO2_adsorb_rate_constant, lattice, time_since_last_CO, time_since_last_H2O, time_since_last_CO2, number_of_CO, number_of_H2O, number_of_CO2):
 

    changed_i=[]
    changed_j=[]   
             
    if number_of_CO !=0 and CO_adsorb_rate_constant !=0:
        if time_elapsed - time_since_last_CO >= 1 / CO_adsorb_rate_constant: #if suffcient time has passed, evalute the sites for CO adsorption
            list_of_CO_suitable_sites=[] #free sites not next to CO, O or H
            for i in range(d):
                for j in range(d):
                    if lattice[i][j] == 0:
                        super_surround_list = superSurround(lattice, i, j, d)
                        if 1 not in super_surround_list and 3 not in super_surround_list and 5 not in super_surround_list:
                            list_of_CO_suitable_sites.append([i,j])
            if len(list_of_CO_suitable_sites) != 0:
                r = random.randint(0, len(list_of_CO_suitable_sites)-1) #randomly pick one of the suitable sites.
                site = list_of_CO_suitable_sites[r]
                site_i = site[0]
                site_j = site[1]
                lattice[site_i][site_j] = 1
                time_since_last_CO = time_elapsed
                if replenish_reactant_gas==False:
                    number_of_CO-=1
                changed_i.append(site_i)
                changed_j.append(site_j)
                all_lattice_populations[0]+=1
                process_count[14]+=1
                    
    if number_of_H2O != 0 and H2O_adsorb_rate_constant !=0:
        if time_elapsed - time_since_last_H2O >= 1 / H2O_adsorb_rate_constant:
            list_of_H2O_suitable_sites=[]
            for i in range(d):
                for j in range(d):
                    if lattice[i][j] == 0:
                        super_surround_list = superSurround(lattice, i, j, d)
                        if 5 not in super_surround_list:
                            list_of_H2O_suitable_sites.append([i,j])
            if len(list_of_H2O_suitable_sites) != 0:
                r = random.randint(0, len(list_of_H2O_suitable_sites)-1)
                site = list_of_H2O_suitable_sites[r]
                site_i = site[0]
                site_j = site[1]
                lattice[site_i][site_j] = 2
                time_since_last_H2O = time_elapsed
                if replenish_reactant_gas==False:
                    number_of_H2O-=1
                changed_i.append(site_i)
                changed_j.append(site_j)
                all_lattice_populations[1]+=1
                process_count[16]+=1
                    
    if number_of_CO2 != 0 and CO2_adsorb_rate_constant !=0:
        if supress_reverse_product_reaction==False: #If we allow this product gas to undergo a reverse adsorption...
            if time_elapsed - time_since_last_CO2 >= 1 / CO2_adsorb_rate_constant:
                list_of_CO2_suitable_sites=[] #free sites not nect to CO, O or H
                for i in range(d):
                    for j in range(d):
                        if lattice[i][j] == 0:
                            list_of_CO2_suitable_sites.append([i,j]) #No lateral interactions with CO2
                if len(list_of_CO2_suitable_sites) != 0:
                    r = random.randint(0, len(list_of_CO2_suitable_sites)-1)
                    site = list_of_CO2_suitable_sites[r]
                    site_i = site[0]
                    site_j = site[1]
                    lattice[site_i][site_j] = 9
                    time_since_last_CO2 = time_elapsed
                    number_of_CO2-=1
                changed_i.append(site_i)
                changed_j.append(site_j)
                all_lattice_populations[6]+=1
                process_count[23]+=1
                
    return time_since_last_CO, time_since_last_H2O, time_since_last_CO2, lattice, number_of_CO, number_of_H2O, number_of_CO2, changed_i, changed_j
#####################################################

#########################################################
def desorption(time_elapsed, lattice, time_since_last_desorb_CO, time_since_last_desorb_H2O, time_since_last_desorb_CO2, number_of_CO, number_of_H2O, number_of_CO2):
    
    changed_i=[]
    changed_j=[]
    
    if time_elapsed - time_since_last_desorb_CO >= 1 / CO_desorb_rate_constant :
        #### need a list of CO positions
        CO_pos_list = particle_position_list(1)
        if len(CO_pos_list) !=0:
            r = random.randint(0, len(CO_pos_list)-1)
            position_co_to_desorb = CO_pos_list[r]
            lattice[position_co_to_desorb[0]][position_co_to_desorb[1]] = 0
            time_since_last_desorb_CO = time_elapsed
            if replenish_reactant_gas==False:
                number_of_CO+=1
            changed_i.append(position_co_to_desorb[0])
            changed_j.append(position_co_to_desorb[1])
            all_lattice_populations[0]-=1
            process_count[15]+=1
                
    if time_elapsed - time_since_last_desorb_H2O >= 1 / H2O_desorb_rate_constant:
        #### need a list of H2O positions
        H2O_pos_list = particle_position_list(2)
        if len(H2O_pos_list) != 0:
            r = random.randint(0, len(H2O_pos_list)-1)
            position_H2O_to_desorb = H2O_pos_list[r]
            lattice[position_H2O_to_desorb[0]][position_H2O_to_desorb[1]] = 0
            time_since_last_desorb_H2O = time_elapsed
            
            if replenish_reactant_gas==False:
                number_of_H2O+=1
            changed_i.append(position_H2O_to_desorb[0])
            changed_j.append(position_H2O_to_desorb[1])
            all_lattice_populations[1]-=1
            process_count[17]+=1

    #Unlike adsorption, CO2 is always allowed to desorb, as this is not a reverse reaction of this product gas
    if time_elapsed - time_since_last_desorb_CO2 >= 1 / CO2_desorb_rate_constant :
        #### need a list of CO2 positions
        CO2_pos_list = particle_position_list(9)
        if len(CO2_pos_list) !=0:
            r = random.randint(0, len(CO2_pos_list)-1)
            position_co2_to_desorb = CO2_pos_list[r]
            lattice[position_co2_to_desorb[0]][position_co2_to_desorb[1]] = 0
            time_since_last_desorb_CO2 = time_elapsed
            number_of_CO2+=1
            changed_i.append(position_co2_to_desorb[0])
            changed_j.append(position_co2_to_desorb[1]) 
            all_lattice_populations[6]-=1
            process_count[22]+=1
                
    return time_since_last_desorb_CO, time_since_last_desorb_H2O, time_since_last_desorb_CO2, lattice, number_of_CO, number_of_H2O, number_of_CO2, changed_i, changed_j
#################################################################

#####################################################################
def check_reversibility(reaction_number, forwards):
    
    reaction_number = int(reaction_number)
    
    reaction_history[(reaction_number)] +=1
    if forwards == 1:
        execution_history[(reaction_number)] += 1
    else:
        execution_history[(reaction_number)] -= 1
########################################################################

#######################################################################   
#Determines which reaction channel the selected process belongs to, for stiffness scaling purposes. 
def check_reaction_number(i):
    if i<14:
        if i % 2 == 0:
            
            reaction_number = i / 2
            forwards = 1
        else:
            reaction_number = (i-1) / 2
            forwards = 0
        return int(reaction_number), forwards
    else:
        if i == 14:
            reaction_number = 7
        if i == 15:
            reaction_number = 8
        if i == 16:
            reaction_number = 9
        if i == 17:
            reaction_number = 10
        
        return reaction_number, 1
########################################################################
                

##########################################################################
#Counts the number of each species on the surface
def speciesCount(species):
    #every 50 ms, calculate the coverages and print them to a file
    count = 0
    for i in range(d):
        for j in range(d):
            if lattice[i][j] == species:
                count += 1
    #coverage = count / (d*d)
    
    return count
############################################################################

sweep = 0

temp_run = 0
run = 0

while temp_run< number_of_temperatures:
    print("")
    print("Temperature = ", temperature)
    run = 0
    if temperature > 649: #increase the strength of the stiffness scaling at higher temperatures to speed up the simulation
        scaling_factor = scaling_factor * 10

    while run < number_of_runs:
        imades_saved = 0 
        lattice = [[0 for i in range(d)] for j in range(d)] #Starts off completely empty
    
        #Put CO at one corner and H2O at the furthest point away to get things going but get as close to starting empty as possible (without actually doing so, as this would mean nothing could ever happen)
        lattice[0][0]=1
        lattice[round(d/2)][round(d/2)]=2
        
        #List of activation barriers for all reaction and diffusion steps (in the same order as before)
        activation_barrier_matrix=[1.320,0.240,0.880,1.410,0.320,0.350,0.810,0.018,0.690,0.439,0.720,1.471,0.540,0.580,0.09,0.4,0.017,0]
        #For entry 7 (reverse reaction of H2) and 9 (reverse of CO2) which require these species straight from the gas phase, the rate constants should be treated like those for adsorption in terms of dependence on number of gas molecules
    
        #Rate constant list now varies with temperature
        rate_constant_matrix = [pre_exp_factor * np.exp(-(activation_barrier_matrix[index]*e)/(kb*temperature)) for index in range(len(activation_barrier_matrix))]
        if supress_reverse_product_reaction==True:
            rate_constant_matrix[7]=0
            rate_constant_matrix[9]=0
        else:
            #Reverse H2 reaction
            rate_constant_matrix[7] = (number_of_H2/(number_of_sites*height_of_gas)) * np.sqrt((kb*temperature)/(2*np.pi*H2_mass)) * np.exp(-activation_barrier_matrix[7]*e/(kb*temperature))
            #Reverse CO2 reaction
            rate_constant_matrix[9] = (number_of_CO2/(number_of_sites*height_of_gas)) * np.sqrt((kb*temperature)/(2*np.pi*CO2_mass)) * np.exp(-activation_barrier_matrix[9]*e/(kb*temperature))
        
        #Counts the number of times each process occurs
        #Key: regular reactions, CO and H2O ads/des, diffusion, CO2 des/ads (same order as in Chutia paper)
        process_count=[0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,  0,0, 0,0,  0,0,0,0,  0,0] 
        #This keeps track of the number of each species on the lattice over time for use in tracking coverage. I figured this would be more efficient than needing to run through the whole surface each sweep for each species and save a lot of crying and screaming (from Nimbus as well as us)
        all_lattice_populations=[speciesCount(1),speciesCount(2),0,0,0,0,0] #Number of adsorbed: CO, H2O, OH, H, O, COOH, CO2(srf)
        all_gas_amounts=[number_of_CO, number_of_H2O, number_of_H2, number_of_CO2]
        
        number_of_lattice_species=len(all_lattice_populations)
        number_of_gas_species=len(all_gas_amounts)
        
        unequilibriated = [1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0] # all reactions start off unequilibriated, diffusion is always equilibriated
        execution_history = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        reaction_history = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        sufficiently_executed = [0,0,0,0,0,0,0,1,1,1,1]
        
        time_elapsed = 0
        
        all_time_elapsed = [0]
        
        #Adsorption barriers for each gas species
        CO_adsorb_barrier = 0*e #Specified to be 0
        H2O_adsorb_barrier = 0*e #Specified to be 0
        CO2_adsorb_barrier = 0*e #This is assumed based on the other adsorption barriers
        
        #Rate constants for each gas, with the latter varying with the number of gas molecules of each present
        CO_adsorb_rate_constant = (number_of_CO/((number_of_sites)*height_of_gas)) * np.sqrt((kb*temperature)/(2*np.pi*CO_mass)) * np.exp(-CO_adsorb_barrier/(kb*temperature))
        H2O_adsorb_rate_constant = (number_of_H2O/((number_of_sites)*height_of_gas)) * np.sqrt((kb*temperature)/(2*np.pi*H2O_mass)) * np.exp(-H2O_adsorb_barrier/(kb*temperature))
        CO2_adsorb_rate_constant = (number_of_CO2/((number_of_sites)*height_of_gas)) * np.sqrt((kb*temperature)/(2*np.pi*CO2_mass)) * np.exp(-CO2_adsorb_barrier/(kb*temperature)) #No CO2 to begin with, so it cannot adsorb (rate const is 0 as number of CO2 is 0). This will always stay at 0 if supress_reverse_product_reaction is True
        
        #Most favourable adsorption energies
        CO_adsorb_energy = -1.924*e #bridging (B) adsorption
        H2O_adsorb_energy = -0.147*e #bridging as well
        CO2_adsorb_energy = 0*e
        #I am assuming desorption is straight from the most favoured adsorption mode (bridging for both CO, H2O)rather than going through another step (such as changing to a less favoured mode)
        CO_desorb_barrier = -CO_adsorb_energy
        H2O_desorb_barrier = -H2O_adsorb_energy
        CO2_desorb_barrier = -CO2_adsorb_energy #Barrier specified to be 0
        
        #Gas desorption rate constants (CO, H2O are reactants; CO2 is a product)
        CO_desorb_rate_constant = pre_exp_factor*np.exp(-CO_desorb_barrier/(kb*temperature))
        H2O_desorb_rate_constant = pre_exp_factor*np.exp(-H2O_desorb_barrier/(kb*temperature))
        CO2_desorb_rate_constant = pre_exp_factor*np.exp(-CO2_desorb_barrier/(kb*temperature))
        
        
        #Simulation times at which an adsorption of each species took place
        time_since_last_CO = 0
        time_since_last_H2O = 0
        time_since_last_CO2 = 0
        
        number_of_CO = round((starting_CO_fraction*starting_pressure*(d**2)*(Pd_atom_spacing**2)*height_of_gas)/(kb*temperature))
        number_of_H2O = round((starting_H2O_fraction*starting_pressure*(d**2)*(Pd_atom_spacing**2)*height_of_gas)/(kb*temperature))
        number_of_CO2 = 0
        number_of_H2 = 0
        #Simulation time at which a desorption last took place
        time_since_last_desorb_CO = 0
        time_since_last_desorb_H2O = 0
        time_since_last_desorb_CO2 = 0
        

        time_step = 0
        
        reaction_count = [0,0,0,0,0,0,0,0,0,0,0]
        number_of_times_scaled = [0,0,0,0,0,0,0,0,0,0,0]
        last_scale_sweep = 0
        unscale_count = 0
        last_scale_sweep = [0,0,0,0,0,0,0,0,0,0,0]
        Co_coverage_list = []
        time_of_coverage = []
        
        #####################################
        #INITIAL PERIODIC AVERAGING STUFF
        averaged_time=[0] #Reduced list of time values
        
        averaged_coverages=[[all_lattice_populations[coverage_index]/number_of_sites for coverage_index in range(number_of_lattice_species)]]
        averaged_gas_amounts=[[all_gas_amounts[gas_species_index] for gas_species_index in range(number_of_gas_species)]]
        
        
        previous_current_interval=0
        current_period_times=[] #All time values \\ \\ \\ \\
        current_period_coverages=[[] for tracked_species in all_lattice_populations] #All coverage points lying between current bounds to be averaged once the next one is crossed
        current_period_gas_amounts=[[],[],[],[]] #All numbers of gas molecules (CO, H2O, H2, CO2)
        
        #####################################
        
        all_changed_i=[]
        all_changed_j=[]
        
        
        #############################################
        #############################################
        #THE ACTUAL SIMULATION RUN
        
        while time_elapsed < max_simulation_time:
            
            #initialise a bunch of arrays
            rate_constant_for_time = []
        
            reactions_that_happen_in_sweep = []
            position_of_processes_that_can_happen = []
            
            rates_of_processes_that_can_happen = []
            processes_that_can_happen_in_sweep = []
            position_of_processes_that_can_happen_in_sweep = []
            
            all_changed_i=[]
            all_changed_j=[]
            
            for i in range(d):
                for j in range(d):
                    #get a list of particles around the current lattice point
                    particles_around_lattice_point = surround(lattice, i, j, d)
                    processes_that_can_happen = []
                    
                    # generate a list of the processes that can happen given what is around, and their rates
                    processes_that_can_happen,  position_of_processes_that_can_happen = whichProcess(lattice[i][j], particles_around_lattice_point, number_of_CO2, number_of_H2, i, j, position_of_processes_that_can_happen )
                    
                    for process in processes_that_can_happen:
                        position_of_processes_that_can_happen_in_sweep.append([i,j])
                        processes_that_can_happen_in_sweep.append(process)
                        rates_of_processes_that_can_happen.append(rate_constant_matrix[process])
        
                    ##########################
                    # At this point we have a list of all the processes that can happen in this configuration, and their positions
                    # Now at the end of the lattice, use KMC to pick one
                    
            for m in rate_constant_matrix:
                for entry in processes_that_can_happen_in_sweep: #checks the reaction number of the posible processes
                    reaction_number, forwards = check_reaction_number(entry)
            
            changed_i, changed_j = [], []
            if len(rates_of_processes_that_can_happen) != 0:
                #Determines which process happens and its position on the lattice
                selected_process_list_position = KMC_Sweep(rates_of_processes_that_can_happen, position_of_processes_that_can_happen)
                selected_process = processes_that_can_happen_in_sweep[selected_process_list_position]
                position_of_selected_process = position_of_processes_that_can_happen_in_sweep[selected_process_list_position]
                
                #Determines the time to increase by
                timestep_rng = random.uniform(0, 1)
                if sum(rates_of_processes_that_can_happen) != 0:
                        time_step = -(1/ sum(rates_of_processes_that_can_happen))*np.log(timestep_rng) 
                time_elapsed += time_step
                  
                reaction_number, forwards = check_reaction_number(selected_process) #determines which reaction happened and in which direction (forwards or backwards)
               
                all_lattice_populations,lattice,changed_i,changed_j,number_of_H2,number_of_CO2 = siteChange(all_lattice_populations,selected_process, processes_that_can_happen, i, j , lattice[i][j], particles_around_lattice_point, number_of_CO2, number_of_H2, position_of_selected_process)                
                
                ##########################
                # stiffness scaling (oh boy)
                # reactions have happened now (yay)
                
                #add one to the reaction counter, which counts total forwards and backwards
                reaction_count[reaction_number] += 1
                
                #work out number of forwards reactions - number of backwards
                check_reversibility(reaction_number, forwards) #updates execution history
                
                #Stiffness scaling parameters
                rescale_step = 5
                check_equilibraited_step = 5
                sufficently_executed_threshold = 5
                
                reaction_number = int(reaction_number)
                
                if unequilibriated[reaction_number] == 1:
                
                    #then an unequilibriated reaction has happened, so unscale
                    unscale_count +=1
                    unequilibriated = [1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0]
                    execution_history = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                    reaction_count = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                    sufficiently_executed = [0,0,0,0,0,0,0,1,1,1,1]
                
                    
                    rate_constant_matrix = [pre_exp_factor * np.exp(-(activation_barrier*e)/(kb*temperature)) for activation_barrier in activation_barrier_matrix]
                    if supress_reverse_product_reaction==True:
                        rate_constant_matrix[7]=0
                        rate_constant_matrix[9]=0
                    else:
                        #Reverse H2 reaction
                        rate_constant_matrix[7] = (number_of_H2/(number_of_sites*height_of_gas)) * np.sqrt((kb*temperature)/(2*np.pi*H2_mass)) * np.exp(-activation_barrier_matrix[7]*e/(kb*temperature))
                        #Reverse CO2 reaction
                        rate_constant_matrix[9] = (number_of_CO2/(number_of_sites*height_of_gas)) * np.sqrt((kb*temperature)/(2*np.pi*CO2_mass)) * np.exp(-activation_barrier_matrix[9]*e/(kb*temperature))
                
                ################################################
                #check if quasi equilibriated
                if reaction_count[reaction_number] > check_equilibraited_step and execution_history[reaction_number] < delta * reaction_count[reaction_number] or reaction_number > 6:
                    # the above indicates that the forwards and backwards happen at near enough the same frequency, and allows diffsuion to always be scaled
                    unequilibriated[reaction_number] = 0 # flag this reaction as QE
                    # now need to check whether this reaction is sufficently executed
                    # this is because for low reaction counts, the above will always run
                    if reaction_count[reaction_number] > sufficently_executed_threshold:
                        # flag as sufficently executed
                        sufficiently_executed[reaction_number] = 1
                
                else: # this means that the forwards and back arent happening at the same rate, so not QE. Hence unscale
                    #flag as uneq and not sufficently exec
                    unequilibriated[reaction_number] = 1
                    sufficiently_executed[reaction_number] = 1
                
                #Scale the relevant reaction rates
                if (sweep % rescale_step) == 0:
                   # for k in range(len(rate_constant_matrix)):
                        k = selected_process
                        scale_reaction_number = reaction_number
                        if unequilibriated[scale_reaction_number] == 0 and sufficiently_executed[scale_reaction_number] == 1 or scale_reaction_number>6 :
                            #scale the rates if they are uneq and sufficently exec
                            number_of_times_scaled[reaction_number] +=1
                            if scale_reaction_number < 7:
                                              
                                if forwards == 1:
                                    rate_constant_matrix[k] = rate_constant_matrix[k] / scaling_factor
                                    rate_constant_matrix[k+1] = rate_constant_matrix[k + 1] / scaling_factor
                                else:
                                    rate_constant_matrix[k] = rate_constant_matrix[k] / scaling_factor
                                    rate_constant_matrix[k-1] = rate_constant_matrix[k - 1] / scaling_factor
                        
                            if scale_reaction_number > 6:
                                rate_constant_matrix[k] = rate_constant_matrix[k] / scaling_factor
            
            #If no reactions can occur (length of list_of_processes_that_can_happen is 0), likely due to the lattice being full, this should allow the time to continue by the last amount such that desorption will eventually happen and things can continue
            else:
                time_elapsed += time_step
        
            if changed_i!=[]:
                for single_changed_i in changed_i:
                    all_changed_i.append(single_changed_i)
            if changed_j!=[]:
                for single_changed_j in changed_j:
                    all_changed_j.append(single_changed_j)
        
           
        
        #ALL DESORPTION HAPPENS HERE
            time_since_last_desorb_CO, time_since_last_desorb_H2O, time_since_last_desorb_CO2, lattice, number_of_CO, number_of_H2O, number_of_CO2, changed_i, changed_j  = desorption(time_elapsed, lattice, time_since_last_desorb_CO, time_since_last_desorb_H2O, time_since_last_desorb_CO2, number_of_CO, number_of_H2O, number_of_CO2)
           
            if replenish_reactant_gas==False: #If we are even modelling changing reactant gas quantities...
                if time_since_last_desorb_CO==time_elapsed: #If CO desorbed (which it didn't)...
                    #Calculate the new adsorption rate constant for CO now that the number of molecules in the gas has changed
                    CO_adsorb_rate_constant = (number_of_CO/((d**2)*height_of_gas)) * np.sqrt((kb*temperature)/(2*np.pi*CO_mass)) * np.exp(-CO_adsorb_barrier/(kb*temperature))
                if time_since_last_desorb_H2O==time_elapsed:    
                    H2O_adsorb_rate_constant = (number_of_H2O/((d**2)*height_of_gas)) * np.sqrt((kb*temperature)/(2*np.pi*H2O_mass)) * np.exp(-H2O_adsorb_barrier/(kb*temperature))
        
            if supress_reverse_product_reaction==False: #If we let product gases undergo reverse reactions, including allowing CO2 to adsorb onto the lattice...
                if time_since_last_desorb_CO2==time_elapsed: #If CO2 desorbed...
                    #Calculate the new adsorption rate constant for CO2 now that the number of molecules in the gas has changed
                    CO2_adsorb_rate_constant = (number_of_CO2/((d**2)*height_of_gas)) * np.sqrt((kb*temperature)/(2*np.pi*CO2_mass)) * np.exp(-CO2_adsorb_barrier/(kb*temperature))
        
            if changed_i!=[]:
                for single_changed_i in changed_i:
                    all_changed_i.append(single_changed_i)
            if changed_j!=[]:
                for single_changed_j in changed_j:
                    all_changed_j.append(single_changed_j)
            
        #ALL ADSORPTION HAPPENS HERE
            time_since_last_CO, time_since_last_H2O, time_since_last_CO2, lattice, number_of_CO, number_of_H2O, number_of_CO2, changed_i, changed_j = adsorption(time_elapsed, CO_adsorb_rate_constant, H2O_adsorb_rate_constant, CO2_adsorb_rate_constant, lattice, time_since_last_CO, time_since_last_H2O, time_since_last_CO2, number_of_CO, number_of_H2O, number_of_CO2)
            
            if replenish_reactant_gas==False: #If we are even modelling changing reactant gas quantities...    
                if time_since_last_CO==time_elapsed: #If CO adsorbed...
                    #Calculate the new adsorption rate constant for CO now that the number of molecules in the gas has changed
                    CO_adsorb_rate_constant = (number_of_CO/((d**2)*height_of_gas)) * np.sqrt((kb*temperature)/(2*np.pi*CO_mass)) * np.exp(-CO_adsorb_barrier/(kb*temperature))
                if time_since_last_H2O==time_elapsed:    
                    H2O_adsorb_rate_constant = (number_of_H2O/((d**2)*height_of_gas)) * np.sqrt((kb*temperature)/(2*np.pi*H2O_mass)) * np.exp(-H2O_adsorb_barrier/(kb*temperature))        
            
            if supress_reverse_product_reaction==False: #If we let product gases undergo reverse reactions, including allowing CO2 to adsorb onto the lattice...
                if time_since_last_CO2==time_elapsed: #If CO2 adsorbed...
                    #Calculate the new adsorption rate constant for CO2 now that the number of molecules in the gas has changed
                    CO2_adsorb_rate_constant = (number_of_CO2/((d**2)*height_of_gas)) * np.sqrt((kb*temperature)/(2*np.pi*CO2_mass)) * np.exp(-CO2_adsorb_barrier/(kb*temperature))
        
            if changed_i!=[]:
                for single_changed_i in changed_i:
                    all_changed_i.append(single_changed_i)
            if changed_j!=[]:
                for single_changed_j in changed_j:
                    all_changed_j.append(single_changed_j)
        
            
           
            
            #Change all gas species totals to the new values
            all_gas_amounts=[number_of_CO, number_of_H2O, number_of_H2, number_of_CO2]
            
            #Update reverse product reactions if appropriate
            if supress_reverse_product_reaction==False:
                #Reverse H2 reaction
                rate_constant_matrix[7] = (number_of_H2/(number_of_sites*height_of_gas)) * np.sqrt((kb*temperature)/(2*np.pi*H2_mass)) * np.exp(-activation_barrier_matrix[7]*e/(kb*temperature))
                #Reverse CO2 reaction
                rate_constant_matrix[9] = (number_of_CO2/(number_of_sites*height_of_gas)) * np.sqrt((kb*temperature)/(2*np.pi*CO2_mass)) * np.exp(-activation_barrier_matrix[9]*e/(kb*temperature))
            
            #PERIODIC AVERAGING OF DATA OCCURRING DURING SIMULATION:
            
            current_interval=int(time_elapsed/period) #Current time interval between bounds we are currently in
            if current_interval>previous_current_interval:
                #print("New averaged data")
                averaged_time.append(np.mean(current_period_times)) #Append average time
                averaged_coverages.append([np.mean(current_period_coverages[coverage_index]) for coverage_index in range(number_of_lattice_species)]) #Append average coverages within this period
                averaged_gas_amounts.append([np.mean(current_period_gas_amounts[gas_species_index]) for gas_species_index in range(number_of_gas_species)]) #Append average gas populations
                
                current_period_times=[]
                current_period_coverages=[[] for tracked_species in all_lattice_populations] #Reset list of data points being averaged
                current_period_gas_amounts=[[],[],[],[]]
            
            for coverage_index in range(number_of_lattice_species): 
                current_period_coverages[coverage_index].append( all_lattice_populations[coverage_index]/number_of_sites ) #Adds current data point to the list of raw data points within this time range to be averaged when a bound is crossed
            for gas_species_index in range(number_of_gas_species):
                current_period_gas_amounts[gas_species_index].append( all_gas_amounts[gas_species_index])
            current_period_times.append(time_elapsed)  
        
            previous_current_interval=current_interval   
            
           
            if save_image==True:
                if images_saved < len(image_saving_times):
                    if time_elapsed >= image_saving_times[images_saved]:
                        print("At time ",image_saving_times[images_saved],"s:")
                        latticePlot(lattice, all_changed_i, all_changed_j,False,True)
                        images_saved +=1
            
            
            sweep +=1
            
        ###################
        #Plot the lattice at the given time
        latticePlot(lattice, all_changed_i, all_changed_j,False,True)
        
        #Does one more data and time averaging in the final period (that has not been completed so cannot be done in the for loop)
        averaged_coverages.append([np.mean(all_lattice_populations[coverage_index])/number_of_sites for coverage_index in range(number_of_lattice_species)]) #Append average coverage within this period
        averaged_time.append(np.mean(current_period_times)) #Append average time        
        averaged_gas_amounts.append([np.mean(current_period_gas_amounts[gas_species_index]) for gas_species_index in range(number_of_gas_species)]) #Append average gas populations
    
        print("end of run ", run)
        
        #Saving a bunch of data to well named text files
        f = open(f"T={temperature}K,P={np.round(starting_pressure/1e5,2)}bar,x(CO)={starting_CO_fraction},run={repeat_run} time_values.txt", "a")
        f.writelines( "Time" "\n")
        for averaged_points in range(len(averaged_time)):
            f.writelines(str(averaged_time[averaged_points]) +"\n")
        f.writelines("\n")
        f.close()
    
        f = open(f"T={temperature}K,P={np.round(starting_pressure/1e5,2)}bar,x(CO)={starting_CO_fraction},run={repeat_run} coverage.txt", "a")    
        f.writelines("CO" + " " + "H2O" + " " + "OH" + " " + "H" + " " + "O" + " " + "COOH" + " " + "CO2" "\n")
        for averaged_points in range(len(averaged_time)):
            for coverage_index in range(number_of_lattice_species):
                f.writelines(str(round(averaged_coverages[averaged_points][coverage_index],5)) + " ")
            f.writelines("\n")
        f.writelines("\n")
        f.close()
    
        f = open(f"T={temperature}K,P={np.round(starting_pressure/1e5,2)}bar,x(CO)={starting_CO_fraction},run={repeat_run} gases.txt", "a")    
        f.writelines("CO" + " " + "H2O" + " " + "H2" + " " + "CO2" "\n")
        for averaged_points in range(len(averaged_time)):
            for gas_species_index in range(number_of_gas_species):
                f.writelines(str(averaged_gas_amounts[averaged_points][gas_species_index]) + " ")
            f.writelines("\n")
        f.writelines("\n")
        f.close()
    
        f=open(f"T={temperature}K,P={np.round(starting_pressure/1e5,2)}bar,x(CO)={starting_CO_fraction},run={repeat_run} processes.txt", "a")
        for selected_process_count in process_count:
            f.writelines(str(selected_process_count) + " " )
        f.writelines("\n")
        f.close()
        print(time_elapsed)
        run +=1
    temperature += temperature_increment
    temp_run +=1
    
############################################################
############################################################    

#Does one more data and time averaging in the final period (that has not been completed so cannot be done in the for loop)
averaged_coverages.append([np.mean(all_lattice_populations[coverage_index])/number_of_sites for coverage_index in range(number_of_lattice_species)]) #Append average coverage within this period
averaged_time.append(np.mean(current_period_times)) #Append average time        
averaged_gas_amounts.append([np.mean(current_period_gas_amounts[gas_species_index]) for gas_species_index in range(number_of_gas_species)]) #Append average gas populations

#Tells you how many sweeps of the lattice occured, just for interest
print("number of sweeps ", sweep)

