import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy import signal

pi = 3.14159265359
u = .045 # units P; blood viscosity ranges from 3.5 to 5.5 cp
R_c = 0.0004 # units cm; radius of the capillary
D_n = .00023 # units cm; diameter of nicotine
R_n = D_n/2 # units cm; radius of nicotine
f = 6*pi*u*R_n
K = 1.38065*10**(-16) # Boltzman Constant; g*cm^2*K^-1*s^-2
T = 310 # units K; body temp assumed 37 degrees celsius

Dij = (K*T)/f
print(Dij)

#Short contact time validation
v = 1.46/10 #cm/s; average velocity in the capillaries
v_max = 2*v #units cm/s; assuming pousielle flow
z_cap = .06 #cm; length of the capillary
short_contact = (z_cap*Dij) / (v_max*(R_c**2))
print(short_contact)

#Flux using short contact time solution
L = R_c*2 #cm; diameter of the capillary
Pe = L*v / Dij 

# Define Flux equation
def calculate_flux(Dij, Co, R_c, Pe, z_cap):
    return (((0.6783 * Dij * Co) / R_c) * (Pe * (R_c / z_cap))**(1/3))

# Create an array of Co values
c1 = 0 #mol/mL; min concentration of nicotine in alveoli
c2 = (1.12*10**-6)/1000 #mol/mL; max concentration of nicotine in alveoli
Co_values = np.linspace(c1, c2, 100)  # Adjust the range and number of points as needed

# Calculate flux for each Co value
flux_values = [calculate_flux(Dij, Co, R_c, Pe, z_cap) for Co in Co_values]

# Plot Co vs. Flux
plt.plot(Co_values, flux_values, color = 'c')
plt.xlabel('Initial Concentration (mol/mL)')
plt.ylabel('Flux (mol/cm^2*s)')
plt.title('Flux vs. Initial Concentration')
plt.grid(True)
plt.show()

#Calculate Uptake; uptake is the integral of Area*Flux with respect to time
# Define uptake equation
def calculate_uptake(Co1, t):
    return (((((0.6783 * Dij * Co1) / R_c) * (Pe * (R_c / z_cap))**(1/3))*t)*A)

A = 100*100 #cm; average SA of alveoli in human lung
t_values = np.linspace(1.0, 480.0, 480)
t = 480 #seconds; 8 minutes
Co1_values = np.linspace(c1, c2, 100) #change these values to be more relevant... concentration should decrease so flux should decrease as initial Co decreases

# Create a meshgrid for the variables
concentration_mesh, time_mesh = np.meshgrid(
    Co1_values, t_values, indexing='ij'
)

uptake_co_values = [calculate_uptake(Co1, t) for Co1 in Co_values]

# Plot Co vs. Flux
plt.plot(Co_values, uptake_co_values, color = 'm')
plt.xlabel('Initial Concentration (mol/mL)')
plt.ylabel('Uptake (mol)')
plt.title('Uptake vs. Initial Concentration')
plt.grid(True)
plt.show()

# Calculate uptake for each combination of variables
uptake_values = calculate_uptake(concentration_mesh, time_mesh)

# Create a 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot Uptake as a function of Time and Concentration
surface = ax.plot_surface(concentration_mesh, time_mesh, uptake_values, cmap='viridis')
ax.set_xlabel('Concentration (mol/mL)')
ax.set_ylabel('Time (s)')
ax.set_zlabel('Uptake (mol)')

# Add a color bar
fig.colorbar(surface, ax=ax, label='Uptake (mol)')

plt.title('Uptake vs. Concentration and Time')
plt.show()

NUM_ZEROS = 46000
max_molarity = 1.12*10**-6 #units mol/Liter

window = (max_molarity)*signal.windows.general_gaussian(2000, p=1.5, sig=800) #time in ms
offtime = np.zeros(NUM_ZEROS)
co_cigarette = np.concatenate([window, offtime, window, offtime, window, offtime, window, offtime, window, offtime, window, offtime, window, offtime, window, offtime, window, offtime, window, offtime, window])
time_cigarette = np.linspace(0, len(co_cigarette), len(co_cigarette))

window2 = (max_molarity)*signal.windows.general_gaussian(2000, p=1.7, sig=500) #time in ms
offtime2 = np.zeros(2000)
co_onepuff = np.concatenate([offtime2, window2, offtime2])
time_onepuff = np.linspace(0, len(co_onepuff), len(co_onepuff))

# Plot Concentration vs. Time 
plt.plot(time_cigarette, co_cigarette)
plt.xlabel('Time (ms)')
plt.ylabel('Concentration (M)')
plt.legend()
plt.title('Initial Concentration of Nicotine in Alveoli Over One Cigarette')
plt.grid(True)
plt.show()

# Plot Concentration vs. Time for One Puff
plt.plot(time_onepuff, co_onepuff)
plt.xlabel('Time (ms)')
plt.ylabel('Concentration (M)')
plt.legend()
plt.title('Initial Concentration of Nicotine in Alveoli Over One Puff')
plt.grid(True)
plt.show()

Dij_ms = Dij / 1000 #units cm^2/ms

# Calculate flux for each Co value based on puff duration and duration between puffs
flux_values_cig = [calculate_flux(Dij_ms, Co, R_c, Pe, z_cap) for Co in co_cigarette]

# Plot Flux vs. Time
plt.plot(time_cigarette, flux_values_cig, color = 'c')
plt.xlabel('Time (ms)')
plt.ylabel('Flux (mol/cm^2*s)')
plt.title('Flux vs. Time Over One Cigarette')
plt.grid(True)
plt.show()

uptake_co_values = [calculate_uptake(Co1, t) for Co1 in co_cigarette]

# Plot Co vs. Flux
plt.plot(time_cigarette, uptake_co_values, color = 'm')
plt.xlabel('Time (ms)')
plt.ylabel('Uptake (mol)')
plt.title('Uptake vs. Time Over One Cigarette')
plt.grid(True)
plt.show()

# Calculate flux for each Co value based on puff duration and duration between puffs
flux_values_onepuff = [calculate_flux(Dij_ms, Co, R_c, Pe, z_cap) for Co in co_onepuff]

# Plot Flux vs. Time
plt.plot(time_onepuff, flux_values_onepuff, color = 'c')
plt.xlabel('Time (ms)')
plt.ylabel('Flux (mol/cm^2*s)')
plt.title('Flux vs. Time Over One Puff')
plt.grid(True)
plt.show()

uptake_onepuff = [calculate_uptake(Co1, t) for Co1 in co_onepuff]

# Plot Co vs. Flux
plt.plot(time_onepuff, uptake_onepuff, color = 'm')
plt.xlabel('Time (ms)')
plt.ylabel('Uptake (mol)')
plt.title('Uptake vs. Time Over One Puff')
plt.grid(True)
plt.show()


# Create a meshgrid for the variables
concentration_mesh_onepuff, time_mesh_onepuff = np.meshgrid(
    co_onepuff, time_onepuff, indexing='ij'
)

# Calculate uptake for each combination of variables
uptake_values_onepuff = calculate_uptake(concentration_mesh_onepuff, time_mesh_onepuff)

# Create a 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot Uptake as a function of Time and Concentration
surface = ax.plot_surface(concentration_mesh_onepuff, time_mesh_onepuff, uptake_values_onepuff, cmap='viridis')
ax.set_xlabel('Concentration (mol/mL)')
ax.set_ylabel('Time (ms)')

# Add a color bars
fig.colorbar(surface, ax=ax, label='Uptake (mol)')

plt.title('Uptake vs. Concentration and Time Over One Cigarette')
plt.show()


# Surface Area Modeling

NUM_ZEROS = 46000
max_molarity = 1.12*10**-6 #units mol/Liter
max_SA = 1000000 # cm^2; 100% SA
min_SA = 510000 # cm^2; exhalation
SA_increment = (max_SA-min_SA)*1/3
max_co_onepuff = (1.12*10**-6)/1000 #mol/mL; max concentration of nicotine in alveoli

window = (max_SA)*signal.windows.general_gaussian(2000, p=1.7, sig=400) #time in ms
window1 = (SA_increment+min_SA)*signal.windows.general_gaussian(2000, p=1.7, sig=400) #time in ms
window2 = ((2*SA_increment)+min_SA)*signal.windows.general_gaussian(2000, p=1.7, sig=400) #time in ms
offtime = np.zeros(NUM_ZEROS)
SA_cigarette = np.concatenate([window, offtime, window1, offtime, window2])
time_cigarette1 = np.linspace(0, len(SA_cigarette), len(SA_cigarette))
t = 480 #s

# concentration for 3 puffs
window3 = (max_co_onepuff)*signal.windows.general_gaussian(2000, p=1.5, sig=800) #time in ms
offtime = np.zeros(NUM_ZEROS)
co_cigarette3 = np.concatenate([window3, offtime, window3, offtime, window3])
time_cigarette3 = np.linspace(0, len(co_cigarette3), len(co_cigarette3))

# Plot SA vs. Time 
plt.plot(time_cigarette1, SA_cigarette)
plt.xlabel('Time (ms)')
plt.ylabel('Surface Area (cm^2)')
plt.axhline(y=510000, color='r', linestyle='--', linewidth=2)
plt.axhline(y=1000000, color='r', linestyle='--', linewidth=2)
plt.legend()
plt.title('Surface Area of Nicotine in Alveoli vs. Time over 3 Puffs')
plt.grid(True)
plt.show()

# Uptake versus SA and constant Co
def calculate_uptake(SA, t, Co):
    return (((((0.6783 * Dij_ms * Co) / R_c) * (Pe * (R_c / z_cap))**(1/3))*t)*SA)

uptake_values3 = [calculate_uptake(SA, t, max_co_onepuff) for SA in SA_cigarette]

# Plot Co vs. Flux
plt.plot(time_cigarette1, uptake_values3)
plt.xlabel('Time (ms)')
plt.ylabel('Uptake (mol)')
plt.title('Uptake vs. Time Over 3 Puffs')
plt.grid(True)
plt.show()
