"""
Diffusion Models for LPM of the Circulation for BME 307 Final Project
Jackie Herzberg, Kishen Mitra, Maya Parekh, Darienne Rogers, Arnav Singh

These functions provide the diffusion-over-time models for use in the LPM.
Some adjustments made so that they all accept the same inputs:
- all functions have been modified to accept as inputs [t, this_conc]
- all functions have been modified such that the time input is in SECONDS
- all functions return mmol rather than mol since the volume function is in mL
"""
import numpy as np
from scipy import signal

### Nicotine patch w/ semi-infinite approach (Maya & Kishen)
# Constants
patch_amount = 14 #mg
nicotine_molar_mass = 162.23 #g/mol
thickness = 125*10**-4 #cm
patch_area = 3 * 4.2 # cm^2
patch_volume = patch_area * thickness # cm^3 (mL)
c0p = ((patch_amount) / nicotine_molar_mass) / patch_volume # mg/(g/mol)/cm^3 = mmol/mL
c1p = 0
Dij_memb_p = 1.07*10**-10 #cm^2/s

# Flux equation (units are mmol per cm^2)
def patch_flux_calc(t, c1=c1p, c0=c0p, Dij_memb=Dij_memb_p):
    if t < 10:
        return 0 # We do this bc of float errors
    else:
        flux = -1 * (Dij_memb/(np.pi*t))**(1/2) * (c1 - c0)
        return (flux)

# Mass transfer equation (Flux * area // units are: mmol/cm^2 * cm^2 = mmol)
def patch_diffusion(t, conc):
    return patch_area * patch_flux_calc(t)

# Uptake equation
def uptake_calc(conc, t):
  return (2 * area * conc * (Dij_memb * t_range / np.pi)**(1/2))

### Zyn unsteady diffusion (Arnav)
# Constants
mass_nicotine_mg = 6  # mass in mg
mass_nicotine_g = mass_nicotine_mg / 1000  # convert mg to g
volume_pouch_mm3 = 14 * 28 * 1  # volume in mm^3, this is an example value
volume_pouch_cm3 = volume_pouch_mm3 / 10e3 # convert mm3 to cm3
molecular_weight_nicotine = 162.23  # molecular weight in g/mol
C = (mass_nicotine_g / molecular_weight_nicotine) / volume_pouch_cm3 # mol/mL
L_um = 285.04 # thickness in um
L_m = L_um / 1e6 ## convert um to m
D_cm = 0.071e-7 ## Deff in cm^2/s
D_m = D_cm / 1e4 ## cm^2/s to m^2/s
D_zyn = 0.071e-9  # Diffusion coefficient in m^2/s
L = L_m   # Half-thickness of the gum in meters
C0z = C     # Initial concentration of nicotine
t = 1800    # Total time in seconds for 30 minutes

# Unsteady diffusion concentration gradient
def unsteady_diffusion_rectangular(y, t, D, L, C0):
    theta = 1 - 2 * sum((-1)**n * np.cos((n + 1/2) * np.pi * y / L) *
                         np.exp(-(n + 1/2)**2 * np.pi**2 * D * t / L**2) for n in range(1000))
    C = C0 * theta
    return C

def concentration_gradient(y, t, D, L, C0):
    # Using the analytical solution for the concentration profile
    # and differentiating with respect to position 'y'
    dC_dy = -2 * sum(
        ((-1)**n * (n + 1/2) * np.pi / L * np.sin((n + 1/2) * np.pi * y / L) *
        np.exp(-(n + 1/2)**2 * np.pi**2 * D * t / L**2)) for n in range(100)
    )
    return dC_dy

# Function to calculate flux using Fick's first law: J = -D * dC/dy
def calculate_zyn_flux(D, gradient):
    return -D * gradient

# Mass transfer across membrane:
def zyn_diffusion(t, conc):
    grad_interface = concentration_gradient(0, t, D_zyn, L, C0z)
    print(grad_interface)
    flux = calculate_zyn_flux(D_zyn, grad_interface)
    print(flux)
    return .014*.028*flux 

### Cigarettes and Vapes
# Constants
u = .045 # units P; blood viscosity ranges from 3.5 to 5.5 cp
R_c = 0.0004 # units cm; radius of the capillary
D_n = .00023 # units cm; diameter of nicotine
R_n = D_n/2 # units cm; radius of nicotine
f = 6*np.pi*u*R_n
K = 1.38065*10**(-16) # Boltzman Constant; g*cm^2*K^-1*s^-2
T = 310 # units K; body temp assumed 37 degrees celsius
Dij = (K*T)/f
v = 1.46/10 #cm/s; average velocity in the capillaries
v_max = 2*v #units cm/s; assuming pousielle flow
z_cap = .06 #cm; length of the capillary
L = R_c*2 #cm; diameter of the capillary
Pe = L*v / Dij

# Model cigarette "puffing"
NUM_ZEROS = 46000
max_molarity = 1.12*10**-6 #units mol/Liter
window = (max_molarity)*signal.windows.general_gaussian(2000, p=1.5, sig=800) #time in ms
offtime = np.zeros(NUM_ZEROS)
co_cigarette = np.concatenate([window, offtime, window, offtime, window, offtime, window, offtime, window, offtime, window, offtime, window, offtime, window, offtime, window, offtime, window, offtime, window])
time_cigarette = np.linspace(0, len(co_cigarette), len(co_cigarette))

window2 = (max_molarity)*signal.windows.general_gaussian(2000, p=1.7, sig=500) #time in ms
offtime2 = np.zeros(2000)
co_onepuff = np.concatenate([offtime2, window2, offtime2])
time_onepuff = np.linspace(0, len(co_onepuff), len(co_onepuff)+1)

# Define Flux equation
def calculate_flux(Dij, Co, R_c, Pe, z_cap):
    return (((0.6783 * Dij * Co) / R_c) * (Pe * (R_c / z_cap))**(1/3))

# Define uptake equation
def calculate_uptake(Co1, t):
    return (((((0.6783 * Dij * Co1) / R_c) * (Pe * (R_c / z_cap))**(1/3))*t)*A)