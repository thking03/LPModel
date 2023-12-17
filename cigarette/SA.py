import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy import signal
# Surface Area Modeling

NUM_ZEROS = 46000
max_molarity = 1.12*10**-6 #units mol/Liter
max_SA = 1000000 # cm^2; 100% SA
min_SA = 510000 # cm^2; exhalation
SA_increment = (max_SA-min_SA)*1/3
delta_SA = max_SA - min_SA

window = (max_SA)*signal.windows.general_gaussian(2000, p=1.7, sig=400) #time in ms
window1 = (SA_increment+min_SA)*signal.windows.general_gaussian(2000, p=1.7, sig=400) #time in ms
window2 = ((2*SA_increment)+min_SA)*signal.windows.general_gaussian(2000, p=1.7, sig=400) #time in ms
offtime = np.zeros(NUM_ZEROS)
SA_cigarette = np.concatenate([window, offtime, window1, offtime, window2])
time_cigarette = np.linspace(0, len(SA_cigarette), len(SA_cigarette))

# Plot SA vs. Time 
plt.plot(time_cigarette, SA_cigarette)
plt.xlabel('Time (ms)')
plt.ylabel('Surface Area (cm^2)')
plt.axhline(y=510000, color='r', linestyle='--', linewidth=2)
plt.axhline(y=1000000, color='r', linestyle='--', linewidth=2)
plt.legend()
plt.title('Absolute Surface Area of Nicotine in Alveoli vs. Time over 3 Puffs')
plt.grid(True)
plt.show()

# change in SA graph
window3 = (delta_SA)*signal.windows.general_gaussian(2000, p=1.7, sig=400) #time in ms
window4 = (SA_increment)*signal.windows.general_gaussian(2000, p=1.7, sig=400) #time in ms
window5 = (2*SA_increment)*signal.windows.general_gaussian(2000, p=1.7, sig=400) #time in ms
offtime = np.zeros(NUM_ZEROS)
SA_cigarette2 = np.concatenate([window3, offtime, window4, offtime, window5])
time_cigarette2 = np.linspace(0, len(SA_cigarette), len(SA_cigarette))

# Plot SA vs. Time 
plt.plot(time_cigarette2, SA_cigarette2)
plt.xlabel('Time (ms)')
plt.ylabel('Change in Surface Area (cm^2); 0 = 51 cm^2')
plt.legend()
plt.title('Change in Surface Area of Nicotine in Alveoli vs. Time over 3 Puffs')
plt.grid(True)
plt.show()