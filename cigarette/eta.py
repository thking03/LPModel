#imports and values
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad

vmax = 0.292 #cm/s
Dij = 4.388*(10**(-10)) #cm^2/s
R = 4*(10**-4) #cm
y = np.linspace(0,1.15*(10**-4), 100) #y in cm, radius of capillary?
z = np.linspace(0, 600*(10**-4), 100) #z in cm, length of capillary
infinity = float('inf')
a=0
b=infinity

#Creating Eta Array (for x-axis)
AIDA = y/(((9*R*Dij*z)/(2*vmax))**(1/3)) #yes, spelled wrong
EtaArray = AIDA
#print(EtaArray) #checks out, values increase from 0 to about 2. This is the correct scale of x-axis

#setting up integral function (relative concentration)
def integrand(z):
    return (np.exp(-z**3))*1.1198

def integrate_function(a, b):
    result = quad(integrand, a, b)
    return result

result = integrate_function(a, b)
CoList = []

for i in AIDA:  #for loop to create array (list) of y-axis valules
    value = integrate_function(i, b)
    CoList.append(value)


integrate_function(0, 100)


#defining integral bounds
#a = EtaArray #integral goes from Eta to infinity
#b = infinity

#plotting function: relative concentration vs. EtaArray
plt.plot(EtaArray, CoList)
plt.title('Concentration vs. Eta')
plt.xlabel('Eta')
plt.ylabel('Relative Concentration (Ci/Co)')
# Show the plot
plt.show()

