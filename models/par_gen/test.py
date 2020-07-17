import math

# Fixed parameters
N = 200     # Number of participants
Tmax = 120  # total duration of the simulaion

# Tunable parameters
Tinf = 30   # infectious time
I0 = 1      # number of initial cases
Itot = 150  # total number of cases
cr = 0.005    # per capita contact rate c, 
            # cr x N is the number of contacts per unit of time an infectious individual make

S0 = N - I0
Send = N - Itot

z = Itot/math.log(S0/Send)
print("gamma/beta = ", z)
print("R0 = ", S0/z)

Imax = N - z + z * math.log(z) - z * math.log(S0)
print("IMax = ", Imax)

beta = 1/(z*Tinf)
print("beta = ", beta)

p = beta/cr  # probability that a contact with a susceptible individual results in transmission
print("Probability of infection = ", p)