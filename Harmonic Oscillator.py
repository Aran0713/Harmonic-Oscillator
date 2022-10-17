#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapz

# N = number of points
N = 2000
# h = step size
h = 8.0 / N
# L = x goes from -L to L
L = 4
# Constants
q = 1.0 / (h**2)
u = -0.5*q

# Arrays
V = np.empty([N+1])
M = np.empty([N+1,N+1])
O = np.empty([N+1])
n = np.empty([50])
calceig = np.empty([50])
theeig = np.empty([50])
erroreig = np.empty([50])

# txt files
#f1 = open('doubleSiO2eigenvalues.txt','w')

# Filling up values for x
x = np.linspace(-L,L,N+1)

# Filling up values for array V
for i in range (0,N+1):
    V[i] = 0.5*x[i]**2

# Filling up values for matrix M
for i in range (0,N+1):
    for j in range (0,N+1):
        if j == i:
            M[i][j] = V[i] + q
        elif j == i-1 or j == i+1:
            M[i][j] = u
        else:
            M[i][j] = 0
            
# Diagonalizing matrix M
spectrum,v=np.linalg.eigh(M)
    
# Array for calculatee and theoretical eignevalues
for i in range(0,50):
    n[i] = i
    calceig[i] = spectrum[i]
    theeig[int(n[i])] = n[i] + 0.5
    
# Array for error in eigenvalues
for i in range(0,50):
    erroreig[i] = 100*((theeig[i] - calceig[i]) / theeig[i])
for j in range(0,50):
    if erroreig[j] < 0:
        erroreig[j] = -erroreig[j]

# Defining function for calculated wavefunction
def calcpsi(j):
    for i in range(0,N+1):
        O[i] = v[i][j-1]
    return O

# Defining function for normalized calculated wavefunction
def psi(j):
    A = 1 / np.sqrt(trapz(calcpsi(j)**2,x))
    for i in range (0,N+1):
        O = A * calcpsi(j)
    return O

#################### Plots ####################

# Plot of potential
plt.plot(x,V,color='k')
for  i in range (0,7):
    plt.hlines(spectrum[i],-L,L,color="gray",linestyle='--')
    plt.text(0.8*L,1.05*spectrum[i],"n = "+str(i))
    plt.plot(x,(psi(i+1)**2)+spectrum[i],lw='0.8')
plt.xlabel('x')
plt.ylabel('V(x)')
plt.title('Plot of Potential for Harmonic Oscillator')
plt.show()

# Plot of energy
plt.plot(n,calceig,'r',label="Calculated")
plt.plot(n,theeig,color='b',linestyle='--',label="Theoretical, E = n + 0.5")
plt.xlabel('n')
plt.ylabel('Energy, E')
plt.title('Plot of Energy for Harmonic Oscillator, L = '+str(L)+', N = '+str(N))
plt.legend()
plt.show()

# Plot of error in energy
plt.plot(n,erroreig,'b')
plt.xlabel('n')
plt.ylabel('% error')
plt.title('Plot of error of eigenvalues, L = '+str(L)+', N = '+str(N))
plt.show()