#!/usr/bin/env python
# coding: utf-8

"""
Author: Savannah Wright
Date: 2023-07-13
Description: Calculates the binding energy of nuclei based on nuclear data from the International Atomic Energy Association (IAEA)
"""

import sys
import pandas as pd
import numpy as np
import matplotlib as plt
from scipy.optimize import curve_fit


# Functions
def bindenergy(X,avol,asrf,aclb,asym,apair):
    (A,Z) = X
    BE = (avol * A) - (asrf * (A**(2/3))) - (((aclb * Z * (Z-1)))/(A**(1/3))) - (asym * (((A-(2*Z))**2)/A)) - (apair / (A**(3/4)))
    BE = BE * 1000
    BE = BE / A
    return BE
    
def fit(x,a,b,c):
    return a * np.exp(b * x) + c


# .csv Read-In Section
read_file = pd.read_csv(r"/home/guest/Desktop/Research/new_data.txt")
read_file.to_csv (r"/home/guest/Desktop/Research/new_data.csv", index=None)


# Data and Arrays
elemental_data = pd.read_csv("/home/guest/Desktop/Research/new_data.csv", sep="\t", names=["Neutrons", "Protons", "Element", "Binding Energy", "Error"])
np.set_printoptions(suppress=True)
neutron_arr = elemental_data["Neutrons"].to_numpy()
proton_arr = elemental_data["Protons"].to_numpy()
binding_energy_arr = elemental_data["Binding Energy"].to_numpy()
error_arr = elemental_data["Error"].to_numpy()

A = np.add(neutron_arr,proton_arr)
Z = proton_arr
X = (A,Z)


#POPT, PCOV
popt, pcov = curve_fit(bindenergy, X, binding_energy_arr, sigma=error_arr)
print(popt)
perr = np.sqrt(np.diag(pcov))
print(perr)

#Adjustable plot parameters for comparing each additional parameter to the original data
plt.xlabel("A")
plt.ylabel("Binding Energy in keV")
plt.grid()
plt.scatter(A, binding_energy_arr, s=1, color='black', alpha=0.4)
#plt.scatter(A, bindenergy((A,Z),popt[0],0,0,0,0), color='red', s=1, alpha=0.5)
#plt.scatter(A, bindenergy((A,Z),popt[0],popt[1],0,0,0), color='orange', s=1, alpha=0.5)
#plt.scatter(A, bindenergy((A,Z),popt[0],popt[1],popt[2],0,0), color='green', s=1, alpha=0.5)
#plt.scatter(A, bindenergy((A,Z),popt[0],popt[1],popt[2],popt[3],0), color='cyan', s=1, alpha=0.5)
plt.scatter(A, bindenergy((A,Z),popt[0],popt[1],popt[2],popt[3],popt[4]), color='orange', s=1, alpha=0.3)
plt.legend(["Data-IAEA","SEMF"])
plt.savefig("SEMF.png")
plt.show()





