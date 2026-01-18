#!/usr/bin/env python
# coding: utf-8

"""
Author: Savannah Wright
Date: 2023-10-19
Description: Compressible Liquid Drop Model (CLDM) program for comparing results to IAEA data. Equation used was taken from The Astrophysical Journal, 170:299-317, 1971 December 1 and tailored for research.
"""

import glob
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import time
import sys
import matplotlib.pyplot as plt

#File import directory
directory = '/home/guest/anaconda3/bin/Work/LSL/*.csv' 

#Universal Constants
hbarc2_over_2m = 20.74
e_sq = 1.4399644
p = 3

#True Set Data for popt 'nuclei' (IAEA)
elemental_data = pd.read_csv("new_data1.csv", sep=',', names=["Neutrons", "Protons", "Element", "Binding Energy", "Error"])
np.set_printoptions(suppress=True)
neutron_arr = elemental_data["Neutrons"].to_numpy()
proton_arr = elemental_data["Protons"].to_numpy()
binding_energy_arr = elemental_data["Binding Energy"].to_numpy()
error_arr = elemental_data["Error"].to_numpy()

A = np.add(neutron_arr,proton_arr)
Z = proton_arr
nuclei = (A,Z)

#Functions
def wsenergy(nuclei_w, n):    
    (Aw,Zw) = nuclei_w
    I = 1-((2*Zw)/Aw)
    delta = zeta * I
    x = (1-delta)/2
    W = a1*n - (a2*n*(x*x+((1.000-x)*(1.000-x)))) #adapted FORTRAN code
    W = W + a3*(n**(alpha3+1.000))
    W = W + a43*(n**(alpha4+1.000))
    W = W - (a4*(n**(alpha3+1.000))*(x*x+((1.000-x)*(1.000-x))))
    W = W - (a44*(n**(alpha4+1.000))*(x*x+((1.000-x)*(1.000-x))))
    temp = a5*(n**(5.000/3.000))*((x**(5.000/3.000))+((1.000-x)**(5.000/3.000)))
    temp = temp + (a6*(n**(5.000/3.000))*((x**(8.000/3.000))+((1.000-x)**(8.000/3.000))))
    W = W + ((3.000/5.000)*((3.000*np.pi*np.pi)**(2.000/3.000))*temp)
    Wkin = 5.742468 * hbarc2_over_2m * (n**(2/3)) * (x**(5/3) + (1.000 - x)**(5/3))
    W = W + Wkin
    return W

def binding_energy(Xb,n0,Ccoul,b,sigma,n1,zeta):
    (Ab,Zb) = Xb
    I = 1-((2*Zb)/Ab)
    delta = zeta * I
    x = (1-delta)/2
    n = n0 + (n1 * (I**2))
    r_n = (((3*Ab)/((4*np.pi)*n))**(1/3))
    BE1 = n * wsenergy(Xb, n)
    BE2 = ((4 / 5) * np.pi * e_sq) * ((n * x) ** 2) * (r_n ** 2) * Ccoul
    BE3 = ((3 * sigma) / r_n) * (((2 ** (p+1) + b)) / (((1 / (x ** p))) + b + (1 / ((1 - x) ** p))))
    BE = -(BE1 + BE2 + BE3) / n 
    BE = BE * 1000
    return BE

#Assigning column names to imported *.csv
df_names = ['Filename','J','L','Ksym','n0','Ccoul','b','sigma','n1','zeta','sigma-delta']
sk_df = pd.DataFrame(columns=df_names)
sk_df.to_csv('LSR.csv', index=False)

#Start time to check program efficiency
start = time.time()
#For Loop
for fname in glob.glob(directory):
#Read-In
    df = pd.read_csv(fname, sep=",")
    Jcol = df.insert(0, "J", str(fname[48:53]))
    Lcol = df.insert(1, "L", str(fname[55:60]))
    Ksymcol = df.insert(2, "Ksym", str(fname[62:67]))
    col1 = df.iloc[:,[3]].values
    col2 = df.iloc[:,[4]].values
    col3 = df.iloc[:,[5]].values
    col4 = df.iloc[:,[6]].values
    col5 = df.iloc[:,[7]].values
    SDcol = df.insert(8, "sigma-delta",float())
    #binding_energy((16,8),n0,Ccoul,b,sigma,n1,zeta,sigma-delta)
#Optimal Parameters
    for row in df:
        try:
            #Variables
            Jcol = float(col1[0])
            Lcol = float(col2[0])
            Ksymcol = float(col3[0])
            t0 = col1[1]
            x0 = col1[2]
            b4 = col1[3]
            t1 = col2[1]
            x1 = col2[2]
            b4p = col2[3]
            t2 = col3[1]
            x2 = col3[2]
            alpha3 = col3[3]
            t3 = col4[1]
            x3 = col4[2]
            alpha4 = col4[3]
            t4 = col5[1]
            x4 = col5[2]
            mass = col5[3]
            zeta = 0.75
        #FORTRAN CALCULATIONS    
            a1 = (1.000/4.000)*t0*(2.000+x0)
            a2 = (1.000/4.000)*t0*(2.000*x0+1.000)
            a3 = (1.000/24.000)*t3*(2.000+x3)
            a4 = (1.000/24.000)*t3*(2.000*x3+1.000)
            a43 = (1.000/24.000)*t4*(2.000+x4)
            a44 = (1.000/24.000)*t4*(2.000*x4+1.000)
            a5 = (1.000/8.000)*(t1*(2.000+x1)+(t2*(2.000+x2)))
            a6 = (1.000/8.000)*(t2*(2.000*x2+1.000)-(t1*(2.000*x1+1.000))) 
            popt, pcov = curve_fit(binding_energy, nuclei, binding_energy_arr)
            sk1 = popt[0]
            sk2 = popt[1]
            sk3 = popt[2]
            sk4 = popt[3]
            sk5 = popt[4]
            sk6 = popt[5]
            sk7 = 96 / (16 + sk3)
        except:
            result_file = open("LSR_BadFiles.txt", "a")
            result_file.write(fname[35:] + "\n")
            result_file.close()
            continue 
             
    sk_df = pd.DataFrame([[fname[35:], Jcol, Lcol, Ksymcol, sk1, sk2, sk3, sk4, sk5, sk6, sk7]])
    sk_df.to_csv('LSR.csv', mode='a', header=None, index=False)

stop = time.time()
result = stop - start
print(result)
print(popt)
perr = np.sqrt(np.diag(pcov))
print(perr)

#Scipy curvefit
popt, pcov = curve_fit(binding_energy, nuclei, binding_energy_arr)
print(popt)


#Plot for comparing original data from IAEA to CLDM
plt.xlabel('A')
plt.ylabel('Binding Energy in keV')
plt.grid()
plt.scatter(A, binding_energy_arr, s=1, color='black', alpha=1)
plt.scatter(A, binding_energy((A,Z), popt[0], popt[1], popt[2], popt[3], popt[4], popt[5]), s=1, alpha=0.5, color='red')
plt.legend(['Data - IAEA','CLDM'])
plt.savefig("Data v. CLDM.png")
plt.show()

