#!/usr/bin/env python
# coding: utf-8

"""
Author: Savannah Wright
Date: 2025-10-13
Description: Calculates the shear modulus of a neutron star crust using a distribution of nuclei and a single nucleus at different densities. Equation used for OCP calculations taken from The Astrophysical Journal, 375:679-686,1991 July 10.
"""

import numpy as np
import pandas as pd


#Functions for calculating shear modulus for MCP and OCP
def MCP_Distribution_calculations(Zmin, Zmax, Amin, Amax, DeltZ, DeltA, nb):
    MCP_const = (0.1106 * (((4*np.pi)/3)**(1/3)) * (nb**(4/3)))/4
    MCP_Z_vars = ((Zmax)**3) - ((Zmin)**3)
    MCP_A_vars = np.cbrt(-1/Amax)- np.cbrt(-1/Amin)
    MCP_distribution_calc = (MCP_const * MCP_Z_vars * MCP_A_vars) / (DeltZ * DeltA)
    return MCP_distribution_calc

#All MCP calculations will be divided by OCP calculations for ratio
def OCP_calculations(Zavg, Aavg, nb):
    OCP_const = 0.1106 * (((4*np.pi)/3)**(1/3)) * (nb**(4/3))
    OCP_vars = OCP_const * Aavg**(-4/3) * Zavg**2
    OCP_calc_avg = OCP_const * OCP_vars
    return OCP_calc_avg


#Read-In for .csv data for shear modulus at varied density calculations
df = pd.read_csv(".5MeV_Hoa_Results.csv", sep=',', 
names=["nb","Tmelt","Zavg","Aavg","DeltZ","DeltA","Zmax","Zmin","Amax","Amin","rho6"])
np.set_printoptions(suppress=True)
nb_arr = df["nb"].to_numpy()
Tmelt_arr = df["Tmelt"].to_numpy()
Zavg_arr = df["Zavg"].to_numpy()
Aavg_arr = df["Aavg"].to_numpy()
DeltZ_arr = df["DeltZ"].to_numpy()
DeltA_arr = df["DeltA"].to_numpy()
Zmax_arr = df["Zmax"].to_numpy()
Zmin_arr = df["Zmin"].to_numpy()
Amax_arr = df["Amax"].to_numpy()
Amin_arr = df["Amin"].to_numpy()


#Fields for computations (except Tmelt)
nb = nb_arr
Tmelt = Tmelt_arr
Zavg = Zavg_arr
Aavg = Aavg_arr
DeltZ = DeltZ_arr
DeltA = DeltA_arr
Zmax = Zmax_arr
Zmin = Zmin_arr
Amax = Amax_arr
Amin = Amin_arr

#For appending calculations
MCP_Avg_results = []
OCP_Avg_results = []


for i in range(len(df)):
    Zmin = Zmin_arr[i]
    Zmax = Zmax_arr[i]
    Amin = Amin_arr[i]
    Amax = Amax_arr[i]
    DeltA = DeltA_arr[i]
    DeltZ = DeltZ_arr[i]
    Zavg = Zavg_arr[i]
    Aavg = Aavg_arr[i]
    nb = nb_arr[i]
    Tmelt = Tmelt_arr[i]
    
    #Calculations
    MCAvg = MCP_Distribution_calculations(Zmin, Zmax, Amin, Amax, DeltZ, DeltA, nb)
    OCAvg = OCP_calculations(Zavg, Aavg, nb)
    
    #Append to lists
    MCP_Avg_results.append(MCAvg)
    OCP_Avg_results.append(OCAvg)

#Add results to DataFrame
df['MCP_Average'] = MCP_Avg_results
df['OCP_Average'] = OCP_Avg_results


#Add ratio columns
df['Ratio: MCPAvg_OCAvg'] = df['MCP_Average'] / df['OCP_Average']
df.to_csv('Shear Modulus-Varied Density Results.csv', index=False)



