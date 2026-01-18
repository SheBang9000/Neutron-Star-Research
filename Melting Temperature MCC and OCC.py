#!/usr/bin/env python
# coding: utf-8

"""
Author: Savannah Wright
Date: 2024-10-15
Description: Calculates the melting temperature of a neutron star crust at a specified density using both a distribution of nuclei and single nucleus. Equation used for OCP calculations taken from https://doi.org/10.12942/lrr-2008-10. 
"""

import numpy as np
import pandas as pd

#Functions
def Tmelt_avg(Amin,Amax,Zmin,Zmax,rho6,DeltZ,DeltA):
    TA = ((3/2)*(Amax**(2/3))) - ((3/2)*(Amin**(2/3)))
    TZ = ((1/3)*(Zmax**3)) - ((1/3)*(Zmin**3))
    Tmelt_calc = ((1.3)*(10**5)) * (rho6**(1/3)) * TA * TZ
    Tmelt_calc = Tmelt_calc/(DeltA*DeltZ)
    return Tmelt_calc

def Tmelt_sing_Median(Zavg,Aavg,rho6):
    Tmelt_calc = (1.3*10**5)*(Zavg**2)*((rho6/Aavg)**(1/3))
    return Tmelt_calc

def Tmelt_sing_Peak(PeakZ,PeakA,rho6):
    Tmelt_calc = (1.3*10**5)*(PeakZ**2)*((rho6/PeakA)**(1/3))
    return Tmelt_calc

def Tmelt_sing_Orig(OCPZ,OCPA,rho6):
    Tmelt_calc = (1.3*10**5)*(OCPZ**2)*((rho6/OCPA)**(1/3))
    return Tmelt_calc


#File read-in and conversions
Tmelt_Init = pd.read_csv("FWHM_Dist_P_init.csv", sep = ',',
names = ["Temp","nb","rho6","gamma","Amin","Amax","Zmin","Zmax","DeltZ","DeltA","Zavg","Aavg","PeakZ","PeakA"])
np.set_printoptions(suppress=True)
nb_arr = Tmelt_Init["nb"].to_numpy()
rho6_arr = Tmelt_Init["rho6"].to_numpy()
gamma_arr = Tmelt_Init["gamma"].to_numpy()
Amin_arr = Tmelt_Init["Amin"].to_numpy()
Amax_arr = Tmelt_Init["Amax"].to_numpy()
Zmin_arr = Tmelt_Init["Zmin"].to_numpy()
Zmax_arr = Tmelt_Init["Zmax"].to_numpy()
DeltA_arr = Tmelt_Init["DeltZ"].to_numpy()
DeltZ_arr = Tmelt_Init["DeltA"].to_numpy()
Zavg_arr = Tmelt_Init['Zavg'].to_numpy()
Aavg_arr = Tmelt_Init['Aavg'].to_numpy()
PeakZ_arr = Tmelt_Init["PeakZ"].to_numpy()
PeakA_arr = Tmelt_Init["PeakA"].to_numpy()


#Fields to be appended
Tmelt_avg_distribution_results = []
Tmelt_sing_Median_results = []
Tmelt_sing_Peak_results = []


for i in range(len(Tmelt_Init)):
    nb = nb_arr[i]
    rho6 = rho6_arr[i]
    gamma = gamma_arr[i]
    Amin = Amin_arr[i]
    Amax = Amax_arr[i]
    Zmin = Zmin_arr[i]
    Zmax = Zmax_arr[i]
    DeltA = DeltA_arr[i]
    DeltZ = DeltZ_arr[i]
    Zavg = Zavg_arr[i]
    Aavg = Aavg_arr[i]
    PeakZ = PeakZ_arr[i]
    PeakA = PeakA_arr[i]
    
    #Calculations
    Tmelt_avg = Tmelt(Amin,Amax,Zmin,Zmax,rho6,DeltA,DeltZ)
    TsM = Tmelt_sing_Median(Zavg,Aavg,rho6)
    TsP = Tmelt_sing_Peak (PeakZ,PeakA,rho6)
    
    #Append to lists
    Tmelt_distribution_results.append(Tmelt_avg)
    Tmelt_sing_Median_results.append(TsM)
    Tmelt_sing_Peak_results.append(TsP)
    
#Add results to DataFrame
Tmelt_Init['Tmelt'] = Tmelt_distribution_results
Tmelt_Init['Tmelt_sing_Median'] = Tmelt_sing_Median_results
Tmelt_Init['Tmelt_sing_Peak'] = Tmelt_sing_Peak_results

#Add ratio to Columns
Tmelt_Init['Median Ratio'] = Tmelt_Init['Tmelt'] / Tmelt_Init['Tmelt_sing_Median']
Tmelt_Init['Peak Ratio'] = Tmelt_Init['Tmelt'] / Tmelt_Init['Tmelt_sing_Peak']

#File save
Tmelt_Init.to_csv('FWHM_Distribution_P_Calc.csv', index=False





