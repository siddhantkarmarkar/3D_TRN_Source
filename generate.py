# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 23:23:00 2020

@author: HP
"""
import numpy as np
import matplotlib.pyplot as plt 
from CoolProp.CoolProp import PropsSI
import inputs as I
import initialize
import find
import BC

def resistance(P,Tyz,Pg,Tg): #generates resistance matrix
    a = I.a
    b = I.b
    c = I.c
    R = np.zeros([a,b,c,3,2])


    for i in range(0,a):

        for j in range(0,b):
            for k in range(0,c):

                if((np.mod(j,4) == 1 and np.mod(k,2) == 1) or (np.mod(j,4) == 3 and np.mod(k,2) ==1 )): # Fluid channel

                    htc = find.k_liquid(Tg[i,j,k],Pg[i,j,k],i,j,k)
                    htc = htc*find.Nu(Tg[i,j,k],Pg[i,j,k],i,j,k)
                    htc = htc/I.d # Convective Heat transfer coefficient - phi                                            
                    R[i,j,k,1,0] = a/(htc*I.d*I.L)#Rfy0
                    R[i,j,k,1,1] = a/(htc*I.d*I.L)#Rfy1
                    R[i,j,k,2,0] = a/(htc*I.d*I.L)#Rfz0
                    R[i,j,k,2,1] = a/(htc*I.d*I.L)#Rfz1                                             
                    
            
                if(np.mod(j,2) == 0 and np.mod(k,2) == 1): # Type 1 Solid node
                    deltax = I.L/a
                    deltay = I.t
                    deltaz = I.d
                    R[i,j,k,0,0] = deltax/(2*I.k_solid*deltay*deltaz)# Rwx1
                    R[i,j,k,0,1] = deltax/(2*I.k_solid*deltay*deltaz)# Rwx2
                    R[i,j,k,1,0] = deltay/(2*I.k_solid*deltax*deltaz)# Rwy1
                    R[i,j,k,1,1] = deltay/(2*I.k_solid*deltax*deltaz)# Rwy2
                    R[i,j,k,2,0] = deltaz/(2*I.k_solid*deltay*deltax)# Rwz1
                    R[i,j,k,2,1] = deltaz/(2*I.k_solid*deltay*deltax)# Rwz2
                        
                if(np.mod(j,2) == 0 and np.mod(k,2) == 0): # Type 2 Solid node
                    deltax = I.L/a
                    deltay = I.t
                    deltaz = I.w
                    R[i,j,k,0,0] = deltax/(2*I.k_solid*deltay*deltaz)# Rwx1
                    R[i,j,k,0,1] = deltax/(2*I.k_solid*deltay*deltaz)# Rwx2
                    R[i,j,k,1,0] = deltay/(2*I.k_solid*deltax*deltaz)# Rwy1
                    R[i,j,k,1,1] = deltay/(2*I.k_solid*deltax*deltaz)# Rwy2
                    R[i,j,k,2,0] = deltaz/(2*I.k_solid*deltay*deltax)# Rwz1
                    R[i,j,k,2,1] = deltaz/(2*I.k_solid*deltay*deltax)# Rwz2
                    
                if(np.mod(j,2) == 1 and np.mod(k,2) == 0): # Type 3 Solid node
                    deltax = I.L/a
                    deltay = I.d
                    deltaz = I.w
                    R[i,j,k,0,0] = deltax/(2*I.k_solid*deltay*deltaz)# Rwx1
                    R[i,j,k,0,1] = deltax/(2*I.k_solid*deltay*deltaz)# Rwx2
                    R[i,j,k,1,0] = deltay/(2*I.k_solid*deltax*deltaz)# Rwy1
                    R[i,j,k,1,1] = deltay/(2*I.k_solid*deltax*deltaz)# Rwy2
                    R[i,j,k,2,0] = deltaz/(2*I.k_solid*deltay*deltax)# Rwz1
                    R[i,j,k,2,1] = deltaz/(2*I.k_solid*deltay*deltax)# Rwz2  
        
    return(R)

