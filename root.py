# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 15:46:14 2020

@author: HP
"""
import numpy as np
import matplotlib.pyplot as plt 
from CoolProp.CoolProp import PropsSI
import inputs as I
import initialize
import find
import update
import generate


Tg,Txy,Tyz,Tzx = initialize.temperature()
# print("Tg")
# print(".................")
# print(Tg[:,:,1]- 273.15)
# print(".................")
# print(Tg[:,:,1]- 273.15)
# print(".................")
# print(Tg[:,:,2]- 273.15)

# print(".................")
# print("Tyz")
# print(".................")
# print(Tyz[:,:,1]- 273.15)
# print(".................")
# print(Tyz[:,:,1]- 273.15)
# print(".................")
# print(Tyz[:,:,2]- 273.15)

Pg, P,Pgstat,Pstat = initialize.pressure()
# print("Pg")
# print(".................")
# print(Pg[:,:,1]/1e5)

# print(".................")
# print("P")
# print(".................")
# print(P[:,:,1]/1e5)

Pg, P = update.pressure(P,Pg,Pgstat,Pstat,Tg)

# print("Pg")
# print(".................")
# print(Pg[:,:,1]/1e5)

# print(".................")
# print("P")
# print(".................")
# print(P[:,:,1]/1e5)

# print(".................")
# print("Rhot")
R = generate.resistance(P,Tyz,Pg,Tg)

ordinate = I.delta_x*np.arange(0,I.a+1)
ordinateg = I.delta_x*np.arange(0,I.a) + I.delta_x*0.5

thmean = np.zeros(np.size(ordinate))
tcmean = np.zeros(np.size(ordinate))
phmean = np.zeros(np.size(ordinate))
pcmean = np.zeros(np.size(ordinate))

res = 10
count = 0
limiter = 0
while(res>1e-3 and (limiter<6)):
    print("iteration = ",count)
    A,C = update.linear_equation_system(R,Tg,Tyz,Pg,P)
    print("inverting matrix")
    sol = np.linalg.inv(A)
    print("multiplyig matrix")
    target = np.matmul(sol,C)
    target4 = target
    A4 = A
    Tg,Tyz = update.temperature(target,Tg,Tyz)

  
    
    plt.plot(ordinate,Tyz[:,1,1] -273)
    plt.plot(ordinate,Tyz[:,3,1] -273)

    sumdhhot = 0
    sumdhcold =0
    for height in range(2,I.b -1,4):
        # print("height = ",height)
        for width in range(1,I.c,2):
            # print("width = ",width)
            sumdhhot = sumdhhot + abs(PropsSI("H","T",Tyz[0,height -1,width],"P",P[0,height-1,width],"co2") - PropsSI("H","T",Tyz[I.a,height -1,width],"P",P[I.a,height -1,width],"co2"))
            sumdhcold =sumdhcold+ abs(PropsSI("H","T",Tyz[0,height +1,width],"P",P[0,height +1,width],"co2") - PropsSI("H","T",Tyz[I.a,height +1,width],"P",P[I.a,height +1,width],"co2"))    


    res = abs((sumdhhot/sumdhcold) - 1)
    print("res = ",res)
    count = count + 1
    print("updating R")
    R = generate.resistance(P,Tyz,Pg,Tg)
    Pg, P = update.pressure(P,Pg,Tg,Pstat,Tg)

    limiter = limiter + 1
plt.show()


for height in range(2,I.b -1,4):
    # print("height = ",height)
    for width in range(1,I.c,2):
        # print("width = ",width)
        plt.plot(ordinate,Tyz[:,height+1,width] -273,'b')
        plt.plot(ordinate,Tyz[:,height-1,width] -273,'orange')
        thmean = thmean + Tyz[:,height-1,width]
        tcmean = tcmean + Tyz[:,height+1,width]
        phmean = thmean + P[:,height-1,width]
        pcmean = tcmean + P[:,height+1,width]

plt.show()

thmean = thmean/(I.b_ch*I.c_ch)
tcmean = tcmean/(I.b_ch*I.c_ch)

plt.plot(ordinate,thmean-273,'orange')
plt.plot(ordinate,tcmean-273,'b')
plt.show()

plt.plot(ordinate,phmean/1e5,'b')
plt.show()
print("Thout = ",thmean[I.a]- 273.15)
print("deltap hot = ",(phmean[0]- phmean[I.a])/1e5)