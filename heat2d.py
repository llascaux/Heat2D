#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import time
start_time = time.time()

plt.style.use('dark_background')

#size of room
room_x = 10 # (m)
room_y = 10 # (m)

#Constants
sigma = 0.1 #standard deviation of gaussian distribution
rho = 1.225 #denisty of air (kg/m^3)
cp = 1000 #specific heat capacity of air (J/kg*K)
k = 10**(-2) #Thermal conductivity of air (W/m*K)
h = 1 #heat transfer coefficient (W/m^2*K)
delta_x = 0.001 #small change in x (m)
delta_y = 0.001 #small chage in y (m)
lamda = delta_x*h/k 
T_out = 200 #outside temp (K)
w_h = 500 #power of space heater (W)
L = 2 #dimension of heater (m)
m = L/delta_x
total_time = 40 #repetitions of simulation
source_power = w_h/(m**2 * delta_x**3)
source_x = room_x/2 #x-coordinate of source
source_y = room_y/2 #y-coordinate of source
b_temp_1 = 280 #temp of boundary 1 (K)
b_temp_2 = 280 #temp of boundary 2 (K)
b_temp_3 = 280 #temp of boundary 3 (K)
room_temp = 293 #temp of room (K)
delta_t = 1/2 * (delta_x*delta_x)/(2*k/(rho*cp)) #time interval (s)

#Temperature values
T = np.full((room_x,room_y), room_temp)


fig, ax = plt.subplots(figsize=(room_x,room_y))

#Uncomment for equilibrium temp at various points
# t1 = []
# t2 = []
# t3 = []
# ts = np.empty(100, total_time)

#Uncomment for animation
# ims = []

for t in range(total_time):
    # t1.append(T[int(room_x/2),int(room_y/2)])
    # t2.append(T[int(9),int(9)])
    # t3.append(T[int(3),int(4)])
    for  j in range(1,room_y-1):
        for i in range(1,room_x-1):
            #Equation
            Eq_1 = (T[i-1,j]-2*T[i,j]+T[i+1,j])/(delta_x*delta_x)
            Eq_2 = (T[i,j-1]-2*T[i,j]+T[i,j+1])/(delta_y*delta_y)
            T[i,j] = T[i,j] + delta_t*(k/(rho*cp))*(Eq_1+Eq_2) + delta_t*source_power*(np.exp(-(delta_x*delta_x*(i-source_x)**2)/sigma - (delta_y*delta_y*(j-source_y)**2)/sigma))/(rho*cp)  
            
    for i in range(0,room_x):
        T[i, room_y-1] = (T[i, room_y-1] + lamda*T_out)/(1+lamda)
        T[i, 0] = b_temp_1
    
    for j in range(0,room_y):
        T[room_x-1, j] = b_temp_2
        T[0, j] = b_temp_3
        
         
    plt.clf()
    plt.pcolormesh(T, cmap=plt.cm.jet, vmin=250, vmax=320, shading='gouraud')    
    plt.colorbar()
    plt.pause(0.001)
    
    #Uncomment for animation
    # im =  ax.pcolormesh(T, cmap=plt.cm.jet, vmin=250, vmax=320, shading='gouraud')
    # ims.append([im])

#Uncomment for animation
# ani = animation.ArtistAnimation(fig, ims)
# ani.save('out.gif', writer='imagemagick')
# plt.close()
# ani

#Uncomment for equilibrium temp at various points
# plt.title("Equilibrium Temperature of various position")
# plt.xlabel("Number of times the loop ran")
# plt.ylabel("Temperature (K)")
# plt.xticks(range(0, 100, 2))
# plt.yticks(range(200, 400, 1))
# plt.grid()
# plt.plot(t1)
# plt.plot(t2)
# plt.plot(t3)
# plt.legend(["[5,5]-position","[9,9]-position","[3,4]-position"])

#Uncomment for default plot
print(time.time()-start_time)
plt.show()

