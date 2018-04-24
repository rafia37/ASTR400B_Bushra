import numpy as np
import matplotlib.pyplot as plt
import sys, pdb
sys.path.insert(0, '../HW2/') #Making python search the HW2 directory for ReadFile program
from ReadFile import Read


MW_data = Read('../../MW_000.txt')[2]
M31_data = Read('../../M31_000.txt')[2]

#use CenterOfMass.py to calculate M31's galactic center at Snapshot 0 (xcen0, ycen0, zcen0)
#Collect indices of disk particles that are within 7 & 9 kpc of (xcen0, ycen0, zcen0)
#Last two steps will give the initial 
#Loop over all desired the snaphots
#-> convert simulation coordinates to the com coordinates using CenterOfMass.py
#-> calculate position and kinematics of desired particles (this answers questions 2 & 3)
#Learn to calculate escape velocity of a galaxy
#if some particles have velocities higher than the escape velocity of M31, they are unbound
#Calculate the fraction of unbound particles (answers question 4)

def magnitude(x, y, z = 0):
    return np.sqrt(x**2 + y**2 + z**2)

#creating M31 center of mass object using disk particles at snapshot 0
com = CenterOfMass("../../M31_000.txt", 2)

#Present day com position & velocity of M31
xcom, ycom, zcom = com.COM_P(1)
vxcom, vycom, vzcom = com.COM_V(xcom, ycom, zcom, 15)

#Getting position & velocity of all particle in COM reference frame
x  = com.x - xcom
y  = com.y - ycom
z  = com.z - zcom
vx = com.vx - vxcom
vy = com.vy - vycom
vz = com.vz - vzcom

r = magnitude(x - xcom, y - ycom, z - zcom) #An array that contains distance of all particles from center of mass

eqt_vel = magnitude(vx, vy) #An array that contains equatorial velocity of all particles from center of mass


#R = 8.29 kpc: rmin = R-0.1R = 7.46kpc, rmax = R+0.1R = 9.12kpc
#V = 239 km/s; vmin = V-0.1V = 215.1 km/s, vmax = V+0.1V = 262.9
mask = np.where((r > 7.46) & (r < 9.12) & (eqt_vel > 215.1) & (eqt_vel < 262.9) & (np.abs(vz) < 30))

