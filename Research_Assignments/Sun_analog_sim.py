import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pdb 
from ReadFile import Read
from CenterOfMass import CenterOfMass


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
com = CenterOfMass("HighRes/M31_000.txt", 2)



def get_com_param(com_object):
    
    #Present day com position & velocity of M31
    xcom, ycom, zcom = com_object.COM_P(0.1, 4.0)
    vxcom, vycom, vzcom = com_object.COM_V(xcom, ycom, zcom, 15)
    
    #Getting position & velocity of all particle in COM reference frame
    x  = com_object.x - xcom
    y  = com_object.y - ycom
    z  = com_object.z - zcom
    vx = com_object.vx - vxcom
    vy = com_object.vy - vycom
    vz = com_object.vz - vzcom
    
    param = np.array([x, y, z, vx, vy, vz])

    r = magnitude(x, y, z).value #An array that contains distance of all particles from center of mass

    eq_vel = magnitude(vx, vy).value #An array that contains equatorial velocity of all particles from center of mass
    
    return r, eq_vel, vz, param


r, eq_vel, vz, param = get_com_param(com)


#R = 8.29 kpc: rmin = R-0.1R = 7.46kpc, rmax = R+0.1R = 9.12kpc
#V = 239 km/s; vmin = V-0.1V = 215.1 km/s, vmax = V+0.1V = 262.9
mask1 = np.where((r > 7.46) & (r < 9.12) & (eq_vel > 215.1) & (eq_vel < 262.9) & (np.abs(vz.value) < 30))
mask2 = np.where((r > 7.46) & (r < 9.12) & (eq_vel > 215.1) & (eq_vel < 262.9))
mask3 = np.where((r > 7.46) & (r < 9.12))


print(len(mask1[0]), len(mask2[0]), len(mask3[0]))

"""
sunx  = param[0][mask1]
suny  = param[1][mask1]
sunz  = param[2][mask1]
sunvx = param[3][mask1]
sunvy = param[4][mask1]
sunvz = param[5][mask1]


fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, sharey = True)

im1 = ax1.scatter(param[0], param[1], c = eq_vel, cmap = 'magma')
im2 = ax2.scatter(sunx, suny, c = eq_vel[mask1], cmap = 'magma')
fig.colorbar(im2)
ax1.set_title('X-Y position of all disk \n particles of M31')
ax2.set_title('X-Y position of sun analogs \n in M31\'s disk')
ax1.set_ylabel('Y-position (kpc)')
ax1.set_xlabel('X-position (kpc)')
ax2.set_xlabel('X-position (kpc)')

fig.savefig('masking1.png')
"""

time = np.array([])
for i in np.arange(0, 800, 5):
    
    ilbl = '000' + str(i)
    ilbl = ilbl[-3:]
    folder = 'HighRes/'
    fname = folder + 'M31_' + ilbl + '.txt'
    
    NewCOM = CenterOfMass(fname, 2)
    t = NewCOM.time/1000.0    #converting Myr to Gyr
    print(t)
    time = np.append(time, t)
    r, eq_vel, vz, param = get_com_param(NewCOM)

pdb.set_trace()