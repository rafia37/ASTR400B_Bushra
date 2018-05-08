import numpy as np
from astropy.io import ascii
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pdb 
from ReadFile import Read
from CenterOfMass import CenterOfMass

# Choose (by setting to True) the type of analysis/plotting you want to do
masking           = False     #Makes a plot of the selected (masked) particles
interesting_times = False     #Makes the seperation plot with selected times
snapshots         = False     #Position mapping at selected times
histograms        = False     #Particle distribution at selected times
time_series       = True      #Position and velocity of 5 particles as a function of time

#A function to calculate magnitude of 2D or 3D vectors
def magnitude(x, y, z = 0):
    return np.sqrt(x**2 + y**2 + z**2)

#creating M31 center of mass object using disk particles at snapshot 0
com = CenterOfMass("HighRes/M31_000.txt", 2)

#A function that calculates and returns necessary parameters from COM object
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

#Storing parameters fr snapshot 0
r0, eq_vel0, vz0, param0 = get_com_param(com)


#Creating masks to select particles according to the following criterion:
#R = 8.29 kpc: rmin = R-0.1R = 7.46kpc, rmax = R+0.1R = 9.12kpc
#V = 239 km/s; vmin = V-0.1V = 215.1 km/s, vmax = V+0.1V = 262.9
mask1 = np.where((r0 > 7.46) & (r0 < 9.12) & (eq_vel0 > 215.1) & (eq_vel0 < 262.9) & (np.abs(vz0.value) < 30))
mask2 = np.where((r0 > 7.46) & (r0 < 9.12) & (eq_vel0 > 215.1) & (eq_vel0 < 262.9))
mask3 = np.where((r0 > 7.46) & (r0 < 9.12))
print(len(r0), len(mask1[0]), len(mask2[0]), len(mask3[0]))

#Snapshots I'm interested in
int_times = [140, 275, 335, 410, 425, 545]



## The 5 different kinds of plotting/analysis that can be done with this code
if masking:
    
    #parameters for sun-like particles
    sunx  = param0[0][mask1]
    suny  = param0[1][mask1]
    sunz  = param0[2][mask1]
    sunvx = param0[3][mask1]
    sunvy = param0[4][mask1]
    sunvz = param0[5][mask1]


    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, sharey = True)

    im1 = ax1.scatter(param0[0], param0[1], c = eq_vel0, cmap = 'magma')
    im2 = ax2.scatter(sunx, suny, c = eq_vel0[mask1], cmap = 'magma')
    fig.colorbar(im2)
    ax1.set_title('Present X-Y position of all \n disk particles of M31')
    ax2.set_title('Present X-Y position of sun \n analogs in M31\'s disk')
    ax1.set_ylabel('Y-position (kpc)')
    ax1.set_xlabel('X-position (kpc)')
    ax2.set_xlabel('X-position (kpc)')

    fig.savefig('masking1.png')
    plt.close()


    
    
if interesting_times:
    
    #Using orbit file from homework 6 to make the separartion plot
    mw = ascii.read('../HW6/Orbit_MW.txt')
    m31 = ascii.read('../HW6/Orbit_M31.txt')
    sep = magnitude(mw['x']-m31['x'], mw['y']-m31['y'], mw['z']-m31['z'])
    
    plt.plot(m31['t'], sep)
    plt.plot(m31['t'][int(int_times/5)], sep[int(int_times/5)], '*') #marking selected times
    plt.xlabel('Time (Gyr)')
    plt.ylabel('Separation (Kpc)')
    plt.title('Separation Between Milky Way And Andromeda \n With Points Of Interest')
    
    plt.savefig('interesting_times.png')
    plt.close()
    

    
if snapshots:
    fig, ax_tmp = plt.subplots(nrows=2, ncols=3, sharex = True, sharey = True)
    
    #looping over seleted snapshots
    for i, t in enumerate(int_times):
        
        #creating COM object, regenerating parameters and masked parameters for this snapshot
        com = CenterOfMass('HighRes/M31_%s.txt' % str(t), 2)
        r, eq_vel, vz, param = get_com_param(com)
        sunx  = param[0][mask1]
        suny  = param[1][mask1]
        sunz  = param[2][mask1]
        sunvx = param[3][mask1]
        sunvy = param[4][mask1]
        sunvz = param[5][mask1]
        
        ax = ax_tmp[0][i] if i<3 else ax_tmp[1][i-3]
        ax.plot(param[0], param[1], 'c.', alpha = 0.5)  #All particles
        ax.plot(sunx, suny, 'm.', alpha = 0.5)          #Masked particles
        ax.set_title('Time = %.2f Gyr' % (com.time/10000.0).value, fontsize=7)
        ax.set_ylabel('Y-position (kpc)')
        ax.set_xlabel('X-position (kpc)')
    
    fig.tight_layout()
    fig.savefig('snapshots.png')
    plt.close()

    
    
if histograms:
    fig, ax_tmp = plt.subplots(nrows=2, ncols=3, sharey = True)
    
    #looping over seleted snapshots
    for i, t in enumerate(int_times):
        com = CenterOfMass('HighRes/M31_%s.txt' % str(t), 2)
        
        #creating COM object, regenerating parameters and masked parameters for this snapshot
        r, eq_vel, vz, param = get_com_param(com)
        ax = ax_tmp[0][i] if i<3 else ax_tmp[1][i-3]
        tmp = ax.hist(r[mask1], bins = 20) #Histogram using masked particles
        ax.axvline(x=8.29, color = 'k')
        ax.set_title('Time = %.2f Gyr' % (com.time/10000.0).value, fontsize=7)
        ax.set_ylabel('distribution')
        ax.set_xlabel('Radius (kpc)')
        
        data = np.swapaxes(np.array([tmp[0], tmp[1][:-1]]), 0, 1)
        np.savetxt('hist%i.txt' % i, data, fmt = ['%.2f', '%.2f']) #Saving histogram data
    fig.tight_layout()
    fig.savefig('histograms.png')
    plt.close()

    
    
if time_series:

    #selecting 5 particles randomly from the masked particles
    ind   = np.random.randint(0, len(r0), 5)
    v0    = magnitude(eq_vel0, vz0.value)
    ini_r = r0[ind]
    ini_v = v0[ind]
    
    #initializing arrays to hold information for those 5 particles
    snaps = np.arange(0, 550, 10)
    time  = np.zeros([len(snaps)])
    pos   = np.zeros([len(snaps),5])
    vel   = np.zeros([len(snaps),5])
    
    #looping over every 5 snapshots
    for i, s in enumerate(snaps):
        #Generating proper file name
        ilbl = '000' + str(s)
        ilbl = ilbl[-3:]
        folder = 'HighRes/'
        fname = folder + 'M31_' + ilbl + '.txt'
        
        #Flling up initialized arrays with current values
        com = CenterOfMass(fname, 2)
        r, eq_vel, vz, param = get_com_param(com)
        v = magnitude(eq_vel, vz.value)
        print((com.time/10000.0).value, r[ind], v[ind])
        time[i] = (com.time/10000.0).value    #converting Myr to Gyr
        pos[i]  = r[ind]
        vel[i]  = v[ind]
    
    #Saving the generated arrays
    np.savetxt(time, 'time5.txt', fmt = ['%.2f'])
    np.savetxt(pos, 'pos5.txt', fmt = ['%.2f', '%.2f', '%.2f', '%.2f', '%.2f'])
    np.savetxt(vel, 'vel5.txt', fmt = ['%.2f', '%.2f', '%.2f', '%.2f', '%.2f'])
    
    fig, (ax1, ax2) = plt.subplots(1, 2, sharex = True)
    
    for i in range(5):
        ax1.plot(time, pos[:, i], label = 'r0 = %.2f' % ini_r[i])
        ax2.plot(time, vel[:, i], label = 'r0 = %.2f' % ini_v[i])
    
    ax1.legend(loc = 'best')
    ax1.ylabel('Position (kpc)')
    ax2.ylabel('Velocity(km/s)')
    ax2.xlabel('Time (Gyr)')
    ax2.legend(loc = 'best')
    fig.suptitle('Position & Velocity Of 5 particles as a function of time')
    plt.close()
        
      