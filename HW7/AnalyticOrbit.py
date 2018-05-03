"""
ASTR 400B Homework 7 
M33 analytic orbit
Rafia Bushra
Date Submitted: 5/2/18
"""


import numpy as np
import astropy.units as u
import pdb
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from astropy.io import ascii




class M33AnalyticOrbit:
    
    def __init__(self, filename):
        
        #assigning the filename as an object property
        self.fname = filename
        
        #Initializing COM position and COM velocity vectore of M33 with respect to M31 from HW4 results
        self.x, self.y, self.z    = -98.53, -120.01, -127.76
        self.vx, self.vy, self.vz = -29.13, 173.81, 93.53
        
        #Radius and mass of M31'a disk
        self.rd    = 5
        self.Mdisk = 0.12e12
        
        #Radius and mass of M31'a bulge
        self.rbulge = 1
        self.Mbulge = 0.019e12
        
        #Radius and mass of M31'a halo
        self.rhalo = 62
        self.Mhalo = 1.921e12
        
        #Gravitational constant in units of kpc**3/M_sun/Gyr
        self.G = 4.498768e-6
        
        
    def HernquistAccel(self, M, ra, x, y, z, i):
        """
        A function that calculates bulge or halo acceleration.
        
        PARAMETERS
        ----------
        M       : Halo mass or bulge mass. Type = float
        ra      : Halo radius or bulge radius. Type = float
        x, y, z : COM coordinates, Type = float
        i       : Component of acceleration. Could be x, y or z. Type = str
        
        RETURNS
        -------
        ai : i component of acceleration. Type = float
        """

        r  = np.sqrt(x**2 + y**2 + z**2)
        a  = -(self.G*M)/(r*((ra+r)**2))  #Acceleration derived from a Herquist profile
        ai = (a*x) if i=='x' else (a*y) if i=='y' else (a*z) #choosing the desired component

        return ai


    def MiyamotoNagaiAccel(self, M, rd, x, y, z, i):
        """
        A function that calculates disk acceleration
        
        PARAMETERS
        ----------
        M       : Disk mass. Type = float
        rd      : Disk radius. Type = float
        x, y, z : COM coordinates, Type = float
        i       : Component of acceleration. Could be x, y or z. Type = str
        
        RETURNS
        -------
        ai : i component of acceleration. Type = float
        """
        
        zd = self.rd/5      #Disk scale height
        R  = np.sqrt(x**2 + y**2)        
        B  = rd + np.sqrt(z**2 + zd**2)  
        a  = -(self.G*M)/((R**2 + B**2)**1.5)   #acceleration derived using Miyamoto-Nagai 1975 profile
        if i=='x':
            ai = a*x
        elif i=='y':
            ai = a*y
        else:
            ai = a*(B/np.sqrt(z**2 + zd**2))*z  #z-component slightly different from x & y components

        return ai


    def M31Accel(self, x, y, z, i):
        """
        A function that calculates total acceleration for a given component
        
        PARAMETERS
        ----------
        x, y, z : COM coordinates, Type = float
        i       : Component of acceleration. Could be x, y or z. Type = str
        
        RETURNS
        -------
        a_tot : i component of total acceleration. Type = float
        """

        a_halo  = self.HernquistAccel(self.Mhalo, self.rhalo, x, y, z, i)  #i-th component of halo acceleration
        a_bulge = self.HernquistAccel(self.Mbulge, self.rbulge, x, y, z, i)  #i-th component of bulge acceleration
        a_disk  = self.MiyamotoNagaiAccel(self.Mdisk, self.rd, x, y, z, i)  #i-th component of disk acceleration

        a_tot = a_halo + a_bulge + a_disk  #i component of total acceleration

        return a_tot


    def LeapFrog(self, dt, x, y, z, vx, vy, vz):
        """
        An integration scheme. Calculates x, y, z, vx, vy, vz after a small time interval dt. Needs to be looped over dt.
        """
        
        #Calculating position coordinates at middle of time step
        x_mid = x + (vx*dt*0.5)
        y_mid = y + (vy*dt*0.5)
        z_mid = z + (vz*dt*0.5)
        
        #Calcuating COM velocity after time dt (full step)
        vx_full = vx + self.M31Accel(x_mid, y_mid, z_mid, 'x')*dt
        vy_full = vy + self.M31Accel(x_mid, y_mid, z_mid, 'y')*dt
        vz_full = vz + self.M31Accel(x_mid, y_mid, z_mid, 'z')*dt

        #Calcuating COM position after time dt (full step)
        x_full = x + 0.5*(vx + vx_full)*dt
        y_full = y + 0.5*(vy + vy_full)*dt
        z_full = z + 0.5*(vz + vz_full)*dt

        return x_full, y_full, z_full, vx_full, vy_full, vz_full

    def OrbitIntegrator(self, t0, dt, tmax):
        """
        A function that loops over LeapFrog function to calculate position and velocity coordinates at regular time intervals dt until a desired time tmax. Parameters along with time are stored in an array and saved in current working directory. 
        """
        
        #initializing an array to store the orbital parameters at each time point
        OrbitParam = np.zeros([int(np.around(tmax/dt))+1, 7])

        #Initializing parameters
        x, y, z = self.x, self.y, self.z
        vx, vy, vz = self.vx, self.vy, self.vz
        t = t0

        #looping over time
        while t<=tmax:
            
            #row index
            row = int(np.around(t/dt))
            
            #filling each row with current parameters
            OrbitParam[row] = t, x, y, z, vx, vy, vz

            #updating parameters
            x, y, z, vx, vy, vz = self.LeapFrog(dt, x, y, z, vx, vy, vz)
            t = t + dt

        np.savetxt(self.fname, OrbitParam, header = 't x y z vx vy vz', comments = '#', fmt = ['%.2f', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f'])

        return self.fname

    
def mag(x, y, z = 0):
    """
    Calculates magnitude of 2D or 3D vectors
    """
    return np.sqrt(x**2 + y**2 + z**2)

#Creating the object instance
M33Orb = M33AnalyticOrbit('M33_orbit.txt')

#uncomment next line if you want to regenerate the orbit file
#M33Orb.OrbitIntegrator(0, 0.1, 10.0)

#Reading in orbit files
M33 = ascii.read(M33Orb.fname)
M31_hw6 = ascii.read('Orbit_M31.txt')
M33_hw6 = ascii.read('Orbit_M33.txt')


#Analytic orbital position and velocity of M33
M33_pos = mag(M33['x'], M33['y'], M33['z'])
M33_vel = mag(M33['vx'], M33['vy'], M33['vz'])

#M31-M33 relative position and velocity
relative_pos = mag(M31_hw6['x'] - M33_hw6['x'], 
                   M31_hw6['y'] - M33_hw6['y'], 
                   M31_hw6['z'] - M33_hw6['z'])
relative_vel = mag(M31_hw6['vx'] - M33_hw6['vx'], 
                   M31_hw6['vy'] - M33_hw6['vy'], 
                   M31_hw6['vz'] - M33_hw6['vz'])

#M33 orbit plot
plt.plot(M33['t'], M33_pos, label = 'M33 orbit')
plt.plot(M31_hw6['t'], relative_pos, label = 'M31-M33 separation')
plt.xlabel('Time (Gyr)')
plt.ylabel('Distance (kpc)')
plt.title('M33 Orbit')
plt.legend(loc = 'best')
plt.savefig('Orbital_Position.png')
plt.cla()

#M33 velocity plot
plt.plot(M33['t'], M33_vel, label = 'M33 Orbital Velocity')
plt.plot(M31_hw6['t'], relative_vel, label = 'M31-M33 Relative Velocity')
plt.xlabel('Time (Gyr)')
plt.ylabel('Velocity (km/s)')
plt.title('M33 Velocity')
plt.legend(loc = 'best')
plt.savefig('Orbital_Velocity.png')


#printing answers to terminal
print('')
print('Answer 2:')
print('The plots are consistent only until the first close encounter.')
print('')
print('Answer 3:')
print('The integrator was made assuming M33 is a point mass, however it\'s mass distribution will have an impact on the integrated orbit')
print('')
print('Answer 4:')
print('We could treat MW and M31 like a single system by finding resultant position and velocity vectors. Then instead of calculation M33\'s orbit and velocity with respect to M31, we could calculate it with respect to the MW-M31 system.')
print('')