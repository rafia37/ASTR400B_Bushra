import numpy as np
import astropy.units as u
import pdb
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from astropy.io import ascii




class M33AnalyticOrbit:
    
    def __init__(self, filename):
        
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
        """

        r  = np.sqrt(x**2 + y**2 + z**2)
        a  = -(self.G*M)/(r*((ra+r)**2))
        ai = (a*x) if i=='x' else (a*y) if i=='y' else (a*z)

        return ai


    def MiyamotoNagaiAccel(self, M, rd, x, y, z, i):
        """
        A function that calculates disk acceleration
        """
        zd = self.rd/5
        R  = np.sqrt(x**2 + y**2)
        B  = rd + np.sqrt(z**2 + zd**2)
        a  = -(self.G*M)/((R**2 + B**2)**1.5)
        if i=='x':
            ai = a*x
        elif i=='y':
            ai = a*y
        else:
            ai = a*(B/np.sqrt(z**2 + zd**2))*z

        return ai


    def M31Accel(self, x, y, z, i):

        a_halo  = self.HernquistAccel(self.Mhalo, self.rhalo, x, y, z, i)
        a_bulge = self.HernquistAccel(self.Mbulge, self.rbulge, x, y, z, i)
        a_disk  = self.MiyamotoNagaiAccel(self.Mdisk, self.rd, x, y, z, i)

        a_tot = a_halo + a_bulge + a_disk

        return a_tot


    def LeapFrog(self, dt, x, y, z, vx, vy, vz):

        x_mid = x + (vx*dt*0.5)
        y_mid = y + (vy*dt*0.5)
        z_mid = z + (vz*dt*0.5)

        vx_full = vx + self.M31Accel(x_mid, y_mid, z_mid, 'x')*dt
        vy_full = vy + self.M31Accel(x_mid, y_mid, z_mid, 'y')*dt
        vz_full = vz + self.M31Accel(x_mid, y_mid, z_mid, 'z')*dt

        x_full = x + 0.5*(vx + vx_full)*dt
        y_full = y + 0.5*(vy + vy_full)*dt
        z_full = z + 0.5*(vz + vz_full)*dt

        return x_full, y_full, z_full, vx_full, vy_full, vz_full

    def OrbitIntegrator(self, t0, dt, tmax):

        OrbitParam = np.zeros([int(np.around(tmax/dt))+1, 7])

        x, y, z = self.x, self.y, self.z
        vx, vy, vz = self.vx, self.vy, self.vz
        t = t0

        while t<=tmax:

            row = int(np.around(t/dt))

            OrbitParam[row] = t, x, y, z, vx, vy, vz

            x, y, z, vx, vy, vz = self.LeapFrog(dt, x, y, z, vx, vy, vz)

            t = t + dt

        np.savetxt(self.fname, OrbitParam, header = 't x y z vx vy vz', comments = '#', fmt = ['%.2f', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f'])

        return self.fname

    
def mag(x, y, z = 0):
    return np.sqrt(x**2 + y**2 + z**2)
        
M33Orb = M33AnalyticOrbit('M33_orbit.txt')
M33Orb.OrbitIntegrator(0, 0.1, 10.0)


M33 = ascii.read(M33Orb.fname)
M31_hw6 = ascii.read('Orbit_M31.txt')
M33_hw6 = ascii.read('Orbit_M33.txt')

M33_pos = mag(M33['x'], M33['y'], M33['z'])
M33_vel = mag(M33['vx'], M33['vy'], M33['vz'])

relative_pos = mag(M31_hw6['x'] - M33_hw6['x'], 
                   M31_hw6['y'] - M33_hw6['y'], 
                   M31_hw6['z'] - M33_hw6['z'])

relative_vel = mag(M31_hw6['vx'] - M33_hw6['vx'], 
                   M31_hw6['vy'] - M33_hw6['vy'], 
                   M31_hw6['vz'] - M33_hw6['vz'])


plt.plot(M33['t'], M33_pos)
plt.plot(M31_hw6['t'], relative_pos)
plt.savefig('Orbital_Position.png')
plt.cla()

plt.plot(M33['t'], M33_vel)
plt.plot(M31_hw6['t'], relative_vel)
plt.savefig('Orbital_Velocity.png')

seperation = np.sqrt(M31_hw6['x'])