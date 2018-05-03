import numpy as np
import astropy.units as u
import pdb
from ReadFile import Read


class M33AnalyticOrbit:
    
    def __init__(self, filename):
        
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
            a  = -(G*M)/((R**2 + B**2)**1.5)
            if i=='x':
                ai = a*x
            elif i=='y':
                ai = a*y
            else:
                ai = a*(B/np.sqrt(z**2 + zd**2))*z
                
            return ai
            