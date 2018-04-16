"""
ASTR 400B Homework 5
Mass distribution and roation curve
Rafia Bushra
Date Submitted: 4/10/18
"""

# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G
from scipy.optimize import curve_fit
import sys, pdb
import matplotlib.pyplot as plt
from ReadFile import Read
from CenterOfMass import CenterOfMass


class MassProfile:
    
    def __init__(self, galaxy, snap):
        
        #Creating filename for desired snapshot
        ilbl = '000' + str(snap)
        ilbl = ilbl[-3:]
        self.filename = '%s_' % galaxy + ilbl + '.txt'
        
        self.time, self.total, self.data = Read(self.filename)
        
        #Stroring galaxy properties as global parameter
        self.gname = galaxy
        self.m  = self.data['m']*u.Msun*1e10
        self.x  = self.data['x']*u.kpc
        self.y  = self.data['y']*u.kpc
        self.z  = self.data['z']*u.kpc
        
    def MassEnclosed(self, ptype, radii):
        
        COM = CenterOfMass(self.filename, ptype)
        
        #center of mass position of the galaxy
        xcom, ycom, zcom = COM.COM_P(1) #tolerance level of 1kpc
        
        # calculating coordinates of each particle in COM's frame of reference
        xnew = COM.x - xcom
        ynew = COM.y - ycom
        znew = COM.z - zcom
        rnew = np.sqrt((xnew**2) + (ynew**2) + (znew**2)) #position vecotr magnitudes

        masses = np.zeros(len(radii)) 
        for i, r in enumerate(radii):
            ind = np.where(rnew.value<r)[0]
            mass = np.sum(COM.m[ind])
            masses[i] = mass.value
        return masses*u.Msun
    
    def MassEnclosedTotal(self, radii):
        halo_mass  = self.MassEnclosed(1, radii)
        disk_mass  = self.MassEnclosed(2, radii)
        bluge_mass = self.MassEnclosed(3, radii)
        total_mass = halo_mass + disk_mass + bulge_mass
        return total_mass
    
    def HernquistMass(self, r, a, M_halo):
        M = (M_halo*(r**2))/((a+r)**2)
        return M
    
    def CircularVelocity(self, ptype, radii):
        G = G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        M = self.MassEnclosed(ptype, radii)
        v = np.sqrt((G*M)/radii)
        return v
    
    def CircularVelocityTotal(self, radii):
        M = self.MassEnclosedTotal(radii)
        v = np.sqrt((G*M)/radii)
        return v
    
    def HernquistVCirc(self, r, a, M_halo):
        M = self.HernquistMass(r, a, M_halo)
        v = np.sqrt((G*M)/r)
        return v
    

# Plotting Interesting quantities
##################################

if __name__ == '__main__':
    MW_MP = MassProfile('MW', 0)
    
    radii      = np.linspace(0.1, 30, 50)
    halo_mass  = MW_MP.MassEnclosed(1, radii)
    disk_mass  = MW_MP.MassEnclosed(2, radii)
    bulge_mass = MW_MP.MassEnclosed(3, radii)
    total_mass = MW_MP.MassEnclosedTotal(radii)
    
    """
    plt.semilogy(radii, halo_mass, label = 'Halo')
    plt.semilogy(radii, disk_mass, label = 'Disk')
    plt.semilogy(radii, bulge_mass, label = 'Bulge')
    plt.semilogy(radii, total_mass, label = 'Total')
    plt.title('Mass distribution of Milky Way as a function of radius')
    plt.xlabel('Radius (kpc)')
    plt.ylabel('log(Mass)')
    plt.legend(loc = 'best')
    plt.show()
    """
    
    M_halo = MW_MP.MassEnclosed(1, [15]) #15kpc radius
    plt.semilogy(radii, halo_mass, label = 'Dark Matter Profile')
    for a in np.linspace(0.01, 1, 9):  
        plt.semilogy(radii, MW_MP.HernquistMass(radii, a, M_halo), label = 'a = %.2f' % a)
    plt.legend(loc = 'best')
    plt.xlabel('Radius (kpc)')
    plt.ylabel('Mass')
    plt.title('Fitting Herquist model to MW dark matter profile')
    plt.show()
    
    