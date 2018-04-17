"""
ASTR 400B Homework 5
Mass distribution and roation curve
Rafia Bushra
Date Submitted: 4/16/18
"""

# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G
G = G.to(u.kpc*u.km**2/u.s**2/u.Msun)
from scipy.optimize import curve_fit
import sys, pdb
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from ReadFile import Read
from CenterOfMass import CenterOfMass


class MassProfile:
    """
    An object class that calculates mass and circular velocity for components, hernquist model and total
    """
    
    def __init__(self, galaxy, snap):
        """
        Initializing the object
        
        PARAMETERS
        ----------
        galaxy : name of galaxy. e.g. "MW", "M31" etc. Type = str
        ptype  : Particle type. Possible values are 'Halo', 'Disk' or 'Bulge'. Type = str
        """
        
        #Creating filename for desired snapshot
        ilbl = '000' + str(snap)
        ilbl = ilbl[-3:]  #making 3 digit numbers
        self.filename = '%s_' % galaxy + ilbl + '.txt'
        
        self.time, self.total, self.data = Read(self.filename) #Reading in simulation file
        
        #Stroring galaxy properties as global parameters
        self.gname = galaxy
        self.m  = self.data['m']*u.Msun*1e10
        self.x  = self.data['x']*u.kpc
        self.y  = self.data['y']*u.kpc
        self.z  = self.data['z']*u.kpc
        
    def MassEnclosed(self, ptype, radii):
        """
        A method that calculates the mass enclosed within a given radius of the galaxy for a given component
        
        PARAMETERS
        ----------
        ptype  : Particle type. Possible values are 'Halo', 'Disk' or 'Bulge'. Type = str
        radii  : Radii within which total mass is to be computed. Type = Array
        
        RETURNS
        -------
        masses : masses at each radii. Type = Array
        """
        
        #creating exception for M33's non-existent bulge
        if (self.gname=='M33') & (ptype==3):
            return np.zeros(len(radii))*u.Msun
        
        COM = CenterOfMass(self.filename, ptype)
        
        #center of mass position of the galaxy
        xcom, ycom, zcom = COM.COM_P(1) #tolerance level of 1kpc
        
        # calculating coordinates of each particle in COM's frame of reference
        xnew = COM.x - xcom
        ynew = COM.y - ycom
        znew = COM.z - zcom
        rnew = np.sqrt((xnew**2) + (ynew**2) + (znew**2)) #position vecotr magnitudes

        masses = np.zeros(len(radii)) #initializing mass array to be returned
        #calculating mass at each radius
        for i, r in enumerate(radii):
            ind = np.where(rnew.value<r)[0] #masking for particles within the radius
            mass = np.sum(COM.m[ind])
            masses[i] = mass.value
        return masses*u.Msun
    
    def MassEnclosedTotal(self, radii):
        """
        A method that calculates the total mass enclosed within a given radius of the galaxy
        
        PARAMETERS
        ----------
        radii  : Radii within which total mass is to be computed. Type = Array
        
        RETURNS
        -------
        total_mass : total mass at each radii. Type = Array
        """
        
        halo_mass  = self.MassEnclosed(1, radii)
        disk_mass  = self.MassEnclosed(2, radii)
        bluge_mass = self.MassEnclosed(3, radii)
        total_mass = halo_mass + disk_mass + bulge_mass
        return total_mass
    
    def HernquistMass(self, r, M_halo, a = False, M = False):
        """
        A method that calculates a theoretical mass profile
        
        PARAMETERS
        ----------
        r      : radius within which mass is to be computed. Type = float/int
        M_halo : Total halo mass of the galaxy. Type = float/int
        a      : constant parameter. Required if you want to calculate M. Type = float/int
        M      : Hernquist mass. Required if you want to calculate a. Type = float/int
        
        RETURNS
        -------
        M : Hernquist mass. Will be returned if a is provided. Type = float/int
        a : constant parameter. Will be returned if M is provided. Type = float/int
        """
        
        #A choice to calculate either M or a since these are two unknowns that depend on each other
        if M == False:
            M = (M_halo*(r**2))/((a+r)**2)
            return M
        else:
            a = np.sqrt((M_halo*(r**2))/M) - r
            return a
        
    
    def CircularVelocity(self, ptype, radii):
        """
        A method that calculates a theoretical mass profile
        
        PARAMETERS
        ----------
        ptype  : Particle type. Possible values are 'Halo', 'Disk' or 'Bulge'. Type = str
        radii  : Radii at which circular velocity is to be computed. Type = Array
        
        RETURNS
        -------
        v : Circular velocity computed using component mass. Type = float/int
        """
        
        M = self.MassEnclosed(ptype, radii)
        v = np.sqrt((G*M)/radii)
        return v
    
    def CircularVelocityTotal(self, radii):
        """
        A method that calculates a theoretical mass profile
        
        PARAMETERS
        ----------
        radii  : Radii at which circular velocity is to be computed. Type = Array
        
        RETURNS
        -------
        v : Circular velocity computed using totsl mass. Type = float/int
        """
        
        M = self.MassEnclosedTotal(radii)
        v = np.sqrt((G*M)/radii)
        return v
    
    def HernquistVCirc(self, r, a, M_halo):
        """
        A method that calculates a theoretical mass profile
        
        PARAMETERS
        ----------
        r      : radius at which circular velocity is to be computed. Type = float/int
        M_halo : Total halo mass of the galaxy. Type = float/int
        a      : constant parameter. Type = float/int
        
        RETURNS
        -------
        v : Circular velocity computed using Hernquist mass. Type = float/int
        """
        
        M = self.HernquistMass(r, M_halo, a)
        v = np.sqrt((G*M)/r)
        return v
    

# Plotting Interesting quantities
##################################

if __name__ == '__main__':
    gname = ['MW', 'M31', 'M33']
    radii = np.linspace(0.1, 30, 70)
    
    #looping over each galaxy
    for g in gname:
        
        MP = MassProfile(g, 0)
        
        #Hernquist Mass Function Parameters
        ind    = np.where(MP.data['type']==1)[0]
        M_halo = np.sum(MP.m[ind])
        a      = MP.HernquistMass(30, M_halo, M = MP.MassEnclosed(1, [30])) 
        
        #Calculating, component, total and Hernquist masses of each galaxy
        halo_mass  = MP.MassEnclosed(1, radii)
        disk_mass  = MP.MassEnclosed(2, radii)
        bulge_mass = MP.MassEnclosed(3, radii)
        total_mass = MP.MassEnclosedTotal(radii)
        hernq_mass = MP.HernquistMass(radii, M_halo, a)
        
        #plotting mass profiles
        plt.semilogy(radii, halo_mass, label = 'Halo Mass')
        plt.semilogy(radii, disk_mass, label = 'Disk Mass')
        plt.semilogy(radii, bulge_mass, label = 'Bulge Mass')
        plt.semilogy(radii, total_mass, label = 'Total Mass')
        plt.semilogy(radii, hernq_mass, 'k--', label = 'Hernquist Mass Profile: a = %.2f' % a)
        plt.title('%s Mass Distribution As A Function Of Radius' % g)
        plt.xlabel('Radius (kpc)')
        plt.ylabel(r'Mass\ [M_\odot]')
        plt.legend(loc = 'best')
        plt.savefig('%s_profile.png' % g)
        plt.clf()
        
        #Calculating, component, total and Hernquist cicular velocities of each galaxy
        halo_vel = MP.CircularVelocity(1, radii)
        disk_vel = MP.CircularVelocity(2, radii)
        bulge_vel = MP.CircularVelocity(3, radii)
        total_vel = MP.CircularVelocityTotal(radii)
        hernq_vel = MP.HernquistVCirc(radii, a, M_halo)
        
        #plotting rotation curves
        plt.semilogy(radii, halo_vel, label = 'Halo')
        plt.semilogy(radii, disk_vel, label = 'Disk')
        plt.semilogy(radii, bulge_vel, label = 'Bulge')
        plt.semilogy(radii, total_vel, label = 'Total')
        plt.semilogy(radii, hernq_vel, 'k--', label = 'Hernquist: a = %.2f' % a)
        plt.title('%s Rotation Curve' % g)
        plt.xlabel('Radius (kpc)')
        plt.ylabel('Circular Velocity (km/s)')
        plt.legend(loc = 'best')
        plt.savefig('%s_rotation_curve.png' % g)
        plt.clf()