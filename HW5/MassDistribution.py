"""
ASTR 400B Homework 5
Mass distribution and roation curve
Rafia Bushra
Date Submitted: 4/10/18
"""

# import modules
import numpy as np
import astropy.units as u
import sys, pdb
import matplotlib.pyplot as plt
sys.path.insert(0, '../HW2/') #Making python search the HW2 directory 
sys.path.insert(1, '../HW4/') #Making python search the HW4 directory 
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
        self.m  = self.data['m']
        self.x  = self.data['x']*u.kpc
        self.y  = self.data['y']*u.kpc
        self.z  = self.data['z']*u.kpc
        
    def MassEnclosed(self, ptype, radii):
        
        COM = CenterOfMass(self.filename, ptype)
        
        #center of mass position of the galaxy
        xcom, ycom, zcom = COM.COM_P(1) #tolerance level of 1kpc
        
        # calculating coordinates of each particle in COM's frame of reference
        x_com = COM.x - xcom
        y_com = COM.y - ycom
        z_com = COM.z - zcom
        r_com = np.sqrt((x_com**2) + (y_com**2) + (z_com**2)) #position vecotr magnitudes

        masses = np.zeros(len(radii)) 
        for i, r in enumerate(radii):
            ind = np.where((self.data['type'] == ptype) & (rcom<r))
            mass = np.sum(self.m[ind])
            masses[i] = mass
        return masses
    
    def MassEnclosedTotal(self, radii):
        halo_mass  = self.MassEnclosed(1, radii)
        disk_mass  = self.MassEnclosed(2, radii)
        bluge_mass = self.MassEnclosed(3, radii)
        total_mass = halo_mass + disk_mass + bulge_mass
        return total_mass
        
        