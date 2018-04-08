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
        xcom, ycom, zcom = COM.COM_P()
        #questions to ask ekta:
        #what should delta be in COM_P
        #what range should I choose for the radius array
        #should calculate com position every time I change radius?
        return masses
        
        