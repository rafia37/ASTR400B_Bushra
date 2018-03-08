# Homework 4 
# Center of Mass Position and Velocity

# import modules
import numpy as np
import astropy.units as u
import sys, pdb
sys.path.insert(0, '../HW2/') #Making python search the HW2 directory for ReadFile program
from ReadFile import Read


class CenterOfMass:
   
    def __init__(self, filename, ptype):
        # read in the file and particle type
        self.time, self.total, self.data = Read(filename)
            
        #create an array to store indexes of particles of desired Ptype
        self.index = np.where(self.data['type'] == ptype)
    
        # store the mass, positions, velocities of only the particles of the given type
        self.m  = self.data['m'][self.index]
        self.x  = self.data['x'][self.index]
        self.y  = self.data['y'][self.index]
        self.z  = self.data['z'][self.index]
        self.vx = self.data['vx'][self.index]
        self.vx = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]
      
    def total_mass(self):
        return np.sum(self.m)*u.Msun*1e10
        
    def COMdefine(self, i, j, k):
        
        M = self.total_mass()
        
        xcom = np.sum(self.m*i)/M
        ycom = np.sum(self.m*j)/M
        zcom = np.sum(self.m*k)/M
            
        return xcom, ycom, zcom
            
    def COM_P(self, delta):
        
        # initializing the coordinates and calculating position vector magnitude
        xcom, ycom, zcom = self.COMdefine(self.x, self.y, self.z)
        rcom  = np.sqrt((xcom**2) + (ycom**2) + (zcom**2))

        # calculating coordinates in COM's frame of reference
        x_com = self.x - xcom
        y_com = self.y - ycom
        z_com = self.z - zcom

        # array of magnitudes of position vectors in the COM frame
        rnew  = np.sqrt((x_com**2) + (y_com**2) + (z_com**2))

        rmax = np.amax(rnew)/2

        rcom_diff = 1000   #kpc
        while (rcom_diff < delta):
            ind = np.where(rnew < rmax)
            xcom2, ycom2, zcom2 = self.COMdefine(x_com[ind], y_com[ind], z_com[ind])
            rcom2 = np.sqrt((xcom**2) + (ycom**2) + (zcom**2))
            rcom_diff = rcom - rcom2
            rmax = rmax/2
         
        return rcom2
    
    def COM_V(self, com_p, lim):
        ind   = np.where(com_p < lim)
        vx1   = self.vx[ind]
        vy1   = self.vy[ind]
        vz1   = self.vz[ind]
        com_vx, com_vy, com_vz = self.COMdefine(vx1, vy1, vz1)
        
        return com_vx, com_vy, com_vz





# EXAMPLE OF USING A CLASS
##########################

# Create a Center of mass object for the MW
MWCOM = CenterOfMass("MW_000.txt", 2)

# Calculate quantities for MW data
MW_mass = MWCOM.total_mass()
print("MW Disk Mass:", MW_mass)
