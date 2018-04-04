"""
ASTR 400B Homework 4 
Center of Mass Position and Velocity
Rafia Bushra
Date Submitted: 4/4/18
"""

# import modules
import numpy as np
import astropy.units as u
import sys, pdb
sys.path.insert(0, '../HW2/') #Making python search the HW2 directory for ReadFile program
from ReadFile import Read


class CenterOfMass:
    """
   An object class that stores simulated parameters of a galaxy and contains methods that calculates parameters such as mass, COM position and COM Velocity
   """

    def __init__(self, filename, ptype):
        """
        Initializing the object
        
        PARAMETERS
        ----------
        filename : name of the file containing galaxy info. Type = str
        ptype    : Particle type. Possible values are 'Halo', 'Disk' or 'Bulge'. Type = str
        """
        
        # read in the file and particle type
        self.time, self.total, self.data = Read(filename)
            
        #create an array to store indexes of particles of desired Ptype
        self.index = np.where(self.data['type'] == ptype)
    
        # store the mass, positions, velocities of only the particles of the given type
        kms = u.km/u.s  #creating km/s unit for velocity
        self.m  = self.data['m'][self.index]*u.Msun*1e10
        self.x  = self.data['x'][self.index]*u.kpc
        self.y  = self.data['y'][self.index]*u.kpc
        self.z  = self.data['z'][self.index]*u.kpc
        self.vx = self.data['vx'][self.index]*kms
        self.vy = self.data['vy'][self.index]*kms
        self.vz = self.data['vz'][self.index]*kms
      
    def total_mass(self):
        """
        A method that calculates & returns the total mass of this object (galaxy)
        """
        return np.sum(self.m)
        
    def COMdefine(self, x, y, z):
        """
        A method that calculates center of mass of a 3D parameter
        
        PARAMETERS
        ----------
        x, y, z : parameters whose center of mass needs to be found. Type = Array
        
        RETURNS
        -------
        xcom, ycom, zcom: Center of mass of given parameters. Type = Quantity
        """
        M = self.total_mass()
        
        xcom = np.sum(self.m*x)/M
        ycom = np.sum(self.m*y)/M
        zcom = np.sum(self.m*z)/M
            
        return xcom, ycom, zcom
            
    def COM_P(self, delta):
        """
        A method that returns refined COM position within a given tolerance level (delta)
        """
        # initializing the coordinates and calculating position vector magnitude
        xcom, ycom, zcom = self.COMdefine(self.x, self.y, self.z)
        rcom  = np.sqrt((xcom**2) + (ycom**2) + (zcom**2))

        # calculating coordinates in COM's frame of reference
        x_com = self.x - xcom
        y_com = self.y - ycom
        z_com = self.z - zcom

        # array of magnitudes of position vectors in the COM frame
        rnew  = np.sqrt((x_com**2) + (y_com**2) + (z_com**2))

        #Defining rmax: half of max 3d separation of particles from COM position
        rmax = np.amax(rnew)/2

        #An initial guess for the loop
        rcom_diff = 1000   #kpc
        
        #Refining COM position to get a converged value
        while (rcom_diff < delta):
            ind = np.where(rnew < rmax)
            xcom2, ycom2, zcom2 = self.COMdefine(x_com[ind], y_com[ind], z_com[ind])
            rcom2 = np.sqrt((xcom2**2) + (ycom2**2) + (zcom2**2))
            rcom_diff = rcom - rcom2
            rmax = rmax/2
        
        return xcom2, ycom2, zcom2
    
    def COM_V(self, com_p, lim):
        """
        A method that returns COM velocity of particles within a certain limit (lim) of the COM position
        """
        ind   = np.where(com_p < lim)  #selecting particles within the desired limit
        #Masking velocity arrays to select only particles within the limit 
        vx1   = self.vx[ind]
        vy1   = self.vy[ind]
        vz1   = self.vz[ind]
        #Calculating COM velocities of those particles
        com_vx, com_vy, com_vz = self.COMdefine(vx1, vy1, vz1)
        
        return com_vx, com_vy, com_vz





# Printing answers to the terminal
##########################

# Creating Center of mass objects for the MW, M31, M33
MWCOM  = CenterOfMass("../../MW_000.txt", 2)
M31COM = CenterOfMass("../../M31_000.txt", 2)
M33COM = CenterOfMass("../../M33_000.txt", 2)


# Calculating quantities for MW data
print('')
print('')
print('Problem 1:')
print('')
print('MW Parameters:')
print('')

MW_mass = MWCOM.total_mass()
print("Disk Mass: %.2f solar mass" % MW_mass.value)

print('')

mw_x, mw_y, mw_z = MWCOM.COMdefine(MWCOM.x, MWCOM.y, MWCOM.z)
print("COM Position (kpc): %.2f, %.2f, %.2f" % (mw_x.value, mw_y.value, mw_z.value))

print('')
mw_vx, mw_vy, mw_vz = MWCOM.COMdefine(MWCOM.vx, MWCOM.vy, MWCOM.vz)
print("COM Velocity (km/s): %.2f, %.2f, %.2f" % (mw_vx.value, mw_vy.value, mw_vz.value))


# Calculate quantities for M31 data
print('')
print('')
print('M31 Parameters:')
print('')

M31_mass = M31COM.total_mass()
print("Disk Mass: %.2f solar mass" % M31_mass.value)

print('')

m31_x, m31_y, m31_z = M31COM.COMdefine(M31COM.x, M31COM.y, M31COM.z)
print("COM Position (kpc): %.2f, %.2f, %.2f" % (m31_x.value, m31_y.value, m31_z.value))

print('')
m31_vx, m31_vy, m31_vz = M31COM.COMdefine(M31COM.vx, M31COM.vy, M31COM.vz)
print("COM Velocity (km/s): %.2f, %.2f, %.2f" % (m31_vx.value, m31_vy.value, m31_vz.value))


# Calculate quantities for M33 data
print('')
print('')
print('M33 Parameters:')
print('')

M33_mass = M33COM.total_mass()
print("Disk Mass: %.2f solar mass" % M33_mass.value)

print('')

m33_x, m33_y, m33_z = M33COM.COMdefine(M33COM.x, M33COM.y, M33COM.z)
print("COM Position (kpc): %.2f, %.2f, %.2f" % (m33_x.value, m33_y.value, m33_z.value))

print('')
m33_vx, m33_vy, m33_vz = M33COM.COMdefine(M33COM.vx, M33COM.vy, M33COM.vz)
print("COM Velocity (km/s): %.2f, %.2f, %.2f" % (m33_vx.value, m33_vy.value, m33_vz.value))


#Calculating separations
print('')
print('')
print('Problem 2:')
print('')
print('Separation between MW & M31:')
mw_m31_pos = np.sqrt((mw_x - m31_x)**2 + (mw_y - m31_y)**2 + (mw_z - m31_z)**2)
print('Position(kpc) : %.2f' % mw_m31_pos.value)
mw_m31_vel = np.sqrt((mw_vx - m31_vx)**2 + (mw_vy - m31_vy)**2 + (mw_vz - m31_vz)**2)
print('Velocity(km/s) : %.2f' % mw_m31_vel.value)

print('')
print('')
print('Problem 3:')
print('')
print('Separation between M33 & M31:')
m33_m31_pos = np.sqrt((m33_x - m31_x)**2 + (m33_y - m31_y)**2 + (m33_z - m31_z)**2)
print('Position(kpc) : %.2f' % m33_m31_pos.value)
m33_m31_vel = np.sqrt((m33_vx - m31_vx)**2 + (m33_vy - m31_vy)**2 + (m33_vz - m31_vz)**2)
print('Velocity(km/s) : %.2f' % m33_m31_vel.value)


#Answer to part 6 question 4
print('')
print('')
print('Problem 4:')
print('')
print('The iterative process is important because we need accurate convergance of the COM position to study their collision where their COM positions will be very close to each other.')
print('')