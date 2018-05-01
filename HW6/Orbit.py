"""
ASTR 400B Homework 6
Rafia Bushra
Date Submitted: 4/30/16
"""

import numpy as np
import astropy.units as u
from astropy.io import ascii
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from ReadFile import Read
from CenterOfMass import CenterOfMass


def OrbitCOM(galaxy, start, end, n):
    """
    A function that generates a text file containing orbital information in 7 columns. The file will have a header 't x y z vx vy vz', which serves as the name of each column. Going from left to right, the columns are time, x-position, y-position, z-position, x-velocity, y-velocity, z-velocity.
    
    PARAMETERS
    ----------
    galaxy: Name of the galaxy, e.g. 'MW'. Type = str
    start : Number of the first snapshot to be read in. Type = int
    end   : Number of the last snapshot to be read in. Type = int
    n     : Interval over which snapshots are iterated. Type = int
    
    RETURNS
    -------
    fileout: Name of the generated file. Type = str
    Generates a text file named with the variable fileout in current working directory.
    """
    
    fileout = "Orbit_%s.txt" % galaxy  #Name of the file to be generated
    size = int(end/n) + 1         #Number of rows for the orbit array
    Orbit = np.zeros([size, 7])   #Initializing the array that will hold orbital parameters
    delta, VolDec = 0.5, 4    #Defining parameters required by the COM_P method
    
    #Looping over snapshots and generating orbital parameters at each snapshot
    for i in np.arange(start, end+1, n):
        
        #Creating filename for desired snapshot
        ilbl = '000' + str(i)
        ilbl = ilbl[-3:]  #making 3 digit numbers
        folder = 'VLowRes/'
        filename = folder + '%s_' % galaxy + ilbl + '.txt'
        
        #Creating the center of mass object
        COM = CenterOfMass(filename, 2)
        
        #Generating orbital parameters for the i-th snapshot
        t = COM.time/1000          #converting Myr to Gyr
        x, y, z = COM.COM_P(delta, VolDec)
        vx, vy, vz = COM.COM_V(x, y, z, 15)
        
        #Assigning values to the proper row of orbit array
        row_ind = int(i/n)
        Orbit[row_ind] = [t.value, x.value, y.value, z.value, vx.value, vy.value, vz.value]
        print('%s i = %i' % (galaxy, i))
        
    #Saves the orbit array to a file in current directory    
    np.savetxt(fileout, Orbit, header='t x y z vx vy vz', comments='#', fmt=['%.2f', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f'])
    
    return fileout

def magnitude(x, y, z = 0):
    """
    A function that calculates the magnitude of 2D or 3D vector(s).
    
    PARAMETERS
    ----------
    x, y, z: Components of the vector(s), z being optional. Type = int/array
    
    RETURNS
    -------
    magnitude of the vector(s) with components x, y, z 
    """
    
    return np.sqrt(x**2 + y**2 + z**2)

#Uncomment the statements below to generate new orbit files
#orb_mw  = OrbitCOM('MW', 0, 800, 5)
#orb_m31 = OrbitCOM('M31', 0, 800, 5)
#orb_m33 = OrbitCOM('M33', 0, 800, 5)

#Reading in the orbit files
mw_data = ascii.read('Orbit_MW.txt')
m31_data = ascii.read('Orbit_M31.txt')
m33_data = ascii.read('Orbit_M33.txt')

#Calculating Separations
mw_m31_sep = magnitude(mw_data['x'] - m31_data['x'], mw_data['y'] - m31_data['y'], mw_data['z'] - m31_data['z'])

m31_m33_sep = magnitude(m31_data['x'] - m33_data['x'], m31_data['y'] - m33_data['y'], m31_data['z'] - m33_data['z'])

#Calculating relative velocities
mw_m31_vel = magnitude(mw_data['vx'] - m31_data['vx'], mw_data['vy'] - m31_data['vy'], mw_data['vz'] - m31_data['vz'])

m31_m33_vel = magnitude(m31_data['vx'] - m33_data['vx'], m31_data['vy'] - m33_data['vy'], m31_data['vz'] - m33_data['vz'])


#plotting and saving MW-M31 separation
plt.plot(mw_data['t'], mw_m31_sep)
plt.xlabel('Time[Gyr]')
plt.ylabel('Speration [kpc]')
plt.title('Spatial separation between MW and M31')
plt.savefig('MW_M31_Sep.png')
plt.cla()

#plotting and saving M31-M33 separation
plt.plot(m31_data['t'], m31_m33_sep)
plt.xlabel('Time[Gyr]')
plt.ylabel('Speration [kpc]')
plt.title('Spatial separation between M31 and M33')
plt.savefig('M31_M33_Sep.png')
plt.cla()

#plotting and saving MW-M31 relative velocity
plt.plot(mw_data['t'], mw_m31_vel)
plt.xlabel('Time[Gyr]')
plt.ylabel('Relative Velocity [km/s]')
plt.title('Relative velocity between MW and M33')
plt.savefig('MW_M31_Vel.png')
plt.cla()

#plotting and saving M31-M33 relative velocity
plt.plot(m31_data['t'], m31_m33_vel)
plt.xlabel('Time[Gyr]')
plt.ylabel('Relative Velocity [km/s]')
plt.title('Relative velocity between M31 and M33')
plt.savefig('M31_M33_Vel.png')


#Printing part 4 answers to terminal
print('')
print('Answer 1:')
print('MW and M31 will experience 3 close encounters in the future.')
print('')
print('Answer 2:')
print('Relative velocity between MW and M31 increases as separation between them decreases with time')
print('')
print('Answer 3:')
print('MW and M31 will merge around 6.5 Gyr from now. Around that time, M31 and M33 experiences a close encounter and has more frequent and closer encounters in the future. Which means M33\'s orbit slowly decays into the MW-M31 merger remnant.')
print('')