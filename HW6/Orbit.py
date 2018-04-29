import numpy as np
import astropy.units as u
from astropy.io import ascii
from astropy.table import Table
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from ReadFile import Read
from CenterOfMass import CenterOfMass


def OrbitCOM(galaxy, start, end, n):
    
    fileout = "Orbit_%s.txt" % galaxy
    size = int(end/n) + 1
    Orbit = np.zeros([size, 7])
    delta, VolDec = 0.5, 4
    
    for i in np.arange(start, end+1, n):
        
        #Creating filename for desired snapshot
        ilbl = '000' + str(i)
        ilbl = ilbl[-3:]  #making 3 digit numbers
        folder = 'VLowRes/'
        filename = folder + '%s_' % galaxy + ilbl + '.txt'
        
        COM = CenterOfMass(filename, 2)
        
        t = COM.time/1000          #converting Myr to Gyr
        x, y, z = COM.COM_P(delta, VolDec)
        vx, vy, vz = COM.COM_V(x, y, z, 15)
        
        row_ind = int(i/n)
        Orbit[row_ind] = [t.value, x.value, y.value, z.value, vx.value, vy.value, vz.value]
        print('%s i = %i' % (galaxy, i))
        
    np.savetxt(fileout, Orbit, header='t x y z vx vy vz', comments='#', fmt=['%.2f', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f'])
    
    return fileout

def magnitude(x, y, z = 0):
    return np.sqrt(x**2 + y**2 + z**2)

#Uncomment the statements below to generate new orbit file
#orb_mw  = OrbitCOM('MW', 0, 800, 5)
#orb_m31 = OrbitCOM('M31', 0, 800, 5)
#orb_m33 = OrbitCOM('M33', 0, 800, 5)

#Reading in the orbit files
mw_data = ascii.read('Orbit_MW.txt')
m31_data = ascii.read('Orbit_M31.txt')
m33_data = ascii.read('Orbit_M33.txt')

mw_m31_sep = magnitude(mw_data['x'] - m31_data['x'], mw_data['y'] - m31_data['y'], mw_data['z'] - m31_data['z'])

m31_m33_sep = magnitude(m31_data['x'] - m33_data['x'], m31_data['y'] - m33_data['y'], m31_data['z'] - m33_data['z'])

plt.plot(mw_data['t'], mw_m31_sep)
plt.xlabel('Time[Gyr]')
plt.ylabel('Speration [kpc]')
plt.title('Spatial separation between MW and M31')
plt.savefig('MW_M31_Sep.png')
plt.cla()

plt.plot(m31_data['t'], m31_m33_sep)
plt.xlabel('Time[Gyr]')
plt.ylabel('Speration [kpc]')
plt.title('Spatial separation between M31 and M33')
plt.savefig('M31_M33_Sep.png')
