import numpy as np
import astropy.units as u
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

orb_mw  = OrbitCOM('MW', 0, 800, 5)
orb_m31 = OrbitCOM('M31', 0, 800, 5)
orb_m33 = OrbitCOM('M33', 0, 800, 5)
