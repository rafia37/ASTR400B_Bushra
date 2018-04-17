import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from Readfile import Read
from CenterOfMass import CenterOfMass


def OrbitCOM(self, galaxy, start, end, n):
    
    fileout = "Orbit_%s.txt" % galaxy
    size = int(end/n) + 1
    Orbit = np.zeros([size, 7])
    delta, VolDec = 0.5, 4
    
    for i in np.arange(start, end+1, n):
        
        #Creating filename for desired snapshot
        ilbl = '000' + str(i)
        ilbl = ilbl[-3:]  #making 3 digit numbers
        filename = '%s_' % galaxy + ilbl + '.txt'
        
        COM = CenterOfMass(filename, 2)
        
    
        
        