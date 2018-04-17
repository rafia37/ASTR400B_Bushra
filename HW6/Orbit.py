import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from Readfile import Read
from CenterOfMass import CenterOfMass


def OrbitCOM(self, galaxy, start, end, n):
    
    fileout = "Orbit_%s.txt" % galaxy
    size = int(end/n) + 1
    Orbit = np.zeros(size)
    
        
        