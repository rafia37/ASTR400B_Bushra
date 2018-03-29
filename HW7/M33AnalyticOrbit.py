import numpy as np
import astropy.units as u
import sys, pdb
sys.path.insert(0, '../HW2/') #Making python search the HW2 directory for ReadFile program
from ReadFile import Read


class M33AnalyticOrbit:
    
    def __init__(self, filename):
        