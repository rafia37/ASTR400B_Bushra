import numpy as np
import astropy.units as u
import sys, pdb
sys.path.insert(0, '../HW2/') #Making python search the HW2 directory for ReadFile program
from ReadFile import Read


def ComponentMass(filename, ptype):
    """
    A function that generates masss of each component of a galaxy
    
    PARAMETERS
    ----------
    filename : Name of the file containing galaxy info. Type = str
    ptype    : Particle type. Possible values are 'Halo', 'Disk' or 'Bulge'. Type = str
    
    RETURNS
    -------
    CompMass : Component mass. Mass of halo, disk or bulge as indicated by ptype in solar mass. Type = Quantity
    """
    data     = Read(filename)[2]
    CompType = 1 if ptype=='Halo' else 2 if ptype=='Disk' else 3 
    ind      = np.where(data['type'] == CompType)
    CompMass = np.sum(data['m'][ind]*1e10*u.M_sun)
    return CompMass

#What this script will do when it's called from terminal
if __name__ == "__main__":
    argv = sys.argv
    
    filename = argv[1]     #filename is the first argument provided in terminal
    ptype    = argv[2]     #secondparticle type is the  argument provided in terminal
    
    print('') # For readability
    print('%s mass of this galaxy is %f solar mass' % (ptype, ComponentMass(filename, ptype).value))
    print('')