import numpy as np
import astropy.units as u


def Read(filename):
    """
    A function that generates, time, number of particles and particle information from a given file.
    
    PARAMETERS
    filename : Name of the file to generate information from. Type = str
    
    Returns
    time : Time in 10 Myr. Type = float
    N    : Number of particles. Type = int
    data : type, mass, position(x, y, z), velocity(vx, vy, vz) of particle. Type = Array
    """
    file = open(filename,'r')

    #Storing time in units of 10Myr
    line1 = file.readline()
    label, value = line1.split()
    time = float(value)*u.Myr
    
    
    #Storing the number of particles
    line2 = file.readline()
    label, value = line2.split()
    N = int(value)
    
    file.close()   #closes file
    
    data = np.genfromtxt(filename, dtype = None, names = True, skip_header = 3) #generating particle info from given file
    
    return time, N, data


