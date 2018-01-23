import numpy as np
from ReadFile import Read
import astropy.units as u
import sys, pdb

def ParticleInfo(filename, ind):
    
    """
    A function that provides information about requested particle.
    
    PARAMETERS
    filename : name of the file containing galaxy info. Type = str
    ind      : index of the particle for which information is requested. Type = int
    
    RETURNS
    x, y, z    : Distance of partcle from center in kpc. Type = Quantity
    vx, vy, vz : Velocity of particle in km/s. Type = Quantity
    mass       : Mass of particle in solar mass. Type = Quantity
    """
    
    time, N, data = Read(filename)  #acquiring particle data from MW_000.txt file
    
    #collecting position information and setting correct unit (Kpc)
    x = data['x'][ind]*u.kpc
    y = data['y'][ind]*u.kpc
    z = data['z'][ind]*u.kpc
    
    #collecting velocity information and setting correct unit (Km/s)
    kms = u.km/u.s
    vx = data['vx'][ind]*kms
    vy = data['vy'][ind]*kms
    vz = data['vz'][ind]*kms
    
    #collecting mass information and setting correct unit (M_sun)
    mass = data['m'][ind]*u.M_sun
    
    return x, y, z, vx, vy, vz, mass

#What this script will do when it's called from terminal
if __name__ == "__main__":
    argv = sys.argv
    
    filename = argv[1]         #filename is the first argument provided in terminal
    ind      = int(argv[2])    #Index of particle is the second argument provided in terminal
    x, y, z, vx, vy, vz, mass = ParticleInfo(filename, ind)   #Collecting information about requested particle
    
    print('')    #For readability in terminal
    print('Information about the %i th particle:' % (ind+1))
    print('')
    print('Distance in (x,y,z) kpc: (%.3f, %.3f, %.3f)' % (x.value, y.value, z.value))
    print('')
    print('Distance in (x,y,z) ly: (%.3f, %.3f, %.3f)' % (x.to(u.lyr).value, y.to(u.lyr).value, z.to(u.lyr).value))
    print('')
    print('Velocity in (vx,vy,vz) km/s: (%.3f, %.3f, %.3f)' % (vx.value, vy.value, vz.value))
    print('')
    print('Mass in solar mass: %f' % mass.value)
    print('')
    
    