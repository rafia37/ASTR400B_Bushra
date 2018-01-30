"""
A script to generate a table containing masss distribution of the local group
"""
import numpy as np
from astropy.table import Table
from astropy.io import ascii
import astropy.units as u
import sys, pdb
sys.path.insert(0, '../HW2/')
from ReadFile import Read
from GalaxyMass import ComponentMass


#File names
mw  = '../../MW_000.txt'
m31 = '../../M31_000.txt'
m33 = '../../M33_000.txt'

GalFiles = [mw, m31, m33]
GalNames = ['Milky Way', 'M31', 'M33']

MassDist = Table(names=('Galaxy Name', 'Halo Mass [1e12 M_sun]', 'Disk Mass [1e12 M_sun]', 'Bulge Mass [1e12 M_sun]', 'Total Mass [1e12 M_sun]', 'f_bar'), dtype = ('S15', 'f8', 'f8', 'f8', 'f8', 'f8'))
# Calculating component masses of each galaxy and appending them to the table 
for i, (fname, name) in enumerate(zip(GalFiles, GalNames)):
    halo  = ComponentMass(fname, 'Halo')/1e12
    disk  = ComponentMass(fname, 'Disk')/1e12
    bulge = ComponentMass(fname, 'Bulge')/1e12
    total = halo + disk + bulge
    f_bar = (disk + bulge)/total
    MassDist.add_row([name, halo, disk, bulge, total, f_bar])

    
# Calculating same variables for the local group
lg_halo  = np.sum(MassDist.columns[1])
lg_disk  = np.sum(MassDist.columns[2])
lg_bulge = np.sum(MassDist.columns[3])
lg_total = np.sum(MassDist.columns[4])
lg_f_bar = (lg_disk + lg_bulge)/lg_total
MassDist.add_row(['Local Group', lg_halo, lg_disk, lg_bulge, lg_total, lg_f_bar])

print('') #Readability purposes
print(MassDist)
print('')

ascii.write(MassDist, 'Table.csv', overwrite = True)