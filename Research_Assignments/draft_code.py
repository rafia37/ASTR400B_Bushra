import numpy as np
import matplotlib.pyplot as plt
import sys, pdb
sys.path.insert(0, '../HW2/') #Making python search the HW2 directory for ReadFile program
from ReadFile import Read


MW_data = Read('../../MW_000.txt')[2]
M31_data = Read('../../M31_000.txt')[2]

#use CenterOfMass.py to calculate M31's galactic center at Snapshot 0 (xcen0, ycen0, zcen0)
#Collect indices of disk particles that are within 7 & 9 kpc of (xcen0, ycen0, zcen0)
#Last two steps will give the initial 
#Loop over all desired the snaphots
#-> convert simulation coordinates to the com coordinates using CenterOfMass.py
#-> calculate position and kinematics of desired particles (this answers questions 2 & 3)
#Learn to calculate escape velocity of a galaxy
#if some particles have velocities higher than the escape velocity of M31, they are unbound
#Calculate the fraction of unbound particles (answers question 4)