#!/usr/env python3 

"""Examples for the Migdal module. P. Cox (2022)"""

from math import log
import numpy as np

from migdal import Migdal

#######################################################

# Initialisation
Ne = Migdal('Ne')


#######################################################
# Orbitals and binding energies
#######################################################
print('\nOrbital binding energies:')
for (orbital, energy) in Ne.orbitals:
    print('{:4}: {:.3e} keV'.format(orbital, energy))


#######################################################
# Exclusive differential ionisation probabilities
#######################################################

print('\nDifferential probabilities:')

Ne.load_probabilities(double=True)

# Array of points to evaluate
pts1 = np.array([[log(1.3), log(0.005)]]) # (ln(En), ln(v))
pts2 = np.array([[log(1.3), log(0.05), log(0.005)]]) # (ln(En1), ln(En2), ln(v))

# Single ionisation
print("dpI1[1s] =", Ne.dpI1(pts1, '1s'))

# Double ionisation
print("dpI2[1s2p] =", Ne.dpI2(pts2, '1s2p'))

# Note that these functions return dP/dE, not dP/dlnE

# For dark matter applications (i.e. recoil velocities v/alpha << 1),
# there is a logarthmic velocity grid available
#
# Ne.load_probabilities(dark_matter=True, double=False)


#######################################################
# Exclusive integrated ionisation probabilities
#######################################################

print('\nIntegrated probabilities (Eth = 1.0 keV):')

Ne.load_probabilities(double=True, integrated=True, e_threshold=1.0)

# Array of velocities to evaluate
pts = np.array([log(0.005)]) # ln(v)

# Single ionisation above threshold
print("pI1[1s]  =", Ne.pI1(pts, '1s'))

# Double ionisation, 1e above threshold
print("pI21[1s2p] =", Ne.pI21(pts, '1s2p'))

# Double ionisation, both electrons above threshold
print("pI2[1s2p]  =", Ne.pI2(pts, '1s2p'))


#######################################################
# Semi-inclusive differential ionisation probabilities
#######################################################

print('\nSemi-inclusive differential probabilities:')

Ne.load_probabilities(inclusive=True, e_threshold=1.0)

# Array of points to evaluate
pts1 = np.array([[log(1.3), log(0.005)]]) # (ln(En), ln(v))

print("dpI1[1s] =", Ne.dpI1(pts1, '1s'))

# Note that this function returns dP/dE, not dP/dlnE

#######################################################
# Semi-inclusive integrated ionisation probabilities
#######################################################

print('\nSemi-inclusive integrated probabilities (Eth = 1.0 keV):')

Ne.load_probabilities(inclusive=True, integrated=True, e_threshold=1.0)

# Array of velocities to evaluate
pts = np.array([log(0.005)]) # ln(v)

print("pI1[1s] =", Ne.pI1(pts, '1s'))

#######################################################