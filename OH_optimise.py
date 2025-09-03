'''print the alloy composition that maximises the catalytic activity to the terminal'''

from scipy.optimize import minimize
import numpy as np
from random import random

# lower and upper bounds on molar fractions
lower, upper = 0.00, 1.00

# physical constants
e = 1.602176565e-19 # J/eV
kB = 1.3806488e-23 # J/K

# absolute temperature (K)
T = 300

# Boltzmann's constant times temperature (eV)
kBT = (kB/e)*T

# optimal adsorption energy (eV)
energyOpt = 1.10

# metals in alloy
metals = ('Ir', 'Pd', 'Pt', 'Rh', 'Ru')

# number of metals
nMetals = len(metals)

# location and path to csv with all possible fingerprints
csv = '../5_pred_histogram/OH_all_slabs.csv'

# load number of each metal, adsorption energies, and multiplicity of fingerprint
counts = np.loadtxt(csv, delimiter=',', dtype='uint8', usecols=range(nMetals))
energies = np.loadtxt(csv, delimiter=',', dtype='float', usecols=nMetals)
mults = np.loadtxt(csv, delimiter=',', dtype='uint16', usecols=nMetals+1)

def activity_neg(f):
	'''return the negative sum of individual fingerprint activity values
	f	list 	metallic molar fractions'''

	for fraction in f:
	
		# enforce bounds on f by returning large value
		if fraction < -0.01 or fraction > 1.01:
			return 1e2
	
	# initiate activity to zero	
	activity = 0
	
	# loop through all fingerprints
	for count, energy, mult in zip(counts, energies, mults):

		# calculate probability factor
		fProd = 1
		for f0, n0 in zip(f, count):
			fProd *= f0**n0
		
		# update the activity with the activity of the current fingerprint
		activity += fProd*mult*np.exp(-abs(energy - energyOpt) / kBT)
	
	# return the negative activity	
	return -activity

# lower and upper bounds on molar fractions
bnds = ((lower, upper),
		(lower, upper),
		(lower, upper),
		(lower, upper),
		(lower, upper))

# constrain the sum of the molar fractions to equal unity
cons = ({'type': 'eq',
		 'fun' : lambda f: sum(f)-1})

# print headers
print('local maxima:')
print('%8s %8s %8s %8s %8s %10s'%(metals[0], metals[1], metals[2], metals[3], metals[4], 'activity'))

# loop through 25 random initial compositions for the optimizer
for j in range(25):
	
	# generate a random initial composition
    lst = [random() for _ in range(nMetals)]
    sumlst = sum(lst)
    f = [elem/sumlst for elem in lst]
    
    # minimize the negative activity (i.e. maximize the activity)
    res = minimize(activity_neg, f, method='SLSQP', bounds=bnds, constraints=cons, tol=1e-7)
    
    # print result to terminal
    print('%8.3f %8.3f %8.3f %8.3f %8.3f %10.6f'
    	  %(res.x[0], res.x[1], res.x[2], res.x[3], res.x[4], -activity_neg(res.x)))