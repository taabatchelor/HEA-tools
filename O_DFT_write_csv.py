'''save a csv file containing DFT O adsorption energies and their fingerprints for both
the test and training set'''

import sys
from ase.db import connect
import numpy as np
import itertools as it

# update system path to be able to import the Slab class
sys.path.append('../model')

# import Slab class
from Slab import Slab

# location and name of database
dbTrain = connect('../DFT_calculations/train.db')
dbTest = connect('../DFT_calculations/test.db')

# molecular DFT reference energies (eV)
refs = {'H2O': -12.261853, 'H2': -6.673272}

# metals in alloy
metals = ('Ir', 'Pd', 'Pt', 'Rh', 'Ru')

# loop through training set, then test set
for i, db in enumerate([dbTrain, dbTest]):
	
	# name output csv-files appropriately 
	if i == 0:
		filename = 'O_train.csv'
	elif i == 1:
		filename = 'O_test.csv'
	
	# initiate adsorption energy and fingerprint lists
	energies, fingerprints = [], []	

	# loop through relaxed slabs
	for row in db.select(relaxed=1, type='O'):
		
		# load atoms object of database entry
		atoms = db.get_atoms(row.id)

		# assign atoms to slab object
		slab = Slab(atoms)

		# omit if slab height has increased more than 10%
		if slab.distorted(threshold=1.10):
			print('distorded slab id: %d'%row.slabId)
			continue

		# omit if O is not in an fcc hollow site
		slab.get_site()
		if slab.site is not 'fcc':
			print('not fcc slab id: %d'%row.slabId)
			continue

		# get DFT energy of slab with O
		ESlabO = row.energy

		# get DFT energy of corresponding slab without O
		ESlab = db.get(relaxed=1, type='slab', slabId=row.slabId).energy
		
		# calculate adsorption energy relative to H2O and H2
		## H2O + * -> O* + H2
		energies.append(ESlabO + refs['H2'] - ESlab - refs['H2O'])
		
		# get the fingerprint of the on-top adsorption site 
		fingerprints.append(slab.features(metals, onTop=False, 
							zones=('ens', 'sf', 'ssf', 'sn', 'ssn')))

	# save as numpy arrays
	energies = np.array([energies])
	fingerprints = np.array(fingerprints)

	# number of features in fingerprint
	nFeat = fingerprints.shape[1]

	# save fingerprints and adsorption energies to a csv file
	output = np.concatenate((fingerprints, energies.T), axis=1)
	np.savetxt(filename, output, fmt=['%d']*nFeat + ['%.5f'], delimiter=',')