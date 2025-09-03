'''write 35 csv files with all possible fingerprints with the number of each element, the O adsorption
energy and the fingerprint multiplicity. This takes around 2.5 hours for this 5-zone surface description of O adsorption. '''

import sys
import numpy as np
import itertools as it

# load LinearRegression class
sys.path.append('../model')
from LinearRegression import LinearRegression
model = LinearRegression()

# location and name of training set
csvTrain = '../DFT_histogram/O_train.csv'

# load training adsorption energies
XTrain = np.loadtxt(csvTrain, dtype='int', delimiter=',', usecols=range(55))
ETrain = np.loadtxt(csvTrain, dtype='float', usecols=-1, delimiter=',')

# calculate linear regression parameters
w = model.weights(XTrain, ETrain, nMetals=5, zoneSizes=(3, 6, 3, 3, 3))

# predict energies of all possible fingerprints using linear regression parameters and save to file
model.all_fingerprints('O_all_slabs.csv', w, nMetals=5, zoneSizes=(3, 6, 3, 3, 3))