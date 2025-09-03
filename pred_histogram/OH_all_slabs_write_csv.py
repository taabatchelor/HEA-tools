'''write csv file with all possible fingerprints with the number of each element, the OH adsorption
energy and the fingerprint multiplicity'''

import sys
import numpy as np
import itertools as it

# load LinearRegression class
sys.path.append('../model')
from LinearRegression import LinearRegression
model = LinearRegression()

# location and name of training set
csvTrain = '../DFT_histogram/OH_train.csv'

# load training adsorption energies
XTrain = np.loadtxt(csvTrain, dtype='int', delimiter=',', usecols=range(15))
ETrain = np.loadtxt(csvTrain, dtype='float', usecols=-1, delimiter=',')

# calculate linear regression parameters
w = model.weights(XTrain, ETrain, nMetals=5, zoneSizes=(1, 6, 3))

# predict energies of all possible fingerprints using linear regression parameters and save to file
model.all_fingerprints('OH_all_slabs.csv', w, nMetals=5, zoneSizes=(1, 6, 3))