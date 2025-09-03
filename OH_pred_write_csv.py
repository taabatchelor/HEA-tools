'''write csv files with predicted OH adsorption energies in the same order as the energies used for training'''

import sys
import numpy as np

# load LinearRegression class
sys.path.append('../model')
from LinearRegression import LinearRegression
model = LinearRegression()

# load training set fingerprints and energies
csvTrain = '../DFT_histogram/OH_train.csv'
XTrain = np.loadtxt(csvTrain, dtype='int', delimiter=',', usecols=range(15))
ETrain = np.loadtxt(csvTrain, dtype='float', delimiter=',', usecols=-1)

# calculate linear regression parameters
w = model.weights(XTrain, ETrain, nMetals=5, zoneSizes=(1, 6, 3))

# load test set fingerprints and energies
csvTest = '../DFT_histogram/OH_test.csv'
XTest = np.loadtxt(csvTest, dtype='int', delimiter=',', usecols=range(15))
ETest = np.loadtxt(csvTrain, dtype='float', delimiter=',', usecols=-1)

# predict energies from fingerprints using linear regression parameters
EPredTrain = model.predicted_energies(XTrain, w)
EPredTest = model.predicted_energies(XTest, w)

# write to csv files
np.savetxt('OH_train_pred.csv', EPredTrain.T, fmt='%.5f', delimiter=',')
np.savetxt('OH_test_pred.csv', EPredTest.T, fmt='%.5f', delimiter=',')