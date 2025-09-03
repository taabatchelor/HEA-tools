import matplotlib.pyplot as plt
import numpy as np
import itertools as it
from math import factorial as fact
import sys

# update system path to be able to import the plot function
sys.path.append('../model')

# import plot function
from plot import applyPlotStyle

# number of metals
nMetals = 5

# number of atoms in zones
zoneSizes = (1, 6, 3)

# location and name of csv file
filename = 'OH_all_slabs.csv'

# load data from file
energies = np.loadtxt(filename, dtype='float', delimiter=',', usecols=5)
mults    = np.loadtxt(filename, dtype='int32', delimiter=',', usecols=6)

def prod(ns):
	'return the product of the numbers given'
	'ns		list of floats		numbers to multiply'
	p = 1
	for n in ns:
		p *= n
	return p

def count(config, nMetals=5):
	'return the count list of the configuration given'
	'config		tuple of ints	metal types as integers {0...nMetals-1}'
	'nMetals	int				number of metals'
	
	counts = [0]*nMetals
	for c in config:
		for m in range(nMetals):
			if c == m:
				counts[m] += 1
				
	return tuple(counts)

def no_of_zone_ensembles(zoneSize, nMetals):
	'return number of ensembles for the zone with the given number of atoms'
	return fact(nMetals+zoneSize-1)/(fact(zoneSize)*fact(nMetals-1))	

def ensemble_from_index(idx, nMetals=5, zoneSizes=(1, 6, 3)):
	'return the adsorption ensemble id corresponding to the index idx'
	'idx		int		fingerprint index  e.g. {1...36750}'
	'nMetals	int		number of metals'
	'zoneSizes	tuple	number of atoms in each zone in the fingerprint'	
	
	# number of zones
	nZones = len(zoneSizes)
	
	# number of ensembles for each zone
	nEns = np.zeros(nZones)
	for i, zoneSize in enumerate(zoneSizes):
		nEns[i] = no_of_zone_ensembles(zoneSize, nMetals)
	
	# ids of each zone
	zoneIds = [0]*nZones
	
	# last two zone ids
	zoneIds[-1] = int(idx % nEns[-1])
	zoneIds[-2] = int((idx-zoneIds[-1])/nEns[-1] % nEns[-2])
	
	# the remaining zone ids
	for i in range(nZones-3, -1, -1): # e.g. nZones=5 -> [2, 1, 0], nZones=3 -> [0]
		sumx = sum(zoneIds[j]*prod(nEns[k] for k in range(j+1, nZones)) for j in range(i+1, nZones-1))
		denominator = prod(nEns[j] for j in range(i+1, nZones))
		zoneIds[i] = int((idx-zoneIds[-1]-sumx)/denominator % nEns[i])

	return zoneIds[0]

nEns = no_of_zone_ensembles(zoneSizes[0], nMetals)

# sort according to ensemble ids

# number of data
nData = len(mults)

# data for each adsorption ensemble
nEach = nData/nEns

# make arrays with column for each adsorption ensemble
energiesArray = np.zeros((nEach, nEns), dtype='float')
multsArray = np.zeros((nEach, nEns), dtype='int')

# define counting array
counters = [0]*nEns

for i in range(nData):
	# load fingerprint from index
	ensId = ensemble_from_index(i, nMetals, zoneSizes)
	
	# put the data in appropriate columns and rows in the arrays
	energiesArray[counters[ensId], ensId] = energies[i]
	multsArray[counters[ensId], ensId] = mults[i]
	counters[ensId] += 1

# define figure
fig, ax = plt.subplots(figsize=(4,3))

# specify bin sizes
start, stop, spacing = 0.0, 1.5, 0.007
bins = np.arange(start, stop, spacing)

# define colors of individual metals
metals  = ('Ir', 'Pd', 'Pt', 'Rh', 'Ru')
colors  = ('green', 'blue', 'orange', 'red', 'magenta')
zorders = (4, 1, 2, 3, 5)

# plot individual ensembles
for i in range(nEns):
	# plot enegies in column i in bins multiplied with the fingerprint multiplicities
	plt.hist(energiesArray[:, i], bins, weights=multsArray[:, i],
			 facecolor=colors[i], ec='black', alpha=0.75,
			 histtype='stepfilled', zorder=zorders[i], label=metals[i])

# plot total
plt.hist(energies, bins, weights=mults,
		 facecolor='grey', ec='black', alpha=0.75,
		 histtype='stepfilled', zorder=0, label='total')

# number of samples
nSamples = sum(mults)

# axis ticks and labeling
applyPlotStyle(ax,
			   xlim = [start, stop],
			   labels = [r'$\Delta E_{pred}$ [eV]', r'frequency $\times 10^{4}$'],
			   xticks = [0.5, 0.1],
			   yticks = [5e4, 1e4],
			   tickformat = ['%.1f','%d'],
			   fontsize = 12,
			   ypow4 = True)

# legend style
h, l = ax.get_legend_handles_labels()
handles, labels = (h[4], h[0], h[3], h[2], h[1], h[5]), (l[4], l[0], l[3], l[2], l[1], l[5])
lgd = ax.legend(handles, labels,
				loc='center',
				handlelength=1.25,
				bbox_to_anchor=(0.5, 1.06),
				handletextpad=0.2,
				columnspacing=0.7,
				ncol=6,
				fancybox=False,
				shadow=False,
				prop={'size': 14},
				frameon=False)
plt.tight_layout()

# save image to file
plt.savefig('OH_pred_histogram.png',
			bbox_extra_artists=(lgd, ),
			bbox_inches='tight',
			dpi=300)