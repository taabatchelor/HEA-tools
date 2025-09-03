'''generates the histogram for the distribution of O adsorption energies. Takes around 20 min.'''

import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
from math import factorial as fact
import sys

# update system path to be able to import the plot function
sys.path.append('../model')

# import plot function
from plot import applyPlotStyle

# number of metals
nMetals = 5

# number of atoms in zones
zoneSizes = (3, 6, 3, 3, 3)

# location and name of csv file (without numbering)
filename = 'O_all_slabs.csv'

nEns = fact(nMetals+zoneSizes[0]-1)/(fact(zoneSizes[0])*fact(nMetals-1))
start, stop, spacing = -0.5, 2.5, 0.0075
binedges = np.arange(start, stop+spacing, spacing)
centerbins = binedges[:-1]+spacing/2
counts = np.zeros(len(binedges)-1, dtype='int64')

# define figure
fig, ax = plt.subplots(figsize=(4,3))


###### color the ensembles in the same way as the DFT histogram

# colors
colors = cm.rainbow(np.linspace(0, 1, nEns))

# load ensemble ids and DFT adsorption energies
ensembles = np.loadtxt('../DFT_histogram/O_train.csv', 
					   dtype='int', delimiter=',', usecols=range(nEns))
energies = np.loadtxt('../DFT_histogram/O_train.csv',
					  dtype='float', delimiter=',', usecols=-1)

# distribute energies according to their adsorption ensemble
energyEns = [[] for _ in range(nEns)]
for ensemble, energy in zip(ensembles, energies):
	([ensId],) = np.where(ensemble==1)
	energyEns[ensId].append(energy)
	
# means for each ensemble (from DFT energies)
means = [0.]*nEns
for i, ensEnergies in enumerate(energyEns):
	means[i] = sum(ensEnergies)/len(ensEnergies)

# ordering based on means
order = [idx for _, idx in sorted(zip(means, range(nEns)))]
orderId = [order.index(idx) for idx in range(nEns)]

##### end of coloring


# plot individual ensembles
for i in range(nEns):

	# label only the first instance
	if i == 0:
		label = 'individual sites'
	else:
		label = ''
		
	parts = filename.rpartition('.')
	fname = parts[0] + '_%d'%i + parts[1] + parts[2]
	energies = np.loadtxt(fname, dtype='float', delimiter=',', usecols=5)
	mults = np.loadtxt(fname, dtype='int64', delimiter=',', usecols=6)
	count, _ = np.histogram(energies, binedges, weights=mults)
	counts += count
	
	# plot individual ensemble
	plt.hist(centerbins, bins=binedges, weights=count, 
			 facecolor=colors[orderId[i]], ec='black', alpha=0.75, 
			 histtype='stepfilled', zorder=1, label=label)
			 
	print('file %d read'%i)

# plot total
plt.hist(centerbins, bins=binedges, weights=counts, 
		 facecolor='grey', ec='black', alpha=0.75,
		 histtype='stepfilled', zorder=0, label='total')

# axis ticks and labeling
applyPlotStyle(ax,
			   xlim = [start, stop],
			   labels = [r'$\Delta E_{pred}$ [eV]', r'frequency $\times 10^{9}$'],
			   xticks = [0.5, 0.1],
			   yticks = [5e9, 1e9],
			   tickformat = ['%.1f','%d'],
			   fontsize = 12,
			   ypow9 = True)

# legend style
lgd = ax.legend(loc='center', 
				handlelength=1.25,
				bbox_to_anchor=(0.5, 1.06),
				handletextpad=0.2,
				columnspacing=0.7,
          		ncol=2,
				fancybox=False,
				shadow=False,
				prop={'size': 14},
				frameon=False)

plt.tight_layout()

# save image to file
plt.savefig('O_pred_histogram.png',
			bbox_extra_artists=(lgd, ),
			bbox_inches='tight',
			dpi=300)