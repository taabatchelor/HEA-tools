import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from math import factorial
import sys

# update system path to be able to import the plot function
sys.path.append('../model')

# import plot function
from plot import applyPlotStyle

# number of metals
nMetals = 5

# adsorption site ensemble size (e.g. 1 for on-top and 3 for (111) hollow sites)
ensSize = 3

# number of on-top adsorbing metals (ensembles)
nEns = factorial(nMetals+ensSize-1)//(factorial(nMetals-1)*factorial(ensSize))

# load ensemble ids and DFT adsorption energies
ensembles = np.loadtxt('O_train.csv', dtype='int', delimiter=',', usecols=range(nEns))
energies = np.loadtxt('O_train.csv', dtype='float', delimiter=',', usecols=-1)

# number of DFT samples
nSamples = len(energies)

# distribute energies according to their adsorption ensemble
energyEns = [[] for _ in range(nEns)]
for ensemble, energy in zip(ensembles, energies):
	([ensId],) = np.where(ensemble==1)
	energyEns[ensId].append(energy)

# means for each ensemble
means = [0.]*nEns
for i, ensEnergies in enumerate(energyEns):
	means[i] = sum(ensEnergies)/len(ensEnergies)

# ordering based on means
order = [idx for _, idx in sorted(zip(means, range(nEns)))]
orderId = [order.index(idx) for idx in range(nEns)]

# define figure
fig, ax = plt.subplots(figsize=(4,3))

# start, end, and spacing of adsorption energies in histogram (eV)
start, stop, spacing = -0.5, 2.5, 0.075
bins = np.arange(start, stop, spacing)

# colors
colors = cm.rainbow(np.linspace(0, 1, nEns))

# plot individual ensembles
for i, es in enumerate(energyEns):
	
	# label only the first instance
	if i == 0:
		label = 'individual sites'
	else:
		label = ''
		
	plt.hist(es, bins,
			 facecolor=colors[orderId[i]], ec='black', alpha=0.75,
			 histtype='stepfilled', zorder=1, label=label)

# plot total
plt.hist(energies, bins,
		 facecolor='grey', ec='black', alpha=0.75,
		 histtype='stepfilled', zorder=0, label='total')

# axis ticks and labeling
applyPlotStyle(ax,
			   xlim = [start, stop],
			   labels = [r'$\Delta E_{DFT}$ [eV]', r'frequency'],
			   xticks = [0.5, 0.1],
			   yticks = [10, 5],
			   tickformat = ['%.1f','%d'],
			   fontsize = 12)

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
plt.savefig('O_DFT_histogram.png',
			bbox_extra_artists=(lgd, ),
			bbox_inches='tight',
			dpi=300)