import matplotlib.pyplot as plt
import numpy as np
from math import factorial
import sys

# update system path to be able to import the plot function
sys.path.append('../model')

# import plot function
from plot import applyPlotStyle

# metals in alloy
metals = ('Ir', 'Pd', 'Pt', 'Rh', 'Ru')

# number of metals
nMetals = len(metals)

# adsorption site ensemble size (e.g. 1 for on-top and 3 for (111) hollow sites)
ensSize = 1

# number of on-top adsorbing metals (ensembles)
nEns = factorial(nMetals+ensSize-1)//(factorial(nMetals-1)*factorial(ensSize))

# load ensemble ids and DFT adsorption energies
ensembles = np.loadtxt('OH_train.csv', dtype='int', delimiter=',', usecols=range(nEns))
energies = np.loadtxt('OH_train.csv', dtype='float', delimiter=',', usecols=-1)

# number of DFT samples
nSamples = len(energies)

# distribute energies according to their adsorption ensemble
energyEns = [[] for _ in range(nEns)]
for ensemble, energy in zip(ensembles, energies):
	([ensId],) = np.where(ensemble==1)
	energyEns[ensId].append(energy)

# define figure
fig, ax = plt.subplots(figsize=(4,3))

# start, end, and spacing of adsorption energies in histogram (eV)
start, stop, spacing = 0.0, 1.5, 0.05
bins = np.arange(start, stop, spacing)

# colors in the order given by "metals
colors = ('green', 'blue', 'orange', 'red', 'magenta')

# layering of histograms (higher numbers in front)
zorders = (4, 1, 2, 3, 5)

# plot individual ensembles
for i, es in enumerate(energyEns):
	plt.hist(es, bins,
			 facecolor=colors[i], ec='black', alpha=0.75,
			 histtype='stepfilled', zorder=zorders[i], label=metals[i])

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
plt.savefig('OH_DFT_histogram.png',
			bbox_extra_artists=(lgd, ),
			bbox_inches='tight',
			dpi=300)