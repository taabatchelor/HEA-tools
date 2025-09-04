'''plot DFT vs predicted OH adsorption energies'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

# location and name of csv files with training and test adsorption energies
csvTrainDFT = '../DFT_histogram/OH_train.csv'
csvTestDFT = '../DFT_histogram/OH_test.csv'
csvTrainPred = 'OH_train_pred.csv'
csvTestPred = 'OH_test_pred.csv'

# load training adsorption energies
ETrainDFT = np.loadtxt(csvTrainDFT, usecols=-1, delimiter=',')
ETrainPred = np.loadtxt(csvTrainPred, usecols=0, delimiter=',')

# load test adsorption energies
ETestDFT = np.loadtxt(csvTestDFT, usecols=-1, delimiter=',')
ETestPred = np.loadtxt(csvTestPred, usecols=0, delimiter=',')

# initiate figure
fig, ax = plt.subplots(figsize=(4, 3))

# start and end of adsorption energies to plot
start, stop = 0.0, 1.5

# plot fontsize
fs = 12

# plot solid diagonal line
ax.plot([start,stop],[start,stop], 'k-', linewidth=1.0,
		label=r'$\Delta E_{\mathrm{pred}} = \Delta E_{\mathrm{DFT}}$')

# plot dashed diagonal lines 0.1 eV above and below solid diagonal line
pm = 0.1
ax.plot([start,stop],[start+pm, stop+pm], 'k--', linewidth=1.0,
label=r'$\pm %.1f \mathrm{eV}$'%pm)
ax.plot([start+pm,stop],[start,stop-pm], 'k--', linewidth=1.0)

## training set ##
# Root mean square deviation (RMSD) between DFT and predicted energies
RMSD = np.sqrt( sum( (ETrainPred - ETrainDFT)**2 ) / len(ETrainPred) )

# show RMSD as text in the plot
txtTrain = 'RMSD (training set)'
ax.text(0.005,0.93,'%s = %.3f eV'%(txtTrain, RMSD), transform=ax.transAxes, fontsize=10)

## test set ##
# Root mean square deviation (RMSD) between DFT and predicted energies
RMSD = np.sqrt( sum( (ETestPred - ETestDFT)**2 ) / len(ETestPred) )

# show RMSD as text in the plot
txtTest = 'RMSD (test set)'
ax.text(0.005,0.93-0.075,'%-21s = %.3f eV'%(txtTest, RMSD), transform=ax.transAxes, fontsize=10)

# number of training and test samples
nTrain = len(ETrainPred)
nTest = len(ETestPred)

# plot training DFT vs predicted adsorption energies
ax.plot(ETrainDFT, ETrainPred, 'b.', markersize=2, label='training set (%d)'%nTrain)

# plot test DFT vs predicted adsorption energies
ax.plot(ETestDFT, ETestPred, 'rx', markersize=3, label='test set (%d)'%nTest)

# set legend style
plt.rc('legend', **{'fontsize':fs})
ax.legend(loc='right', bbox_to_anchor=(1.04, 0.2), fontsize=10,
		  markerscale=1.5, frameon=False, labelspacing=0.1,
		  handlelength=1.5)


# set style of labels
plt.xlabel(r'$\Delta E_{\mathrm{DFT}} \, [\mathrm{eV}]$', fontsize=fs)
plt.ylabel(r'$\Delta E_{\mathrm{pred}} \, [\mathrm{eV}]$', fontsize=fs)
ax.set_xlim((start, stop))
ax.set_ylim((start, stop))

# specify style of ticks
majorTick, minorTick, tickFormat = 0.2, 0.05, '%.1f'

majorLocator = MultipleLocator(majorTick)
majorFormatter = FormatStrFormatter(tickFormat)
minorLocator = MultipleLocator(minorTick)

ax.xaxis.set_major_locator(majorLocator)
ax.xaxis.set_major_formatter(majorFormatter)
ax.xaxis.set_minor_locator(minorLocator)

ax.yaxis.set_major_locator(majorLocator)
ax.yaxis.set_major_formatter(majorFormatter)
ax.yaxis.set_minor_locator(minorLocator)

for xTick,yTick in zip(ax.xaxis.get_major_ticks(), ax.yaxis.get_major_ticks()):
	xTick.label.set_fontsize(fs)
	yTick.label.set_fontsize(fs) 
ax.xaxis.label.set_fontsize(fs)
ax.yaxis.label.set_fontsize(fs)

ax.tick_params(which='both', top=True, right=True)
		
plt.tight_layout()

# save image
plt.savefig('OH_DFT_vs_pred.png', dpi=300)