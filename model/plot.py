import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, FuncFormatter

def pow4(x, pos):
    'The two args are the value and tick position'
    return '%d' % (x*1e-4)

def pow9(x, pos):
    'The two args are the value and tick position'
    return '%d' % (x*1e-9)

def applyPlotStyle(ax,													# axis object
				   xlim = [-0.5, 2.5], 									# [start,stop]
				   labels = [r'$\Delta E_{DFT}$ [eV]', r'frequency'], 	# [x,y]
				   xticks = [0.5, 0.1],									# [major,minor]
				   yticks = [10, 2],									# [major,minor]
				   tickformat = ['%.1f','%d'], 							# [x,y]
				   fontsize = 12,
				   ypow9 = False,
				   ypow4 = False):
				   
	plt.xlabel(labels[0])
	plt.ylabel(labels[1])
	
	majorXLocator = MultipleLocator(xticks[0]) # MAJOR X TICKS
	minorXLocator = MultipleLocator(xticks[1]) # MINOR X TICKS
	majorXFormatter = FormatStrFormatter(tickformat[0])

	majorYLocator = MultipleLocator(yticks[0]) # MAJOR Y TICKS
	minorYLocator = MultipleLocator(yticks[1]) # MINOR Y TICKS
	majorYFormatter = FormatStrFormatter(tickformat[1])	

	ax.xaxis.set_major_locator(majorXLocator)
	ax.xaxis.set_minor_locator(minorXLocator)
	ax.xaxis.set_major_formatter(majorXFormatter)

	ax.yaxis.set_major_locator(majorYLocator)
	ax.yaxis.set_minor_locator(minorYLocator)
	
	if ypow9:
		ax.yaxis.set_major_formatter(FuncFormatter(pow9))
	elif ypow4:
		ax.yaxis.set_major_formatter(FuncFormatter(pow4))
	else:
		ax.yaxis.set_major_formatter(majorYFormatter)

	for xTick in ax.xaxis.get_major_ticks():
		xTick.label.set_fontsize(fontsize)
	for yTick in ax.yaxis.get_major_ticks():
		yTick.label.set_fontsize(fontsize)
	ax.xaxis.label.set_fontsize(fontsize)
	ax.yaxis.label.set_fontsize(fontsize)

	ax.tick_params(which='both',top=True,right=True)
	ax.set_xlim((xlim[0],xlim[1]))
	
	return None