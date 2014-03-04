from pylab import *
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import pickle 
from matplotlib.ticker import NullFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.optimize as opt



MChain = pickle.load(open('MChain.p', 'r'))['Chain']
best = pickle.load(open('bestpar.p', 'r'))['BEST']

MChain = np.array(MChain)
x = MChain[:,0]
y = MChain[:,1]
x = np.array(x)
x *=1.e12
y = np.array(y)
y *=1.e4
#print x.shape, y.shape
#print x
#print y
bs = len(x)
H, xed, yed = np.histogram2d(x,y, bins=25)
#extent = [yed[-1], yed[0], xed[-1], xed[0]]
extent = [ xed[0], xed[-1], yed[-1], yed[0]]
#extent = [ yed[0], yed[-1], xed[0], xed[-1]]
#print H.shape, xed.shape, yed.shape
#print extent
#print H
xbin, ybin = [], []
for i in range(1, len(xed)):
    xbin.append((xed[i] + xed[i-1])/2.)
    ybin.append((yed[i] + yed[i-1])/2.)




from scipy import stats

#xkde = stats.gaussian_kde(x)
#ykde = stats.gaussian_kde(y)

from numpy import *

xm = mean(x)
ym = mean(y)
#xsigma = stats.tstd(x)
#ysigma = stats.tstd(y)
xsigma = np.std(x)
ysigma = np.std(y)
#xpdf = norm(xm, xsigma).pdf()
#ypdf = norm(ym, ysigma).pdf()
xpdf = stats.norm(xm, xsigma)
ypdf = stats.norm(ym, ysigma)

#import sys
#sys.exit(0)

nullfmt   = NullFormatter()         # no labels

# definitions for the axes 
#left, width = 0.1, 0.65
#bottom, height = 0.1, 0.65
#bottom_h = left_h = left+width+0.02

#rect_scatter = [left, bottom, width, height]
#rect_histx = [left, bottom_h, width, 0.2]
#rect_histy = [left_h, bottom, 0.2, height]

# start with a rectangular Figure
plt.figure(1, figsize=(8,8))
cnorm = cm.colors.Normalize(vmax=abs(H).max(), vmin=0.)
cmap = cm.binary

axScatter = plt.subplot(111)
#axScatter = plt.axes(rect_scatter)
#axHistx = plt.axes(rect_histx)
#axHisty = plt.axes(rect_histy)


# no labels
#axHistx.xaxis.set_major_formatter(nullfmt)
#axHisty.yaxis.set_major_formatter(nullfmt)

# the scatter plot:
#axScatter.scatter(x, y, edgecolor=None)
total = np.sum(H)

def onesig(lev):
    return np.sum(H[H>lev])/total - 0.682
def twosig(lev):
    return np.sum(H[H>lev])/total - 0.954
def trisig(lev):
    return np.sum(H[H>lev])/total - 0.997

#print onesig(50000), twosig(10000), trisig(10000)

#levels = opt.bisect(onesig, 1, 100000), opt.bisect(twosig, 1, 50000), opt.bisect(trisig, 1, 10000)
levels = opt.bisect(onesig, 1, 100000), opt.bisect(twosig, 1, 50000)#, opt.bisect(trisig, 1, 10000)
levellabel = {levels[0]:'68%', levels[1]:'95%'}
#sys.exit(0)

#print levels
def fomt(x):
    #return '%s$\sigma$' % (i[0]/2)
    #lbs = ['99.7%%', '95%%', '68%%' ] * 3
    return levellabel[x]

#levels = [500, 1000., 1500]
#axScatter.imshow(H, extent=extent, aspect='auto', interpolation='bicubic', cmap=cm.get_cmap(cmap, len(levels)-1))
#axScatter.imshow(H)
CS = axScatter.contour(xbin,ybin, H, levels=levels)
clabel(CS, inline=1, fmt=fomt, manual=True)
#print H.min(), H.max()
#axScatter.set_aspect(1.)
xlabel(r'$\dot{G}/G [10^{-12}{\rm yr}^{-1}]$')
ylabel('$\kappa_D[10^{-4}]$')

axScatter.axvline(0, linestyle='--', color='k')
axScatter.axhline(0, linestyle='--', color='k')
#axScatter.axvspan((-0.7-7.6)/10, (-0.7+7.6)/10, facecolor='0.5', alpha=0.3, linewidth=2)
axScatter.axvspan((-0.7-7.6)/10, (-0.7+7.6)/10, fill=False, hatch='/', linewidth=2)

axScatter.plot([0.], [0.], 'ro', ms=10)
axScatter.text(0.,0.,'GR',color='r',fontsize=28)

axScatter.set_xlim([-2.,2.])
axScatter.set_ylim([-6,6])
#divider = make_axes_locatable(axScatter)
#axHistx = divider.append_axes("top", 1.2, pad=0.1, sharex = axScatter)
#axHisty = divider.append_axes("right", 1.2, pad=0.1, sharey = axScatter)

#plt.setp(axHistx.get_xticklabels() + axHisty.get_yticklabels(), visible=False)

# now determine nice limits by hand:
#Nbins = 20
#xmax = np.max(np.fabs(x))
#ymax = np.max(np.fabs(y))
#xbinwidth = 2*xmax/Nbins
#ybinwidth = 2*ymax/Nbins
#xlim = ( int(xmax/xbinwidth) + 1) * xbinwidth
#ylim = ( int(ymax/ybinwidth) + 1) * ybinwidth

#axScatter.set_xlim( (-xlim, xlim) )
#axScatter.set_ylim( (-ylim, ylim) )

#xbins = np.arange(-xlim, xlim + xbinwidth, xbinwidth)
#ybins = np.arange(-ylim, ylim + ybinwidth, ybinwidth)
#axHistx.hist(x, bins=xbins, normed=1)
#axHisty.hist(y, bins=ybins, orientation='horizontal', normed=1)
#axHistx.plot(xbins,xkde.evaluate(xbins),'r-')
#axHisty.plot(ykde.evaluate(ybins),ybins,'r-')
#axHistx.plot(xbins,xpdf.pdf(xbins),'g-')
#axHisty.plot(ypdf.pdf(ybins),ybins,'g-')

#axHistx.set_xlim( axScatter.get_xlim() )
#axHisty.set_ylim( axScatter.get_ylim() )

#for t1 in axHistx.get_xticklabels():
    #t1.set_visible(False)
#axHistx.set_yticks([0])
#for t1 in axHisty.get_yticklabels():
    #t1.set_visible(False)
#axHisty.set_xticks([0])

plt.draw()
plt.show()

from round import shortform as SF
from round import TexStyle as TS
print TS((xm,2*xsigma)), TS((ym, 2*ysigma))
secperday = 24*3600
secperyear = secperday * 365.24218967
H0 = 67.80 #+/- 0.77 km/s/Mpc Planck
Mpc_in_km  = 3.08568e19
T_Hubble = Mpc_in_km/H0/secperyear  #year
print T_Hubble
print 'factor of change  in Hubble time:', T_Hubble * (xm+3*xsigma)*1.e-12

