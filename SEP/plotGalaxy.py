import numpy as np
import matplotlib.pyplot as plt
from pylab import rand
import pulsars

data = np.genfromtxt(open('log_arms.out', 'r'), dtype = [('Num', 'i'), ('n', 'i'), ('x', 'f8'),('y', 'f8')])[1:]
UniqNum = np.unique(data['Num'])
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, aspect='equal')
solardist = 8.3

for num in UniqNum:
    x = data['x'][data['Num'] == num]
    y = data['y'][data['Num'] == num]
    ax.plot(x,y,'b-')
    ax.plot([0],[solardist], 'yo', ms=15)
    ax.plot([0],[0], 'k*', ms=15)
    #ax.text(0,0,'GC')
    ax.set_xlabel('x (kpc)')
    ax.set_ylabel('y (kpc)')
    ax.set_xlim((-4.0,6.8))
    ax.set_ylim((-0.2,10.6))

#psrlist=["0407","0437","1045","1125","1216","1455","1600","1640","1643","1709","1711","1713","1732","1745","1751","1804","1810","1855","1918","1933","2016","2019","2033","2129","2229"]
#allpsrs = pulsars.ATNFpulsarlist()
allpsrs = pulsars.allpulsars()
print allpsrs
psrlist = ['J0030+0451', 'J0711-6830', 'J1024-0719', 'J1730-2304', 'J1744-1134', 'J1905+0400', 'B1937+21', 'J1944+0907', 'J2124-3358', 'J2322+2057'] + ['J0437-4715','J0613-0200','J0751+1807', 'J1012+5307', 'J1023+0038', 'J1045-4509', 'J1455-3330', 'J1600-3053', 'J1640+2224', 'J1643-1224', 'J1709+2313', 'J1713+0747', 'J1853+1303', 'B1855+09', 'J1903+0327', 'J1909-3744', 'J1910+1256', 'J1911-1114', 'J1918-0642', 'B1953+29', 'B1957+20', 'J2019+2425', 'J2033+1734', 'J2051-0827', 'J2129-5721', 'J2229+2643', 'J2317+1439']
psrs = [allpsrs[psr] for psr  in psrlist]
for i,n in enumerate(psrs):
    if n.props['ecc'] == None:
        print psrlist[i], n.props
gl = np.array([float(psr.props['gl'])/180.*np.pi for psr in psrs])
gb = np.array([float(psr.props['gb'])/180.*np.pi for psr in psrs])
D = np.array([float(psr.props['dist']) for psr in psrs])
ecc = np.array([float(psr.props['ecc']) for psr in psrs])


x = D * np.cos(gb) * np.sin(gl)
y = solardist - D * np.cos(gb) * np.cos(gl) 
angle = np.arctan(y/x)*180./np.pi

from matplotlib.patches import Ellipse
#ells = [Ellipse(xy=[x[i], y[i]], width=(0.6+0.4*rand())*(1.-0.5), height=0.6*1, angle=90, fill=False, lw=2)  for i in range(len(psrs))]
#ells = [Ellipse(xy=[x[i], y[i]], width=(0.8+0.2*rand())*(1.-0.5), height=0.6*1, angle=rand()*360, fill=False, lw=2)  for i in range(len(psrs))]
#ells = [Ellipse(xy=[x[i], y[i]], width=0.6*1, height=0.6*(1.-0.5), angle=angle[i], fill=False, lw=2)  for i in range(len(psrs))]

#for e in ells:
    #ax.add_artist(e)
    #e.set_alpha(1.)
    #e.set_edgecolor('r')
#plt.plot(x,y, 'b.')
#plt.plot()



plt.show()

