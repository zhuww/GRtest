import numpy as np
import matplotlib.pyplot as plt
from pylab import rand
from datatools.tempo import *
from astropy import coordinates as coord
from tools.Coordinate import *
from psrpbdot import M1
import numpy.random as npr
from cmath import *

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
import numpy.linalg as linalg
from GalacticGeometry import *

G = 6.673e-8
c = 2.99792458e10
PI = np.pi
AU = 1.469e13
Tsun = 4.925490947e-6
secperday = 24*3600
secperyear = secperday*365.24218967
kpc = 3.0857e21

"""load the parfile"""
#pf = PARfile('mcmcresult.par')
#pf = PARfile('1713.final.par')
pf = PARfile('Feb.T1.RN.par')
Eerr = float(pf.E[1])
OMerr = float(pf.OM[1])/180.*np.pi
"""calculate coordinate transfermation matrix using the parfile"""
T, GT = GetTransferMatrix(pf)

def EccArea(ECC, EF, THETA):
    """ calculate the "area" in the four 3-sigma points for given ECC_observed, ECC_forced and the angle between them.
    """
    global Eerr, OMerr
    #THETA[THETA<0]+=(np.pi*2)
    Emax, Emin = ECC +Eerr*3 , ECC-Eerr*3
    Temax, Temin = THETA +OMerr*3, THETA-OMerr*3
    fourpnts = [[(x*np.cos(y)-EF[i]+x*np.sin(y)*1j) for x in (Emax[i], Emin[i]) for y in (Temax[i], Temin[i])] for i in range(len(ECC))]
    areas = []
    for v4 in fourpnts:
        abs4 = [abs(v) for v in v4]
        p4 = np.array([phase(v) for v in v4])
        if p4.min() < 0. and p4.max() - p4.min() > np.pi:
            p4[p4<0.] += 2.*np.pi
        if np.mean(abs4) < Eerr*3 and abs(np.mean(p4))<OMerr*3:
            #print "hei, it's happening ******** "
            value = max((np.log10(min(max(abs4), 0.05)) - np.log10(max(min(abs4),1.e-6))),0)*np.pi*2
        else:
            value = max((np.log10(min(max(abs4), 0.05)) - np.log10(max(min(abs4),1.e-6))),0)*(max(p4)-min(p4))
            if max(p4)-min(p4)>1.:
                print 'what?', value
        areas.append(value)
    areas = np.array(areas)
    return areas

def Pintegrant(M1, M2, PB, F0, ECC, PMRA, PMDEC, PX, SINI, PAASCNODE, OM, w):
    """ Calculate the "probability" for given Delta and timing parameters
    """
    Mtot = M1+M2

    wserr = 0.9e5 #Kogut et al. 1993, Fixsen et al 1996, Hinshaw et al. 2009
    wsolar = 369.e5 + npr.randn()*wserr#See ref below (aaa+13 Planck Team: Aghanim, N. et al. 2013. Planck confirms this using a different method)
    lws, bws = 263.99/180.*np.pi, 48.26/180.*np.pi
    w_s = wsolar * (GT.I * np.matrix((np.cos(bws)*np.cos(lws),np.cos(bws)*np.sin(lws),np.sin(bws))).T)
    ws_NSEW = T * w_s
    D = kpc/PX
    
    wx = w * (np.matrix((-1., 0., 0.)).T)
    wy = PMRA*1.e-3/60./60.*np.pi/180./secperyear * D * (np.matrix((0.,-1.,0.)).T)
    wz = PMDEC*1.e-3/60./60.*np.pi/180./secperyear * D * (np.matrix((0.,0.,1.)).T)

    w_ns = (wx + wy + wz)
    w =  w_ns + ws_NSEW

    incang = np.arcsin(SINI)
    Omgang = PAASCNODE/180.*np.pi
    A = -1./np.tan(incang)
    B = -1./np.tan(Omgang)
    C = 1.
    A_ref = np.matrix((0, -1.* np.sin(Omgang), np.cos(Omgang)))

    n_orb = np.matrix((A, B, C)).T
    n_orb= n_orb/linalg.norm(n_orb)

    w_proj = w - n_orb * (w.T*n_orb) 
    w_leg = linalg.norm(w_proj)
    w_dir = w_proj/w_leg
    w_ang = np.arccos(A_ref * w_dir)

    #EF = lambda w:0.21* M1 * w * (PB**2) * c**2 *F0/24 / PI /G /Mtot/Msun
    EF = lambda w:0.21* M1 * w * (PB**2) *F0/24 / PI /Mtot /c /Tsun

    theta, ef = float(w_ang - OM/180.*np.pi), EF(w_leg) 
    return theta, ef

''' load in the MCMC results for Delta estimation '''
import cPickle as pickle
import sys
#from Coordinate import RA, Dec
secperday = 24*3600

#dic = pickle.load(open('bestpar.p', 'rb'))
#best = dic['BEST']
#plist = dic['parameters'] 
#MChain = pickle.load(open('TinyMChain.p','rb'))
plist = open('pars.txt', 'r').readlines()
pl = [p.strip('\n') for p in plist]
plist = pl
#MChain = pickle.load(open('SmallMChain.p','rb'))
#MChain = pickle.load(open('MChain.p','rb'))
#MarkovChain = MChain['Chain']
#MCMCSize = len(MarkovChain)
MarkovChain = np.loadtxt('chain_1.0.txt')
MCMCSize, cols = MarkovChain.shape
pi = 3.141592653589793
G = 6.673e-8
c = 2.99792458e10
Msun = 1.9882e33
twopi = 6.283185307179586
fac = 1.536e-16 
Tsun = 4.925490947e-6
#if0 = plist.index('F0')
im2 = plist.index('M2')
ipb = plist.index('PB')
#isini = plist.index('SINI')
ii = plist.index('KIN')
ia = plist.index('A1')
iecc = plist.index('ECC')
ipx = plist.index('PX')
#iomega = plist.index('PAASCNODE')
iomega = plist.index('KOM')
iom = plist.index('OM')
ipmra = plist.index('PMRA')
ipmdec = plist.index('PMDEC')
#F0 = np.array([float(p[if0]) for p in MarkovChain])
M2 = MarkovChain[:,im2]
PB = secperday * MarkovChain[:,ipb]
SINI = np.sin(MarkovChain[:,ii]/180.*np.pi)
OM = MarkovChain[:,iom]
a =  MarkovChain[:, ia]
#M1 = (PB/2/pi*np.sqrt(G*(M2*SINI)**3/a**3)-M2)#/Msun
M1 = PB/2/pi*(np.sqrt(Tsun*(M2*SINI)**3/a**3))-M2
M2 = M2#/Msun
ECC = MarkovChain[:,iecc]
PX = MarkovChain[:,ipx]
PAASCNODE = MarkovChain[:,iomega]
PMRA = MarkovChain[:,ipmra]
PMDEC = MarkovChain[:,ipmdec]
#if0 = plist.index('F0')
#im2 = plist.index('M2')
#ipb = plist.index('PB')
#isini = plist.index('SINI')
#ia = plist.index('A1')
#ichisq = plist.index('chisq')
#iecc = plist.index('E')
#ipx = plist.index('PX')
#iomega = plist.index('PAASCNODE')
#iom = plist.index('OM')
#M2 = np.array([float(p[im2]) for p in MarkovChain])
#PB = np.array([float(p[ipb])*secperday for p in MarkovChain])
#SINI = np.array([float(p[isini]) for p in MarkovChain])
#A = np.array([float(p[ia]) for p in MarkovChain])
##print PB, M2, SINI, 
#F0 = np.array([float(p[if0]) for p in MarkovChain])
#M1 = (PB/2/pi*np.sqrt(Tsun*(M2*SINI)**3/A**3)-M2)
#ECC = np.array([float(p[iecc]) for p in MarkovChain])
#PX = np.array([float(p[ipx]) for p in MarkovChain])
#PMRA = np.array([float(p[ipmra]) for p in MarkovChain])
#PMDEC = np.array([float(p[ipmdec]) for p in MarkovChain])
#PAASCNODE = np.array([float(p[iomega]) for p in MarkovChain])
#OM = np.array([float(p[iom]) for p in MarkovChain])
#chisq = [p[ichisq] for p in MarkovChain]
#bestidx = chisq.index(min(chisq))
Dmean = 1/PX.mean()
wy_m = PMRA.mean()*1.e-3/60./60.*np.pi/180./secperyear * Dmean 
wz_m = PMDEC.mean()*1.e-3/60./60.*np.pi/180./secperyear * Dmean 
w_m = np.sqrt(wy_m**2 + wz_m**2) * kpc
#print 'mean proper motion', w_m * kpc
w = np.array([ y for y in npr.normal(0., w_m, MCMCSize)])
F0 = float(pf.F0[0])

#print SINI.mean()
#print M1.mean()
#print M2.mean()
#sys.exit(0)

def Integ(alpha):
    global M1, M2, PB, F0, ECC, PMRA, PMDEC, PX, SINI, PAASCNODE, OM, w
    #pres = np.array([Pintegrant(M1[i], M2[i], PB[i], F0[i], ECC[i], PMRA[i], PMDEC[i], PX[i], SINI[i], PAASCNODE[i], OM[i], w[i]) for i in range(len(w))])
    pres = np.array([Pintegrant(M1[i], M2[i], PB[i], F0, ECC[i], PMRA[i], PMDEC[i], PX[i], SINI[i], PAASCNODE[i], OM[i], w[i]) for i in range(len(w))])
    THETA = 1.* pres[...,0]
    EF = alpha * pres[...,1]
    Area1 = sum(EccArea(ECC, EF, THETA))
    Area2 = sum(EccArea(ECC, EF, np.pi-THETA))
    Areas = Area1 + Area2
    return Areas


alpha = np.arange(5.e-22, 1.e-19, 5.e-22) #ingrid setting
res = np.array([Integ(a) for a in alpha])

np.save(open('alpha3_jae.npy', 'wb'), res)

cdf = [res[0]]
for r in res[1:]:
    s = cdf[-1] + r
    cdf.append(s)
sumres = cdf[-1]
cdf = np.array(cdf)/sumres

for i,c in enumerate(cdf):
    if c > 0.95:
        print alpha[i], c
        break

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(alpha, res, '-')
#ax.hist(delta, bins=50, normed=1)
ax.semilogx(nonposy='clip')
#ax.semilogy(nonposy='clip')
plt.show()
