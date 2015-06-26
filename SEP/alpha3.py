import numpy as np
import matplotlib.pyplot as plt
from pylab import rand
from datatools.tempo import *
from astropy import coordinates as coord
from astropy import constants as const
from tools.Coordinate import *
from psrpbdot import M1
import numpy.random as npr
from cmath import *
import os,sys

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
import numpy.linalg as linalg
from GalacticGeometry import *

PI = np.pi
Tsun = 4.925490947e-6
secperday = 24*3600
secperyear = secperday*365.24218967
c = const.c.cgs.value
kpc = const.kpc.cgs.value
AU = const.au.cgs.value #1.469e13
Msun = const.M_sun.cgs.value #1.9882e33
G = const.G.cgs.value

from optparse import OptionParser
usage = "usage: %prog [options] arg"
parser = OptionParser()
parser.add_option("-f", '--parfile', dest="parfile", help="par file")
(options, args) = parser.parse_args(args=sys.argv[1:])
print options
parfile = options.parfile
"""load the parfile"""
pf = PARfile(parfile)

"""calculate coordinate transfermation matrix using the parfile"""
T, GT = GetTransferMatrix(pf)#, PAASCNODE)

if pf.__dict__['BINARY'] in ['DD', 'T2']:
    E = float(pf.E[0])
    Eerr = float(pf.E[1])
    OMerr = float(pf.OM[1])/180.*np.pi
elif pf.__dict__['BINARY'] in ['ELL1']:
    E = np.sqrt(float(pf.EPS1[0]**2 + pf.EPS2[0]**2))
    Eerr = np.sqrt(float(pf.EPS1[1]**2 + pf.EPS2[1]**2))
    e1 = float(pf.EPS1[0])
    e2 = float(pf.EPS2[0])
    er1 = pf.EPS1[1]/pf.EPS1[0]
    er2 = pf.EPS2[1]/pf.EPS2[0]
    er = np.sqrt(float(er1**2 + er2**2))
    OM = np.arctan2(e1,e2)*180./np.pi % 360
    x = e1/e2
    OMerr = 1./(1. + x**2) * er * x * 180./np.pi

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
            #value = max((np.log10(min(max(abs4), 0.05)) - np.log10(max(min(abs4),1.e-6))),0)*np.pi*2
            value = max((np.log10(min(max(abs4), 0.05)) - np.log10(max(min(abs4),0.e-10))),0)*np.pi*2
        else:
            #value = max((np.log10(min(max(abs4), 0.05)) - np.log10(max(min(abs4),1.e-6))),0)*(max(p4)-min(p4))
            value = max((np.log10(min(max(abs4), 0.05)) - np.log10(max(min(abs4),0.e-10))),0)*(max(p4)-min(p4))
            if max(p4)-min(p4)>2.*np.pi:
                print 'what?', p4, value
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

dic = pickle.load(open('bestpar.p', 'rb'))
best = dic['BEST']
plist = dic['parameters'] 
#MChain = pickle.load(open('TinyMChain.p','rb'))
#MChain = pickle.load(open('SmallMChain.p','rb'))
MChain = pickle.load(open('MChain.p','rb'))
MarkovChain = MChain['Chain']
MCMCSize = len(MarkovChain)
#pi = 3.141592653589793
#G = 6.673e-8
#c = 2.99792458e10
#Msun = 1.9882e33
#twopi = 6.283185307179586
#fac = 1.536e-16 
#Tsun = 4.925490947e-6
ipb = plist.index('PB')
PB = np.array([float(p[ipb])*secperday for p in MarkovChain])
#isini = plist.index('SINI')
#SINI = np.array([float(p[isini]) for p in MarkovChain])

if pf.__dict__['BINARY'] in ['DD', 'T2']:
    iom = plist.index('OM')
    OM = np.array([float(p[iom]) for p in MarkovChain])
elif pf.__dict__['BINARY'] in ['ELL1']:
    ie1 = plist.index('EPS1')
    ie2 = plist.index('EPS2')
    E1 = np.array([float(p[ie1]) for p in MarkovChain])
    E2 = np.array([float(p[ie2]) for p in MarkovChain])
    OM = np.arctan2(E1,E2)*180./np.pi % 360

ia = plist.index('A1')
try:
    ipx = plist.index('PX')
    PX = np.array([float(p[ipx]) for p in MarkovChain])
except:
    PX = 1./np.array(float(pf.Dist[0]) + np.random.rand(MCMCSize)*float(pf.Dist[1]))
ipmra = plist.index('PMRA')
ipmdec = plist.index('PMDEC')
A = np.array([float(p[ia]) for p in MarkovChain])
#if0 = plist.index('F0')
#F0 = np.array([float(p[if0]) for p in MarkovChain])
#iecc = plist.index('E')
#ECC = np.array([float(p[iecc]) for p in MarkovChain])
PMRA = np.array([float(p[ipmra]) for p in MarkovChain])
PMDEC = np.array([float(p[ipmdec]) for p in MarkovChain])
#ia = plist.index('A1')
#a = np.array([float(p[ia]) for p in MarkovChain])
#ipx = plist.index('PX')
#PX = np.array([float(p[ipx]) for p in MarkovChain])

TS99 = lambda pb, p: (pb/p[1])**(1./p[0]) + p[2]
def PbToM2(pb): #based on different theoretical models in Tauris & Savonije 1999
    pars = [(4.5, 1.2e5, 0.12), (4.75, 1.1e5, 0.115), (5.0, 1.e5, 0.11)]
    M2s = np.array([TS99(pb, p) for p in pars])
    return np.random.uniform(M2s.min(), M2s.max())
    
ipb = plist.index('PB')
PB = np.array([float(p[ipb])*secperday for p in MarkovChain])

if 'M2' in plist:
    im2 = plist.index('M2')
    M2 = np.array([float(p[im2]) for p in MarkovChain])
else:
    M2 = np.array([PbToM2(p/secperday) for p in PB])

#print 'M2', M2

if 'SINI' in plist:
    isini = plist.index('SINI')
    SINI = np.array([float(p[isini]) for p in MarkovChain])
elif 'KIN' in plist:
    ii = plist.index('KIN')
    SINI = np.sin(np.array([float(p[ii])/180.*np.pi for p in MarkovChain]))
else:
    COSI = np.random.uniform(0., 1., MCMCSize)
    SINI = np.sqrt(1. - COSI**2)

#print 'SINI', SINI


if 'E' in plist:
    iecc = plist.index('E')
    ECC = np.array([float(p[iecc]) for p in MarkovChain])
elif 'ECC' in plist:
    iecc = plist.index('ECC')
    ECC = np.array([float(p[iecc]) for p in MarkovChain])
elif 'EPS1' in plist:
    ie1 = plist.index('EPS1')
    ie2 = plist.index('EPS2')
    E1 = np.array([float(p[ie1]) for p in MarkovChain])
    E2 = np.array([float(p[ie2]) for p in MarkovChain])
    ECC = np.sqrt(E1**2 + E2**2)


if 'PAASCNODE' in plist:
    iomega = plist.index('PAASCNODE')
    PAASCNODE = np.array([float(p[iomega]) for p in MarkovChain])
elif 'KOM' in plist:
    iomega = plist.index('KOM')
    PAASCNODE = np.array([float(p[iomega]) for p in MarkovChain])
else:
    PAASCNODE = np.random.uniform(0., 360.,  MCMCSize)


#M1 = (PB/2/pi*np.sqrt(G*(M2*SINI)**3/a**3)-M2)/Msun
M1 = PB/2/pi*(np.sqrt(Tsun*(M2*SINI)**3/A**3))-M2
#print 'M1', M1


#im2 = plist.index('M2')
#M2 = np.array([float(p[im2]) for p in MarkovChain])
#M1 = (PB/2/pi*np.sqrt(Tsun*(M2*SINI)**3/A**3)-M2)
#iomega = plist.index('PAASCNODE')
#PAASCNODE = np.array([float(p[iomega]) for p in MarkovChain])


Dmean = 1/PX.mean()
wy_m = PMRA.mean()*1.e-3/60./60.*np.pi/180./secperyear * Dmean 
wz_m = PMDEC.mean()*1.e-3/60./60.*np.pi/180./secperyear * Dmean 
w_m = np.sqrt(wy_m**2 + wz_m**2) * kpc
#print 'mean proper motion', w_m * kpc
w = np.array([ y for y in npr.normal(0., w_m, MCMCSize)])
F0 = float(pf.F0[0])

ichisq = plist.index('chisq')
chisq = [p[ichisq] for p in MarkovChain]
bestidx = chisq.index(min(chisq))

def Integ(alpha):
    global M1, M2, PB, F0, ECC, PMRA, PMDEC, PX, SINI, PAASCNODE, OM, w
    VI = np.arange(M1.size)[np.logical_and(M1 > 1.0, M1 < 2.5)] #valid indices
    #print 'validindices', VI
    #pres = np.array([Pintegrant(M1[i], M2[i], PB[i], F0[i], ECC[i], PMRA[i], PMDEC[i], PX[i], SINI[i], PAASCNODE[i], OM[i], w[i]) for i in VI])
    pres = np.array([Pintegrant(M1[i], M2[i], PB[i], F0, ECC[i], PMRA[i], PMDEC[i], PX[i], SINI[i], PAASCNODE[i], OM[i], w[i]) for i in VI])
    THETA = 1.* pres[...,0]
    EF = alpha * pres[...,1]
    Area1 = sum(EccArea(ECC[VI], EF, THETA))
    Area2 = sum(EccArea(ECC[VI], EF, np.pi-THETA))
    Areas = Area1 + Area2
    return Areas


alpha = np.arange(5.e-22, 1.e-19, 5.e-22) #ingrid setting
res = np.array([Integ(a) for a in alpha])

np.save(open('alpha3.npy', 'wb'), res)

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
