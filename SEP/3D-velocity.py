import numpy as np
import matplotlib.pyplot as plt
from pylab import rand, show
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

#from Arrow3D import Arrow3D


def getGpos(pf):
    ra = RA(pf.RAJ[0])
    dec = Dec(pf.DECJ[0])
    pos = coord.FK5Coordinates(str(ra) +' '+ str(dec))
    l = pos.galactic.l.radians
    b = pos.galactic.b.radians
    return l,b

#pf = PARfile('1713.Dec.mcmc.par')
pf = PARfile('1713.sns.par')
#pf = PARfile('mcmcresult.par')

from GalacticGeometry import *

T, GT = GetTransferMatrix(pf)
pos1713 = coord.FK5Coordinates(str(COORD.RA(pf.RAJ[0]))+' '+str(COORD.Dec(pf.DECJ[0])))
GCpos = coord.FK5Coordinates('17h45m40.04s -29d00m28.1s')
Eerr = float(pf.E[1])
OMerr = float(pf.OM[1])


#def alpha3(M1, M2, PB, F0, ECC, PMRA, PMDEC, PX, SINI, PAASCNODE, OM, w, xi):
#def Pintegrant(M1, M2, PB, F0, ECC, PMRA, PMDEC, PX, SINI, PAASCNODE, OM, w):
def RunTest(M1, M2, PB, F0, ECC, PMRA, PMDEC, PX, SINI, PAASCNODE, OM, wx):
    G = 6.673e-8
    c = 2.99792458e10
    PI = np.pi
    AU = 1.469e13
    Msun = 1.9882e33
    secperday = 24*3600
    secperyear = secperday*365.24218967
    kpc = 3.0857e21
    Mtot = M1+M2

    wserr = 0.9e5 #Kogut et al. 1993, Fixsen et al 1996, Hinshaw et al. 2009
    wsolar = 369.e5 + npr.randn()*wserr#See ref below (aaa+13 Planck Team: Aghanim, N. et al. 2013. Planck confirms this using a different method)
    lws, bws = (263.99 + npr.randn()*0.14)/180.*np.pi, (48.26 + npr.randn()*0.03)/180.*np.pi
    #w_s = wsolar * (GT.I * np.matrix((np.cos(bws)*np.sin(lws),np.cos(bws)*np.cos(lws),np.sin(bws))).T) #wrong? why did I put it this way?
    w_s = wsolar * (GT.I * np.matrix((np.cos(bws)*np.cos(lws),np.cos(bws)*np.sin(lws),np.sin(bws))).T)
    ws_NSEW = T * w_s
    D = kpc/PX
    
    wx = wx * (np.matrix((-1., 0., 0.)).T)
    wy = PMRA*1.e-3/60./60.*np.pi/180./secperyear * D * (np.matrix((0.,-1.,0.)).T)
    wz = PMDEC*1.e-3/60./60.*np.pi/180./secperyear * D * (np.matrix((0.,0.,1.)).T)

    #print 'wy, wz:', linalg.norm(wy)/1.e5, linalg.norm(wz)/1.e5
    print 'wy, wz:', wy[1]/1.e5, wz[2]/1.e5

    w_ns = (wy + wz)
    w =  wx + w_ns + ws_NSEW
    #print ws_NSEW.T/1.e5
    print 'total velocity:', linalg.norm(w)/1.e5

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
    #w_ecc = np.cos(w_ang - OM/180.*np.pi) * w_leg
    print 'w_psr_orb:', w_leg/1.e5

    EF = lambda w:0.21* M1 * w * (PB**2) * c**2 *F0/24 / PI /G /Mtot/Msun

    theta, ef = float(w_ang - OM/180.*np.pi), EF(w_leg) 
    #return theta, ef
    print 'two velocity vectors:', linalg.norm(wx + wy + wz )/1.e5, linalg.norm(ws_NSEW)/1.e5
    print 'theta:', (theta+2*np.pi)*180./np.pi
    return T.I * w_ns/linalg.norm(w_ns), w_s/linalg.norm(w_s)

if __name__ == '__main__':
    G = 6.673e-8
    c = 2.99792458e10
    PI = np.pi
    Msun = 1.9882e33
    secperday = 24*3600
    secperyear = secperday*365.24218967

    pf = model('1713.sns.par')
    M2 = float(pf.M2[0])
    PB = float(pf.PB[0]*secperday)
    F0 = float(pf.F0[0])
    ECC = float(pf.E[0])
    PMRA = float(pf.PMRA[0])
    PMDEC = float(pf.PMDEC[0])
    PX = float(pf.PX[0])
    PAASCNODE = float(pf.PAASCNODE)
    OM = float(pf.OM[0])
    SINI = float(pf.SINI[0])
    A1 = float(pf.A1[0])
    M1 = float(M1(pf))
    w = npr.normal(0., 2900000.)
    print 'M1, M2, PB, F0, ECC, PMRA, PMDEC, PX, SINI, PAASCNODE, OM'
    print M1, M2, PB, F0, ECC, PMRA, PMDEC, PX, SINI, PAASCNODE, OM
    #a,b = RunTest(M1, M2, PB, F0, ECC, PMRA, PMDEC, PX, SINI, PAASCNODE, OM, w) 
    #print a, b*1.6e-19
    w, ws = RunTest(M1, M2, PB, F0, ECC, PMRA, PMDEC, PX, SINI, PAASCNODE, OM, w) 
    #sys.exit(0)

    import numpy as np
    import matplotlib.pyplot as plt

    from matplotlib.patches import Ellipse, Circle
    import mpl_toolkits.mplot3d.art3d as art3d
    from matplotlib.collections import PatchCollection
    #ells = [Ellipse(xy=[x[i], y[i]], width=(0.6+0.4*rand())*(1.-0.5), height=0.6*1, angle=90, fill=False, lw=2)  for i in range(len(psrs))]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    gl, gb = getGpos(pf)
    D = 1./float(pf.PX[0])
    x = D * np.cos(gb) * np.sin(gl)
    y = solardist - D * np.cos(gb) * np.cos(gl) 
    z = np.sin(gb) * D
    PeriAng = float(pf.OM[0])
    Ang = PeriAng + gl*180./np.pi
    e = Ellipse(xy = [x,y], width=(0.9), height=0.5, angle=Ang, fill=True, lw=1)
    #ax.add_artist(ell)
    #p = PatchCollection([e])
    #ax.add_collection3d(p, zs=[0])
    #p = Circle((1,1),1)
    #ax.add_patch(p)
    ax.add_patch(e)
    art3d.pathpatch_2d_to_3d(e, z=z, zdir=(0.5,0.5,0.5))
    e.set_alpha(1.)
    e.set_edgecolor('r')

    data = np.genfromtxt(open('log_arms.out', 'r'), dtype = [('Num', 'i'), ('n', 'i'), ('x', 'f8'),('y', 'f8')])[1:]
    UniqNum = np.unique(data['Num'])
    for num in UniqNum:
        gx = data['x'][data['Num'] == num]
        gy = data['y'][data['Num'] == num]
        ax.plot(gx,gy,0,'b-')
        #ax.plot([0],[solardist], 0.17, 'yo', ms=15) #OLausen & Kaspi
        ax.plot([0],[solardist], 0., 'yo', ms=5) #OLausen & Kaspi
        ax.plot([0],[0], 0, 'k*', ms=15)
    ax.set_xlabel('x (kpc)')
    ax.set_ylabel('y (kpc)')
    ax.set_xlim((-11.,11.))
    ax.set_ylim((-11,11.))
    ax.set_zlim(-11.,11.)

    #Xstar = np.matrix((3, 0, 0)).T
    #Ystar = np.matrix((0, 3, 0)).T
    #Zstar = np.matrix((0, 0, 3)).T
    #X = T.I * Xstar
    #Y = T.I * Ystar
    #Z = T.I * Zstar
    #GX = GT * X 
    #GY = GT * Y 
    #GZ = GT * Z 
    #AX = Arrow3D([x,x+GX[1]], [y,y-GX[0]], [z,z+GX[2]], mutation_scale=20, lw=1, arrowstyle="-|>", color="r")
    #AY = Arrow3D([x,x+GY[1]], [y,y-GY[0]], [z,z+GY[2]], mutation_scale=20, lw=1, arrowstyle="-|>", color="b")
    #AZ = Arrow3D([x,x+GZ[1]], [y,y-GZ[0]], [z,z+GZ[2]], mutation_scale=20, lw=1, arrowstyle="-|>", color="y")
    #ax.add_artist(AX)
    #ax.add_artist(AY)
    #ax.add_artist(AZ)


    Gw = GT * (w * 4)
    Aw = Arrow3D([x,x+Gw[1]], [y,y-Gw[0]], [z,z+Gw[2]], mutation_scale=20, lw=1, arrowstyle="-|>", color="r")
    lws, bws = (263.99 + npr.randn()*0.14)/180.*np.pi, (48.26 + npr.randn()*0.03)/180.*np.pi
    Gws = 4 * (np.matrix((np.cos(bws)*np.cos(lws),np.cos(bws)*np.sin(lws),np.sin(bws))).T)
    #Gws = GT * (ws * 4)
    Aws = Arrow3D([0,0+Gws[1]], [solardist,solardist-Gws[0]], [0,0+Gws[2]], mutation_scale=20, lw=1, arrowstyle="-|>", color="y")
    ax.add_artist(Aw)
    ax.add_artist(Aws)

    #print 'is this zero:', GX.T * Gw


    #wserr = 0.9e5 #Kogut et al. 1993, Fixsen et al 1996, Hinshaw et al. 2009
    #wsolar = 369.e5 + npr.randn()*wserr#See ref below (aaa+13 Planck Team: Aghanim, N. et al. 2013. Planck confirms this using a different method)
    #lws, bws = 263.99/180.*np.pi, 48.26/180.*np.pi
    #w_s = wsolar * (np.matrix((np.cos(bws)*np.sin(lws),np.cos(bws)*np.cos(lws),np.sin(bws))).T)
    #Gw_s = w_s/linalg.norm(w_s) * 4
    #Aw_s = Arrow3D([0,0+Gw_s[1]], [solardist,solardist-Gw_s[0]], [0,0+Gw_s[2]], mutation_scale=20, lw=1, arrowstyle="-|>", color="b")
    #ax.add_artist(Aw_s)



show()
sys.exit(0)



def EccArea(ECC, EF, THETA):
    global Eerr, OMerr
    THETA[THETA<0]+=(np.pi*2)
    Emax, Emin = ECC +Eerr*3 , ECC-Eerr*3
    Temax, Temin = THETA +OMerr*3, THETA-OMerr*3
    fourpnts = [[(x*np.cos(y)-EF[i]+x*np.sin(y)*1j) for x in (Emax[i], Emin[i]) for y in (Temax[i], Temin[i])] for i in range(len(ECC))]
    areas = []
    for v4 in fourpnts:
        abs4 = [abs(v) for v in v4]
        p4 = [phase(v) for v in v4]
        if np.mean(abs4) < Eerr*3 and abs(np.mean(p4))<OMerr*3:
            #print "hei, it's happening"
            value = max((np.log10(min(max(abs4), 0.05)) - np.log10(max(min(abs4),1.e-6))),0)*np.pi*2
        else:
            value = max((np.log10(min(max(abs4), 0.05)) - np.log10(max(min(abs4),1.e-6))),0)*(max(p4)-min(p4))
        #print abs4, p4
        areas.append(value)
    areas = np.array(areas)
    return areas

def Integ(a):
    global M1, M2, PB, F0, ECC, PMRA, PMDEC, PX, SINI, PAASCNODE, OM, w
    pres = np.array([Pintegrant(M1[i], M2[i], PB[i], F0[i], ECC[i], PMRA[i], PMDEC[i], PX[i], SINI[i], PAASCNODE[i], OM[i], w[i]) for i in range(len(w))])
    THETA = 1.* pres[...,0]
    EF = a * pres[...,1]
    Areas = 0.
    Areas += sum(EccArea(ECC, EF, THETA))
    Areas += sum(EccArea(ECC, EF, THETA-np.pi))
    return Areas

''' load in the MCMC results for Delta estimation '''
import cPickle as pickle
import sys
#from Coordinate import RA, Dec
secperday = 24*3600

dic = pickle.load(open('bestpar.p', 'rb'))
best = dic['BEST']
plist = dic['parameters'] 
MChain = pickle.load(open('SmallMChain.p','rb'))
MarkovChain = MChain['Chain']
MCMCSize = len(MarkovChain)
pi = 3.141592653589793
#G = 6.673e-11
#Msun = 1.98892e30
#c = 2.99792458e8
G = 6.673e-8
c = 2.99792458e10
Msun = 1.9882e33
twopi = 6.283185307179586
fac = 1.536e-16 
if0 = plist.index('F0')
im2 = plist.index('M2')
ipb = plist.index('PB')
isini = plist.index('SINI')
ia = plist.index('A1')
ichisq = plist.index('chisq')
iecc = plist.index('E')
ipx = plist.index('PX')
ipmra = plist.index('PMRA')
ipmdec = plist.index('PMDEC')
iomega = plist.index('PAASCNODE')
iom = plist.index('OM')
M2 = np.array([float(p[im2])*Msun for p in MarkovChain])
PB = np.array([float(p[ipb])*secperday for p in MarkovChain])
SINI = np.array([float(p[isini]) for p in MarkovChain])
a = np.array([float(p[ia])*c for p in MarkovChain])
#print PB, M2, SINI, 
F0 = np.array([float(p[if0]) for p in MarkovChain])
M1 = (PB/2/pi*np.sqrt(G*(M2*SINI)**3/a**3)-M2)/Msun
M2 = M2/Msun
ECC = np.array([float(p[iecc]) for p in MarkovChain])
PX = np.array([float(p[ipx]) for p in MarkovChain])
PMRA = np.array([float(p[ipmra]) for p in MarkovChain])
PMDEC = np.array([float(p[ipmdec]) for p in MarkovChain])
PAASCNODE = np.array([float(p[iomega]) for p in MarkovChain])
OM = np.array([float(p[iom]) for p in MarkovChain])
chisq = [p[ichisq] for p in MarkovChain]
bestidx = chisq.index(min(chisq))
w = np.array([ y for y in npr.normal(0., 2900000., MCMCSize)])
#alpha = np.array([alpha3(M1[i], M2[i], PB[i], F0[i], ECC[i], PMRA[i], PMDEC[i], PX[i], SINI[i], PAASCNODE[i], OM[i], w[i], xi[i]) for i in range(len(w))])



alpha = np.arange(1.e-22, 5.e-19, 2.e-22) #ingrid setting
res = np.array([Integ(a) for a in alpha])
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
#ax.logx()
#ax.hist(delta, bins=50, normed=1)
ax.semilogx(nonposy='clip')
plt.show()
