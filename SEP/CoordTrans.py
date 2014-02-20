from psrpbdot import M1
from datatools.tempo import *
from astropy import coordinates as coord
from tools.Coordinate import *
import numpy.linalg as linalg
secperday = 24*3600
solardist = 8.33

def getGpos(pf):
    ra = RA(pf.RAJ[0])
    dec = Dec(pf.DECJ[0])
    pos = coord.FK5Coordinates(str(ra) +' '+ str(dec))
    l = pos.galactic.l.radians
    b = pos.galactic.b.radians
    return l,b

#print 'pf' in locals()
#print 'pf' in globals()
if not 'pf' in locals():
    #pf = PARfile('1713.Dec.mcmc.par')
    pf = PARfile('1713.sns.par')

#print getGpos(pf)
gl, gb = getGpos(pf)
D = 1./float(pf.PX[0])
x = D * np.cos(gb) * np.sin(gl)
y = solardist - D * np.cos(gb) * np.cos(gl) 
angle = np.arctan(y/x)*180./np.pi

GCpos = coord.FK5Coordinates('17h45m40.04s -29d00m28.1s')
pos1713 = coord.FK5Coordinates(str(RA(pf.RAJ[0]))+' '+str(Dec(pf.DECJ[0])))
RA = pos1713.ra.radians
Dec = pos1713.dec.radians
DRA = pos1713.ra.radians - GCpos.ra.radians
DDec = pos1713.dec.radians - GCpos.dec.radians
Omega = float(pf.PAASCNODE)/180.*np.pi
Theta_g = np.pi - np.arctan(np.tan(DRA)/np.sin(DDec))

#print Theta_g*180/np.pi, Omega*180/np.pi

z = np.sin(gb) * D
R0 = solardist
R1 = np.sqrt(R0**2 + (D*np.cos(gb))**2 -2 * R0 * D * np.cos(gb) * np.cos(gl))
lbd = np.arccos((R1**2 + D**2 + z**2 - R0**2)/(2*D*np.sqrt(R1**2 + z**2)))#*180/np.pi
print 'lambda: ', lbd, 'R_PSR: ', R1, 'z: ', z
#rat = np.sqrt(1 - (np.cos(i)*np.cos(lbd) + np.sin(i)*np.sin(lbd)*np.sin(Theta_g - Omega))**2)
#print 'ratio:',rat
Mtot = M1(pf) + pf.M2[0]
print 'M1+M2', Mtot
Kz = lambda z:(2.27*z + 3.68*(1-np.exp(-4.31*z)) ) * 1.e-9 #Galactic acceleration in z direction (cm/s^2)
Omega_G = 27.2 #km s^-1 kpc^-1
#R_G = 8.33 # +/-0.35 kpc
kpcinkm = 3.0857e16
Kr =  Omega_G**2 * R1 / kpcinkm * 1.e5 #Galactic acceleration in radio direction (cm/s^2)
print 'Kz, Kr:', Kr/Kz(z), R1/z
KG = np.sqrt(Kr**2 + Kz(z)**2) 
print 'Galactic acceleration: ',Kz(z), Kr, KG
#print 'Projected Galactic acceleration: ', KG*rat
#print pol_ra,pol_dec 
#sys.exit(0)
""" 
zeta: the angle between the pulsar's radial galactic accelration and the direction of Galactic center.
"""
coszeta = (R0**2 + R1**2 - D**2 + z**2)/R0/R1/2.
print 'coszeta: ', coszeta, np.arccos(coszeta) 

"""
#Transfer the RA-Dec coordiate to the plane of sky coordiate: Two steps:
# 1. rotate to the star
# 2. rotate to the plane of the sky
#two angles: alpha, beta; alpha = Ra - pi/2, beta = dec
"""
RA = pos1713.ra.radians
Dec = pos1713.dec.radians
alpha =  np.pi + RA
beta = Dec
#print 'RA + 180., Dec: ', alpha/np.pi*180, beta/np.pi*180
T1 = np.matrix(((np.cos(alpha), np.sin(alpha), 0.), 
               (0.-np.sin(alpha), np.cos(alpha), 0.), 
               (0., 0., 1)))
T2 = np.matrix(((np.cos(beta), 0., 0.-np.sin(beta)),
               (0., 1., 0.),
               (np.sin(beta), 0., np.cos(beta))))
T = T2*T1
#print T1  
#print T2
#print T
#print T.I
#sys.exit(0)
"""
#Transfer the Galactic plane coordinate to RA-Dec coordinate: Solve the matrix from two equations:
# 1. GC is (cos Rc cos Dc, sin Rc cos Dc, sin Dc)
# 2. Galactic pole is (cos Rp cos Dp, sin Rp cos Dp, sin Dp)
#two angles: alpha, beta; alpha = Ra - pi/2, beta = dec
"""
import numpy.linalg as linalg
GPpos = coord.FK5Coordinates('12h51m26.282s' +' '+ '27d07m42.01s')
Rc, Dc = GCpos.ra.radians, GCpos.dec.radians
Rp, Dp = GPpos.ra.radians, GPpos.dec.radians
Gc = np.matrix((np.cos(Rc)*np.cos(Dc),np.sin(Rc)*np.cos(Dc),np.sin(Dc))).T
Gp = np.matrix((np.cos(Rp)*np.cos(Dp), np.sin(Rp)*np.cos(Dp), np.sin(Dp))).T

alpha = Rc
beta = Dc
GT1 = np.matrix(((np.cos(alpha), np.sin(alpha), 0.), 
               (0.-np.sin(alpha), np.cos(alpha), 0.), 
               (0., 0., 1)))
GT2 = np.matrix(((np.cos(beta), 0., np.sin(beta)),
               (0., 1, 0.),
               (0.-np.sin(beta), 0., np.cos(beta))))
Gp_Gal12 = GT2 * GT1 * Gp
gamma = np.arctan(Gp_Gal12[1]/Gp_Gal12[2])
GT3 = np.matrix(((1., 0., 0.),
                 (0, np.cos(gamma), 0.-np.sin(gamma)),
                 (0, np.sin(gamma), np.cos(gamma))))

GT = GT3 * GT2 * GT1
#GT = GT1 * GT2 * GT3
#Gp_Gal = GT3 * Gp_Gal12
#print Gp_Gal
zeta = np.arccos(coszeta)
#from GalacticGeometry import *

g_r = Kr * (GT.I * np.matrix((np.sin(0.-zeta),np.cos(0.-zeta),0)).T) #X/linalg.norm(X)
g_z = -1. * Kz(z) * Gc
#print g_r, g_z
g = g_r + g_z
#print (g_r.T * g_z)/linalg.norm(g)
#print g
g_NSEW = T * g
#print g_NSEW

"""
#figure out the orbital plane
#define orbital plane as Ax + By + Cz = 0; where A = -1./tg(i), B = -1./tg(Omega), C = 1.
"""
incang = np.arcsin(float(pf.SINI[0]))
Omgang = float(pf.PAASCNODE)/180.*np.pi
A = -1./np.tan(incang)
B = -1./np.tan(Omgang)
C = 1.

n_orb = np.matrix((A, B, C)).T
n_orb= n_orb/linalg.norm(n_orb)


g_proj = g_NSEW - n_orb * (g_NSEW.T*n_orb) 
KG  = float(linalg.norm(g_proj))
print 'Projected Galactic acceleration: ', KG


""" 
The solar system's proper motion in CMB frame is known to be 369+/-0.9 km/s in the direction of (l,b) = 263.99+/-0.14 deg, 48.26+/-0.03 deg (Kogut et al. 1993, Fixsen et al. 1996, Hinshaw et al. 2009)
"""
wsolar = 369. #See ref below (aaa+13 Planck Team: Aghanim, N. et al. 2013. Planck confirms this using a different method)
wserr = 0.9 #Kogut et al. 1993, Fixsen et al 1996, Hinshaw et al. 2009
lws, bws = 263.99/180.*np.pi, 48.26/180.*np.pi

w_s = wsolar * (GT.I * np.matrix((np.cos(bws)*np.sin(lws),np.cos(bws)*np.cos(lws),np.sin(bws))).T)
ws_NSEW = T * w_s
print 'Solar system speed in sky plane frame'
print ws_NSEW

"""
Calculte the projection of g on the orbit plane
"""
#ws_proj = ws_NSEW - n_orb * (ws_NSEW.T*n_orb) 
