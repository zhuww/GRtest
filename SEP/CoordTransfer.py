"""
The angular transformations:
"""
from psrpbdot import M1
from datatools.tempo import *
from astropy import coordinates as coord
import numpy.linalg as linalg
from tools import Coordinate as COORD 

pf = PARfile('1713.Dec.mcmc.par')
pos1713 = coord.FK5Coordinates(str(COORD.RA(pf.RAJ[0]))+' '+str(COORD.Dec(pf.DECJ[0])))
GCpos = coord.FK5Coordinates('17h45m40.04s -29d00m28.1s')
RA = pos1713.ra.radians
Dec = pos1713.dec.radians
alpha =  np.pi + RA
beta = Dec
T1 = np.matrix(((np.cos(alpha), np.sin(alpha), 0.), 
               (0.-np.sin(alpha), np.cos(alpha), 0.), 
               (0., 0., 1)))
T2 = np.matrix(((np.cos(beta), 0., 0.-np.sin(beta)),
               (0., 1., 0.),
               (np.sin(beta), 0., np.cos(beta))))
T = T2*T1

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

print T

print GT

