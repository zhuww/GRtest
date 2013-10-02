#from decimal import *
from math import *
from tools.Coordinate import RA, Dec
from astropy import coordinates as coord

#costants:
c = 2.99792458e10
PI = pi
AU = 1.469e13
Msun = 1.9882e33
secperday = 24*3600

psr = {'P':1/218.81184391573209821,'Pdot':(4.0833313846102338808e-16)/(218.81184391573209821)**2,'Pb':(67.825129921577592219)*secperday, 'M1':(1.4)*Msun, 'M2':(0.3)*Msun}

def Pbdot_Mdot(psr):
    globals().update(psr)
    Ip = float(1.e45)
    c = float(2.99792458e10)
    Pi = float(3.14159265)
    return 8*Pi**2*Ip/c**2/(M1+M2)*Pdot/P**3*Pb


from datatools.tempo import PARfile

pf = PARfile('/home/zhuww/work/1713_data/tempo/1713.PBDOT.par')


def Shlkovskii(pf):
    d = AU * 180 / PI * 3600 / pf.PX[0] * 1000 
    Pb = pf.PB[0] * secperday
    PMRA = pf.PMRA[0] / 3600000 * PI / 180 / secperday / 365
    PMDEC = pf.PMDEC[0] / 3600000 * PI / 180 / secperday / 365
    return float(str((PMRA**2 + PMDEC**2)*Pb*d/c))



def M1(pf):
    G = float(6.673e-11)
    Msun = float(1.98892e30)
    c = float(2.99792458e8)
    m2 = float(pf.M2[0])*Msun
    I = float(pf.KIN[0])/180*PI
    Pb = float(pf.PB[0])*secperday
    a = float(pf.A1[0])*c
    #result = sqrt(930.998*m2**3*Pb**2/a**3) - m2
    return (Pb/2/PI*(sqrt(G*(m2*sin(I))**3/a**3))-m2)/Msun


#from tools.PyATNF import LQatnf
def Pbdot_Gal(pf):
    #name = 'J1713+0747'
    name = pf.PSRJ
    ra = RA(pf.RAJ[0])
    dec = Dec(pf.DECJ[0])
    pos = coord.FK5Coordinates(str(ra) +' '+ str(dec))
    l = pos.galactic.l.radians
    b = pos.galactic.l.radians
    #result = LQatnf(name, ('GL','GB'))
    #l = float(result['GL'][0])/180*PI
    #b = float(result['GB'][0])/180*PI
    Pb = float(pf.PB[0]) * secperday
    d = AU * 180 / PI * 3600 / float(pf.PX[0]) * 1000 
    kpc = 3.08568025e21 #cm
    z_kpc = d*sin(b)/kpc
    a_z = (2.27*z_kpc + 3.681*(1 - exp(-4.31*z_kpc)))*1.e-9 #cm s^-2
    A_z = -1 * a_z *abs(sin(b))/c
    R0 = 8 * kpc # 8.0 +/- 0.4 kpc
    Omega0 = 27.2 * 8 * 1.e5 #27.2 +/- 0.9 kms s^-1
    print Omega0
    beta = float(d/R0) * cos(b) - cos(l)
    A_x = -1/c * cos(b) * (Omega0**2/R0) * (cos(l) + beta/(sin(l)**2 + beta**2))
    return Pb*(A_z + A_x)

def Agal(pf):
    name = pf.PSRJ
    #name = 'J1713+0747'
    ra = RA(pf.RAJ[0])
    dec = Dec(pf.DECJ[0])
    pos = coord.FK5Coordinates(str(ra) +' '+ str(dec))
    l = pos.galactic.l.radians
    b = pos.galactic.l.radians
    #result = LQatnf(name, ('GL','GB'))
    #l = float(result['GL'][0])/180*PI
    #b = float(result['GB'][0])/180*PI
    Pb = float(pf.PB[0]) * secperday
    d = AU * 180 / PI * 3600 / float(pf.PX[0]) * 1000 
    kpc = 3.08568025e21 #cm
    z_kpc = d*sin(b)/kpc
    a_z = (2.27*z_kpc + 3.681*(1 - exp(-4.31*z_kpc)))*1.e-9 #cm s^-2
    #A_z = -1 * a_z *abs(sin(b))/c
    R0 = 8 * kpc # 8.0 +/- 0.4 kpc
    Omega0 = 27.2 * 8 * 1.e5 #27.2 +/- 0.9 kms s^-1
    beta = float(d/R0) * cos(b) - cos(l)
    a_x = -1 * (Omega0**2/R0) * (cos(l) + beta/(sin(l)**2 + beta**2))
    #A_x = -1/c * cos(b) * (Omega0**2/R0) * (cos(l) + beta/(sin(l)**2 + beta**2))
    return a_x, a_z

print Agal(pf) 
print M1(pf)
