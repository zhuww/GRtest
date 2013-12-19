from decimal import *
from math import *
from round import shortform as SF
import tools.Coordinate as COD
from tools.Coordinate import RA, Dec
from astropy import coordinates as coord

#costants:
c = Decimal(2.99792458e10)
PI = Decimal(pi)
AU = Decimal(1.469e13)
Msun = Decimal(1.9882e33)
secperday = 24*3600

infotable = {
"J0437-4715":     {'GL':253.39, 'GB':-41.96, 'PB':5.74104646       },
"J1012+5307":     {'GL':160.35, 'GB':50.86,  'PB':0.60467271355    },
"J1713+0747":     {'GL':28.75,  'GB':25.22,  'PB':67.8251298718    },
"J1738+0333":     {'GL':27.72,  'GB':17.74,  'PB':0.3547907344     },
"J1909-3744":     {'GL':359.73,  'GB':-19.6,  'PB':1.533449450481  }}

#psr = {'P':1/Decimal(218.81184391573209821),'Pdot':Decimal(4.0833313846102338808e-16)/Decimal(218.81184391573209821)**2,'Pb':Decimal(67.825129921577592219)*secperday, 'M1':Decimal(1.4)*Msun, 'M2':Decimal(0.3)*Msun}

def Pbdot_Mdot(psr):
    globals().update(psr)
    Ip = Decimal(1.e45)
    c = Decimal(2.99792458e10)
    Pi = Decimal(3.14159265)
    return 8*Pi**2*Ip/c**2/(M1+M2)*Pdot/P**3*Pb

#print Pbdot_Mdot(psr)

from datatools.tempo import PARfile


def Shlkovskii(pf):
    inc = (90 - Decimal(str(COD.Dec(pf.DECJ).in_unit_degree)))/180*PI
    d = AU * 180 / PI * 3600 / pf.PX[0] * 1000 
    Pb = pf.PB[0] * secperday
    PMRA = pf.PMRA[0] / 3600000 * PI / 180 / secperday / 365
    PMDEC = pf.PMDEC[0] / 3600000 * PI / 180 / secperday / 365 * Decimal(str(sin(inc)))
    #print 'sin: ', sin(inc)
    return float(str((PMRA**2 + PMDEC**2)*Pb*d/c))

#pf.PMRA[0] = Decimal(2.562)
#pf.PMDEC[0] = Decimal(-25.61)
#pf.PX[0] = Decimal(1.22)
#pf.PB[0] = Decimal(0.60467)



def M1(pf):
    G = Decimal(6.673e-11)
    Msun = Decimal(1.98892e30)
    c = Decimal(2.99792458e8)
    m2 = pf.M2[0]*Msun
    if pf.__dict__.has_key('KIN'):
        I = pf.KIN[0]/180*PI
        sini = sin(float(I))
    else:
        sini = pf.SINI[0]
    Pb = pf.PB[0]*secperday
    a = pf.A1[0]*c
    #result = sqrt(930.998*m2**3*Pb**2/a**3) - m2
    return (Pb/2/PI*Decimal(sqrt(G*(m2*Decimal(str(sini)))**3/a**3))-m2)/Msun
    

#print M1(pf)

def Pbdot_GW(pf):
    Pb = pf.PB[0] * secperday
    M2 = pf.M2[0]
    q = M1(pf)/M2
    Tsun = Decimal(4.9225e-6)
    #Pb = Decimal(0.60467)*secperday
    #M2 = Decimal(0.16)
    #q = Decimal(10.5)
    return - 192*PI/5 * (2*PI/Pb)**Decimal(5./3) * (Tsun*M2)**Decimal(5./3) * q / (q+1)**Decimal(1./3)

#print Decimal(2.)**Decimal(1./2)


from tools.PyATNF import Qatnf
def Pbdot_Gal(pf):
    name = pf.PSRJ
    if not name.startswith('J'):name = 'J'+name
    #print name
    #name = 'J1012+5307'
    #name = 'J1713+0747'
    #name = 'J0437-4715'
    #result = Qatnf(name, ('GL','GB', 'PB', 'Dist'))
    #result = infotable[name]
    #l = Decimal(result['GL'])/180*PI
    #b = Decimal(result['GB'])/180*PI
    ra = RA(pf.RAJ[0])
    dec = Dec(pf.DECJ[0])
    pos = coord.FK5Coordinates(str(ra) +' '+ str(dec))
    l = pos.galactic.l.radians
    b = pos.galactic.l.radians
    #print result['GL'], result['GB']
    #print l, b
    Pb = pf.PB[0] * secperday
    #Pb = Decimal(result['PB'][0]) * secperday
    d = AU * 180 / PI * 3600 / pf.PX[0] * 1000 
    kpc = Decimal(3.08568025e21) #cm
    #d = Decimal(result['Dist'][0]) * kpc
    z_kpc = d*Decimal(sin(b))/kpc
    a_z = (Decimal(2.27)*z_kpc + Decimal(3.681)*Decimal(1 - exp(Decimal(-4.31)*z_kpc)))*Decimal(1.e-9) #cm s^-2
    #print a_z
    A_z = -1 * a_z *abs(Decimal(sin(b)))/c
    #print A_z
    R0 = 8 * kpc # 8.0 +/- 0.4 kpc
    #Omega0 = 22000000
    Omega0 = Decimal(27.2 * 8 * 1.e5) #27.2 +/- 0.9 kms s^-1 Feast & Whitelock 1997; Lazaridis et al. 2009
    beta = float(d/R0) * cos(b) - cos(l)
    A_x = -1/c * Decimal(cos(b)) * (Omega0**2/R0) * Decimal(cos(l) + beta/(sin(l)**2 + beta**2))
    #print A_x
    return Pb*(A_z + A_x)

#print Pbdot_GW(pf)
#print Shlkovskii(pf)
#print Pbdot_Gal(pf)
#print Decimal(6.e-13) - (Decimal(Shlkovskii(pf)) + Pbdot_Gal(pf))

def Pbdot_Gdot_cons(pf):
    Pb = pf.PB[0] * secperday
    GdotOG = Decimal(4.e-13) / secperday / 365
    return -2 * GdotOG * Pb
def Pbdot_Gdot_err(pf):
    Pb = pf.PB[0] * secperday
    GdotOG = Decimal(9.e-13) / secperday / 365
    return abs(2 * GdotOG * Pb)

def Pbdot_GD(pf):
    Pb = pf.PB[0] * secperday
    Tsun = Decimal(4.9225e-6)
    M2 = pf.M2[0]
    q = M1(pf)/M2
    return -4*PI**2*Tsun*M2/Pb*q/(q+1)

def GdotOG(pf):
    Pb = pf.PB[0] * secperday
    Sp = Decimal(0.1)*M1(pf)
    M2 = pf.M2[0]
    fac = -2/(1-(1+M2/2/(M1(pf)+M2))*Sp)*Pb /secperday / 365
    #fac = -2/(1)*Pb /secperday / 365
    return Pbdot_res[0]/fac, Pbdot_res[1]/fac, Pbdot_res[2]/fac

if __name__ == '__main__':

    #pf = PARfile('./J0437.par')
    #pf = PARfile('./1713.PBDOT.par')
    #pf = PARfile('./1713.guppi2.PBDOT.par')
    #pf = PARfile('./1713.ext.pbdot.par.t2')
    #pf  = PARfile('./1713.omdot.par.t2')
    #pf  = PARfile('./1713.ext.PBDOT.par')
    #pf  = PARfile('./1738+03.par')
    pf  = PARfile('./1909-3744.par')
    pf.PBDOT[0] = pf.PBDOT[0]*Decimal('1.e-12')
    pf.PBDOT[1] = pf.PBDOT[1]*Decimal('1.e-12')

    #print pf.PMRA, pf.PMDEC

    #print pf.PBDOT[0]
    #print Pbdot_Gal(pf)
    #print Shlkovskii(pf)
    #print Pbdot_Gdot_cons(pf), Pbdot_Gdot_err(pf)
    print 'PBdot(shlkovskii):', Shlkovskii(pf)
    print 'PBdot(Galaxtic):' , Pbdot_Gal(pf)
    print 'PBdot(GW):', Pbdot_GW(pf)
    print "PBdot - PBdot(Galaxtic) - PBdot(shlkovskii) - PBdot(GW), PBdoterr*2:",
    print (pf.PBDOT[0] - Pbdot_Gal(pf) - Decimal(Shlkovskii(pf)) - Pbdot_GW(pf)), pf.PBDOT[1]*2
    print "PBdot , PBdot(Galaxtic) + PBdot(shlkovskii) + PBdot(GW):",
    print pf.PBDOT[0] , Pbdot_Gal(pf) + Decimal(Shlkovskii(pf)) +  Pbdot_GW(pf)
    print "PBdot_excess/PB --> Gdot/G:",
    print SF(( float((pf.PBDOT[0] - Pbdot_Gal(pf) - Decimal(Shlkovskii(pf)) )/ pf.PB[0]/ secperday) , float( pf.PBDOT[1]*2 / pf.PB[0]/ secperday ) ))

    Pbdot_exc_1713 = (pf.PBDOT[0] - Pbdot_Gal(pf) - Decimal(Shlkovskii(pf)) )


    #Pbdot_Gdot = Pbdot_Gdot_cons(pf), Pbdot_Gdot_err(pf)
    #Pbdot_res = Decimal(6.e-13) - (Decimal(Shlkovskii(pf)) + Pbdot_Gal(pf)) - Pbdot_Gdot[0], Decimal(-10.e-13)-Pbdot_Gdot[1], Decimal(2.e-13)+Pbdot_Gdot[1]
    #Pbdot_Gdip = Pbdot_GD(pf)
    #print Pbdot_Gdot
    #print 'result here:', Decimal(6.e-13) - (Decimal(Shlkovskii(pf)) + Pbdot_Gal(pf)), Decimal(-1.e-12), Decimal(2.e-13)
    #print Pbdot_Gdip
    #KDSp2 = Pbdot_res[0]+Pbdot_res[1] / Pbdot_Gdip, Pbdot_res[0]+Pbdot_res[2] / Pbdot_Gdip
    #print KDSp2
    #Sp = Decimal(0.1)*M1(pf)
    #KD = KDSp2[0]/Sp**2, KDSp2[1]/Sp**2
    #print KD


    #print GdotOG(pf)

    #print M1(pf)
    
