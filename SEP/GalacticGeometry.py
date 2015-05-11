from psrpbdot import M1
from datatools.tempo import *
from astropy import coordinates as coord
import numpy.linalg as linalg
from tools import Coordinate as COORD 
from Arrow3D import Arrow3D

G = Decimal(6.673e-8)
c = Decimal(2.99792458e10)
PI = Decimal(np.pi)
AU = Decimal(1.469e13)
Msun = Decimal(1.9882e33)
secperday = 24*3600
solardist = 8.34
#solardist = 7.94


def getGpos(pf):
    ra = COORD.RA(pf.RAJ[0])
    dec = COORD.Dec(pf.DECJ[0])
    pos = coord.SkyCoord(str(ra) +' '+ str(dec))
    l = pos.galactic.l.rad
    b = pos.galactic.b.rad
    return l,b

def GetTransferMatrix(pf):#, PAASCNODE):
    #print getGpos(pf)
    gl, gb = getGpos(pf)
    try:
        D = 1./float(pf.PX[0])
    except:
        D = float(pf.Dist[0])
    x = D * np.cos(gb) * np.sin(gl)
    y = solardist - D * np.cos(gb) * np.cos(gl) 
    angle = np.arctan(y/x)*180./np.pi

    GCpos = coord.SkyCoord('17h45m40.04s -29d00m28.1s')
    pos1713 = coord.SkyCoord(str(COORD.RA(pf.RAJ[0]))+' '+str(COORD.Dec(pf.DECJ[0])))
    RA = pos1713.ra.rad
    Dec = pos1713.dec.rad
    DRA = pos1713.ra.rad - GCpos.ra.rad
    DDec = pos1713.dec.rad - GCpos.dec.rad
    #Omega = float(pf.PAASCNODE)/180.*np.pi
    #Omega = float(PAASCNODE)/180.*np.pi
    Theta_g = np.pi - np.arctan(np.tan(DRA)/np.sin(DDec))
    """
    #Transfer the RA-Dec coordiate to the plane of sky coordiate: Two steps:
    # 1. rotate to the star
    # 2. rotate to the plane of the sky
    #two angles: alpha, beta; alpha = Ra - pi/2, beta = dec
    """
    RA = pos1713.ra.rad
    Dec = pos1713.dec.rad
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
    """
    #Transfer the Galactic plane coordinate to RA-Dec coordinate: Solve the matrix from two equations:
    # 1. GC is (cos Rc cos Dc, sin Rc cos Dc, sin Dc)
    # 2. Galactic pole is (cos Rp cos Dp, sin Rp cos Dp, sin Dp)
    #two angles: alpha, beta; alpha = Ra - pi/2, beta = dec
    """
    GPpos = coord.SkyCoord('12h51m26.282s' +' '+ '27d07m42.01s')
    #print Galactic_pol.galactic.l.rad/np.pi, Galactic_pol.galactic.b.rad/np.pi
    Rc, Dc = GCpos.ra.rad, GCpos.dec.rad
    Rp, Dp = GPpos.ra.rad, GPpos.dec.rad
    Gc = np.matrix((np.cos(Rc)*np.cos(Dc),np.sin(Rc)*np.cos(Dc),np.sin(Dc))).T
    Gp = np.matrix((np.cos(Rp)*np.cos(Dp), np.sin(Rp)*np.cos(Dp), np.sin(Dp))).T
    #Gmat = np.vstack((Gp, Gc, np.matrix((1,1,1))))
    #RHS = np.matrix((0, coszeta, 1))
    #RHS = RHS.T
    #X = linalg.solve(Gmat, RHS)

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

    #sys.exit(0)

    return T, GT


if __name__ == '__main__':
    #pf = PARfile('1713.Dec.mcmc.par')
    pf = PARfile('1713.sns.par')

    gl, gb = getGpos(pf)
    D = 1./float(pf.PX[0])
    x = D * np.cos(gb) * np.sin(gl)
    y = solardist - D * np.cos(gb) * np.cos(gl) 
    angle = np.arctan(y/x)*180./np.pi

    GCpos = coord.SkyCoord('17h45m40.04s -29d00m28.1s')
    pos1713 = coord.SkyCoord(str(COORD.RA(pf.RAJ[0]))+' '+str(COORD.Dec(pf.DECJ[0])))

    RA = pos1713.ra.rad
    Dec = pos1713.dec.rad
    DRA = pos1713.ra.rad - GCpos.ra.rad
    DDec = pos1713.dec.rad - GCpos.dec.rad
    Omega = float(pf.PAASCNODE)/180.*np.pi
    Theta_g = np.pi - np.arctan(np.tan(DRA)/np.sin(DDec))

    T, GT = GetTransferMatrix(pf)
    z = np.sin(gb) * D
    R0 = solardist
    R1 = np.sqrt(R0**2 + (D*np.cos(gb))**2 -2 * R0 * D * np.cos(gb) * np.cos(gl))
    #lbd = np.arccos((R1**2 + D**2 + z**2 - R0**2)/(2*D*np.sqrt(R1**2 + z**2)))#*180/np.pi
    #print 'lambda: ', lbd, 'R_PSR: ', R1, 'z: ', z
    print  'R_PSR: ', R1, 'z: ', z
    #rat = np.sqrt(1 - (np.cos(i)*np.cos(lbd) + np.sin(i)*np.sin(lbd)*np.sin(Theta_g - Omega))**2)
    #print 'ratio:',rat
    Mtot = M1(pf) + pf.M2[0]
    print 'M1+M2', Mtot
    Kz = lambda z:(2.27*z + 3.68*(1-np.exp(-4.31*z)) ) * 1.e-9 #Galactic acceleration in z direction (cm/s^2)
    Omega_G = 27.2 #km s^-1 kpc^-1
    #R_G = 8.33 # +/-0.35 kpc
    kpcinkm = 3.0857e16
    #Kr =  Omega_G**2 * R1 / kpcinkm * 1.e5 #Galactic acceleration in radio direction (cm/s^2)
    Kr =  Omega_G**2 * R0**2 / R1 / kpcinkm * 1.e5 #Galactic acceleration in radio direction (cm/s^2)
    #print 'Kz, Kr:', Kr/Kz(z), R1/z
    KG = np.sqrt(Kr**2 + Kz(z)**2) 
    #print 'Galactic acceleration: ',Kz(z), Kr, KG
    print 'Galactic acceleration: ', KG
    #print pol_ra,pol_dec 
    #sys.exit(0)
    """ 
    zeta: the angle between the pulsar's radial galactic accelration and the direction of Galactic center.
    """
    coszeta = (R0**2 + R1**2 - D**2 + z**2)/R0/R1/2.
    zeta = np.arccos(coszeta)
    print 'coszeta: ', coszeta, np.arccos(coszeta) 

    print 'zeta: ', zeta*180./np.pi
    #g_r = GT.I * Kr * ( np.matrix((np.sin(0.-zeta),np.cos(0.-zeta),0)).T) #X/linalg.norm(X)
    g_r = GT.I * Kr * ( np.matrix((np.cos(0.-zeta),np.sin(0.-zeta),0)).T) #X/linalg.norm(X)
    #g_z = -1. * Kz(z) * Gc
    g_z = GT.I * Kz(z) * (np.matrix((0., 0., -1.)).T) 
    #print g_r, g_z
    g = g_r + g_z
    #print 'GT:', GT
    #print (g_r.T * g_z)/linalg.norm(g)
    #print g
    g_NSEW = T * g
    print 'g:', g.T
    #print 'T:', T
    print 'g_NSEW:', g_NSEW.T
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
    print 'n_orb:', n_orb.T

    """
    Calculte the projection of g on the orbit plane
    """
    g_proj = g_NSEW - n_orb * (g_NSEW.T*n_orb) 
    #print g_proj
    #print  n_orb.T * g_proj
    #print Omgang * 180/np.pi
    KG  = linalg.norm(g_proj)
    g_dir = g_proj/KG
    A_ref = np.matrix((0, -1.* np.sin(Omgang), np.cos(Omgang)))
    #print g_dir, A_ref,A_ref * g_proj 
    #print 'A_ref,', A_ref, 'g_dir', g_dir
    g_ang = np.arccos(A_ref * g_dir)
    print 'g angle on orbit:', g_ang*180/np.pi, A_ref * g_dir
    print 'angle between g and periastron: ', g_ang*180./np.pi ,float(pf.OM[0])
    print 'theta:', g_ang*180./np.pi - float(pf.OM[0])
    ratio = KG/linalg.norm(g)
    print 'projected Galactic acceleration: ', KG, ratio

    import numpy as np
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    Xstar = np.matrix((3, 0, 0)).T
    Ystar = np.matrix((0, 3, 0)).T
    Zstar = np.matrix((0, 0, 3)).T

    X = T.I * Xstar
    Y = T.I * Ystar
    Z = T.I * Zstar

    GX = GT * X 
    GY = GT * Y 
    GZ = GT * Z 

    from matplotlib.patches import Ellipse, Circle
    import mpl_toolkits.mplot3d.art3d as art3d
    from matplotlib.collections import PatchCollection
    #ells = [Ellipse(xy=[x[i], y[i]], width=(0.6+0.4*rand())*(1.-0.5), height=0.6*1, angle=90, fill=False, lw=2)  for i in range(len(psrs))]

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
        #ax.text(0,0,'GC')
    ax.set_xlabel('x (kpc)')
    ax.set_ylabel('y (kpc)')
    #ax.set_xlim((-4.0,6.8))
    #ax.set_ylim((-0.2,10.6))
    #ax.set_zlim(-1.,1.)
    ax.set_xlim((-11.,11.))
    ax.set_ylim((-11,11.))
    ax.set_zlim(-11.,11.)

    AX = Arrow3D([x,x+GX[1]], [y,y-GX[0]], [z,z+GX[2]], mutation_scale=20, lw=1, arrowstyle="-|>", color="r")
    AY = Arrow3D([x,x+GY[1]], [y,y-GY[0]], [z,z+GY[2]], mutation_scale=20, lw=1, arrowstyle="-|>", color="b")
    AZ = Arrow3D([x,x+GZ[1]], [y,y-GZ[0]], [z,z+GZ[2]], mutation_scale=20, lw=1, arrowstyle="-|>", color="y")

    ax.add_artist(AX)
    ax.add_artist(AY)
    ax.add_artist(AZ)

    Gg = GT*g/KG*5
    GA = Arrow3D([x,x+Gg[1]], [y,y-Gg[0]], [z,z+Gg[2]], mutation_scale=20, lw=1, arrowstyle="-|>", color="m")

    ax.add_artist(GA)


    N_ORB = GT * T.I * n_orb * 4
    ANO = Arrow3D([x, x+N_ORB[1]], [y, y-N_ORB[0]], [z,z+N_ORB[2]], mutation_scale=20, lw=1, arrowstyle="-|>", color="g")
    ax.add_artist(ANO)

    plt.show()

    print 'EF:', Decimal(0.009) * Decimal(0.5) * Decimal(KG) * c**2 / G / Mtot / Msun /(2*PI/pf.PB[0]/secperday)**2



