from math import *
import os,sys
from psrpbdot import Pbdot_Gal, Shlkovskii, Decimal, Pbdot_GW
from round import TexStyle as SF

c = 2.99792458e10
PI = pi
AU = 1.469e13
Msun = 1.9882e33
secperday = 24*3600
Tsun = 4.9225e-6

def M2(pf):
    return float(pf.M2[0])*Msun

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


Sp = lambda pf: 0.1*M1(pf)

psr = {'P':1/218.81184391573209821,'Pdot':(4.0833313846102338808e-16)/(218.81184391573209821)**2,'Pb':(67.825129921577592219)*secperday, 'M1':(1.4)*Msun, 'M2':(0.3)*Msun}


from datatools.tempo import tempofit, tempo2fit, touchparfile, uniquename, PARfile

#pf = PARfile('./1713.PBDOT.par')
#pf = PARfile('./1713.guppi2.PBDOT.par')
#pf = PARfile('./1713.ext.pbdot.par.t2')
#pf = PARfile('./1713.ext.par.t2')
#pf = PARfile('./1713.omdot.par.t2')
#pf = PARfile('./1713.ext.PBDOT.par')
#pf = PARfile('./1713.noguppi.par.t2.Jul30')
#pf = PARfile('./1713.noguppi.par.mcmc')
pf = PARfile('./1713.ext.FD.par.t2')
#pf = PARfile('./1713.DMX.final.par.t2')


def Pbdot_exc(psr, GdotOG, KD):
    GdotOG = GdotOG/secperday/365
    M1 = psr['M1']
    M2 = psr['M2']
    q = M1/M2
    return -2.* GdotOG * (1 - (1 + M2/2./(M1+M2))*psr['Sp']) - 4. * PI**2 *Tsun*M2 * KD * psr['Sp']**2 * q/(q+1)/psr['Pb']**2

J1713 = {'M1':M1(pf), 'M2':M2(pf)/Msun, 'Sp':Sp(pf), 'Pb':float(pf.PB[0])*secperday}
#print J1713
#print Pbdot_exc(J1713, 6.e-13, 2.e-4)
#print '1:', Pbdot_exc(J1713, 1.52789743984e-14, 0.00163158650413)
J1012 = {'M2':0.16, 'M1':0.16*10.5, 'Pb':0.60467271355*secperday}
J1012['Sp'] = 0.1*J1012['M1']
#print J1012
#print '2:',Pbdot_exc(J1012, 1.52789743984e-14, 0.00163158650413)
#print -0.4e-14/0.60467271355/secperday
#print pf.PBDOT[0]/pf.PB[0]/secperday


#PbdOPb1713 = 1.5e-13/float(pf.PB[0])/secperday
#PbdOPb1713 = 6.298516276376e-14/float(pf.PB[0])/secperday
if pf.PBDOT[0] > Decimal('1.e-10'):
    pf.PBDOT[0] *= Decimal('1.e-12')
    pf.PBDOT[1] *= Decimal('1.e-12')
Pbdot_exc_1713 = pf.PBDOT[0] - Pbdot_Gal(pf) - Decimal(Shlkovskii(pf)) - Pbdot_GW(pf)
Pbdot_Shl = float(Shlkovskii(pf)) #* 1.e13
Pbdot_Shl_err = float(pf.PX[1]/pf.PX[0]) * Pbdot_Shl
Pbdot_Gal_val = float(Pbdot_Gal(pf))  #* 1.e13
Pbdot_Gal_err = abs(sqrt(4*(0.9/27.2)**2 + (0.4/8.0)**2) * Pbdot_Gal_val)
Pbdot_exc_err = sqrt((float(pf.PBDOT[1]))**2 + Pbdot_Shl_err**2 + Pbdot_Gal_err**2)
print 'PSR ', pf.PSRJ
#print 'Pbdot_obs:', SF((float(pf.PBDOT[0])*1.e13, float(pf.PBDOT[1])*1.e13))
print 'Pbdot_obs:', SF((float(pf.PBDOT[0]), float(pf.PBDOT[1])))
print 'Pbdot_Shl:',  SF((Pbdot_Shl, Pbdot_Shl_err))
print 'Pbdot_Gal:', SF((Pbdot_Gal_val,Pbdot_Gal_err))
#print 'Pbdot_exc:', SF((float(Pbdot_exc_1713)*1.e13, Pbdot_exc_err))
print 'Pbdot_exc:', SF((float(Pbdot_exc_1713), Pbdot_exc_err))
print 'Pbdot_GW:', float(Pbdot_GW(pf)) 
#sys.exit(0)
PbdOPb1713 = float(Pbdot_exc_1713)/float(pf.PB[0])/secperday
PbdOPberr1713 = float(pf.PBDOT[1])/float(pf.PB[0])/secperday

print 'Pbdot/Pb [1713]:', SF((PbdOPb1713, PbdOPberr1713))

sys.exit(0)

PbdOPb1012 = -0.4e-14/0.60467271355/secperday
PbdOPberr1012 = 1.6e-14/0.60467271355/secperday/2 #2sigma 95 error bar needs to be reduced to 1 sigam

#print PbdOPb1012, PbdOPberr1012

#sys.exit(0)

pf = PARfile('./J0437.par')
J0437 = {'M1':M1(pf), 'M2':M2(pf)/Msun, 'Sp':Sp(pf), 'Pb':float(pf.PB[0])*secperday}
PbdOPb0437 = 1.04e-19
#PbdOPb0437 = 3.2e-19
PbdOPberr0437 = 5.7e-19 #2sigma 95 error bar needs to be reduced to 1 sigam

pf = PARfile('./1738+03.par')
pf.PBDOT[0] = pf.PBDOT[0]*Decimal('1.e-15')
pf.PBDOT[1] = pf.PBDOT[1]*Decimal('1.e-15')
J1738 = {'M1':M1(pf), 'M2':M2(pf)/Msun, 'Sp':Sp(pf), 'Pb':float(pf.PB[0])*secperday}
#print J1738['M1']
#sys.exit(0)
PbdOPb1738 = float(pf.PBDOT[0] -Pbdot_Gal(pf)-Decimal(Shlkovskii(pf))-Pbdot_GW(pf) )/float(pf.PB[0])/secperday
#PbdOPb0437 = 3.2e-19
#PbdOPberr1738 = float(pf.PBDOT[1]) /float(pf.PB[0])/secperday #2sigma 95 error bar needs to be reduced to 1 sigam
PbdOPberr1738 = 3.65e-15/ float(pf.PB[0])/secperday #Freire et al .2012


from scipy.stats import chisqprob
from math import *
import os
from copy import *
from numpy.random import normal , uniform ,seed
import numpy as np

from itertools import count
from pylab import *
import pickle 


def probcal(GdotOG, KD):
    #try:
    #chisq = (Pbdot_exc(J1713, GdotOG, KD) - PbdOPb1713)**2/PbdOPberr1713**2 + (Pbdot_exc(J1012, GdotOG, KD) - PbdOPb1012)**2/PbdOPberr1012**2 
    #chisq = (Pbdot_exc(J1012, GdotOG, KD) - PbdOPb1012)**2/PbdOPberr1012**2 + (Pbdot_exc(J0437, GdotOG, KD) - PbdOPb0437)**2/PbdOPberr0437**2
    #chisq = (Pbdot_exc(J1713, GdotOG, KD) - PbdOPb1713)**2/PbdOPberr1713**2 + (Pbdot_exc(J1012, GdotOG, KD) - PbdOPb1012)**2/PbdOPberr1012**2 + (Pbdot_exc(J0437, GdotOG, KD) - PbdOPb0437)**2/PbdOPberr0437**2
    chisq = (Pbdot_exc(J1713, GdotOG, KD) - PbdOPb1713)**2/PbdOPberr1713**2 + (Pbdot_exc(J1012, GdotOG, KD) - PbdOPb1012)**2/PbdOPberr1012**2 + (Pbdot_exc(J0437, GdotOG, KD) - PbdOPb0437)**2/PbdOPberr0437**2 + (Pbdot_exc(J1738, GdotOG, KD) - PbdOPb1738)**2/PbdOPberr1738**2
    #except:
        #print  GdotOG, KD
        #print J1713, J1012
        #print PbdOPb1713, PbdOPb1012
        #print f(J1713, GdotOG, KD)  , f(J1012, GdotOG, KD)  
    return exp(-1.*chisq)


cwd=os.getcwd()
#best = (-0.7e-12, 0.3e-3)
#dict = {}
#dict['BEST'] = best
#pickle.dump(dict, open('%s/bestpar.p' % cwd, 'w'))

from ProgressBar import progressBar
def mcmc(Chain, runtime, mixingtime=1000):
    #mixingtime = 1000
    #runtime = 50000
    pb = progressBar(maxValue = runtime)
    dict = pickle.load(open('%s/bestpar.p' % cwd, 'r'))
    GdotOG, KD = dict['BEST']
    MarkovChain = Chain
    best = (GdotOG, KD)
    p = probcal(GdotOG, KD)
    p0 = p
    n = count()
    m = count()
    while n.next() < mixingtime + runtime:
        NKD = KD + normal(0, 1.e-3)
        NGdotOG = GdotOG + normal(0, 1.e-12)
        p1 = probcal(NGdotOG, NKD)
        c = m.next()
        if p1 > p0:
            if c > mixingtime:
                MarkovChain.append((NGdotOG, NKD ))
                pb(c - mixingtime)
            p0 = p1
            KD, GdotOG = NKD, NGdotOG
            if p1 > p:
                best = (NGdotOG, NKD)
                p = p1
                dict['BEST'] = best
                pickle.dump(dict, open('%s/bestpar.p' % cwd, 'w'))
        else:
            t = uniform(0,1,1)[0]
            if t < p1/p0:
                if c > mixingtime:
                    MarkovChain.append((NGdotOG, NKD))
                    pb(c - mixingtime)
                p0 = p1
                KD, GdotOG = NKD, NGdotOG
            else:
                if c > mixingtime:
                    MarkovChain.append((GdotOG, KD))
    #print  MarkovChain
    print 'Best Gdot/G, Kappa_D:', best
    print '%d new points generated.' %len(MarkovChain)
    #OldChain = pickle.load(open('MChain.p','r'))['Chain']
    #dict['Chain'] = OldChain + MarkovChain
    #dict['Chain'] = MarkovChain
    #os.rmdir(tmpdir)
    #os.chdir(cwd)


class MChain(object):
    def __enter__(self):
        try:
            Chain = pickle.load(open('MChain.p', 'r'))['Chain']
        except:
            Chain = []
        self.Chain = Chain
        #self.cwd = os.getcwd()
        return self.Chain
    def __exit__(self, exc_type, exc_value, exc_tb):
        #os.chdir(self.cwd)
        MarkovChain = self.Chain
        dict = {'Chain':MarkovChain}
        pickle.dump(dict, open('MChain.p', 'w'))

        if exc_type is KeyboardInterrupt:
            print '\nManually Stopped\n'
            return True
        else:
            return exc_type is None
        print '\nFinish running\n' 



if __name__ == '__main__':
    def run(s):
        seed(s) # assigning different initial seed for the random number generator in different threads.
        steps = 5000000
        with MChain() as Chain:
            #print steps, tmpdir
            mcmc(Chain, steps )
        return Chain

    #from multiprocessing import Pool
    #p = Pool(4)
    #p.map(run, range(4))
    run(10)

    #dict = pickle.load(open('bestpar.p', 'r'))
    #best = dict['BEST']
    #MarkovChain = pickle.load(open('MChain.p','r'))['Chain']
    #GdotOG = [x[0] for x in MarkovChain]
    #KD = [x[1] for x in MarkovChain]
    #xlabel('$\dot{G}/G$')
    #ylabel('$\kappa_D$')
    #plot(GdotOG, KD, 'b-')
    #plot([best[0]], [best[1]], 'ro')
    #show()


