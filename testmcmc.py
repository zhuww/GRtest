from datatools.tempo import tempofit, tempo2fit, touchparfile, uniquename, PARfile
from scipy.stats import chisqprob
from math import *
from decimal import *
import os
from copy import *
from numpy.random import normal , uniform ,seed
import numpy as np


#print tempofit('1713.mcmc.par',toafile = '1713.sns.tim')
#print tempo2fit('1713.sns.par.4',toafile = '1713.sns.tim.t2')

#pf = PARfile('1713.mcmc.par')
#pf.write()

#pf.freezeall()
#pf.write('1713.exp.par')
##os.system('tempo -j -f %s %s -no %s' % ('1713.mcmc.par', '1713.sns.tim', 'pulse.tmp'))
#os.system('tempo -j -f %s %s ' % ('1713.mcmc.par', '1713.sns.tim'))
##print tempofit(pf.psrname+'.par', toafile = '1713.sns.tim', pulsefile = 'pulse.tmp')
##print tempofit(pf.psrname+'.par', toafile = '1713.sns.tim')
#pf = PARfile(pf.psrname+'.par')
#while not raw_input() == 'q':
    #new = pf.randomnew()
    #new.freezeall()
    #new.write('1713.exp.par')
    ##print tempofit('1713.exp.par', toafile = '1713.sns.tim', pulsefile='pulse.tmp')
    #print tempofit('1713.exp.par', toafile = '1713.sns.tim')
#pf.freezeall()
#pf.write(pf.psrname+'.par')

#import sys
#sys.exit(0)

#parfile = '1713.mcmc.par'
#toafile = '1713.trim.tim'
parfile = '1713.sns.par'
toafile = '1713.sns.tim'
thetamu = 129.98
x = 32.342420880
mu = 6.0877
twopi = 6.283185307179586
fac = 1.536e-16 * 1.e12
#pf = PARfile(parfile)
def probcal(Omega, cosi, m2, pf):
    pf.write()
    if m2 <= 0 or Omega > 360 or Omega < -360:
        return 0
    sini = sqrt(1 - cosi**2)
    xdot = -1.* fac * x * mu * (cosi/sini) * sin((thetamu-Omega)*twopi/360.)
    touchparfile(parfile, NITS=1, PAASCNODE=Omega, SINI = sini, M2 = m2, XDOT = xdot)
    chisq, dof = tempofit(parfile, toafile = toafile)
    #print chisq, dof
    if chisq >= 99999.:return 0
    #return (dof/chisq)
    #return chisqprob(chisq, dof)
    smallest = 308.28
    #f = 1./sqrt(6.283185307179586)
    #return f*exp(smallest - chisq)*sqrt(chisq - smallest)
    return exp((smallest - chisq)/2.) #Ingrid/Paul?

#print probcal(90, 0.25, 0.28)


from itertools import count
from pylab import *
import pickle 
import os,sys
from tempfile import mkdtemp

#print probcal(iOmega, icosi, im2)

class MChain(object):
    def __enter__(self):
        try:
            Chain = pickle.load(open('MChain.p', 'r'))['Chain']
        except:
            Chain = []
        self.Chain = Chain
        self.cwd = os.getcwd()
        return self.Chain
    def __exit__(self, exc_type, exc_value, exc_tb):
        os.chdir(self.cwd)
        MarkovChain = self.Chain
        dict = {'Chain':MarkovChain}
        pickle.dump(dict, open('MChain.p', 'w'))

        if exc_type is KeyboardInterrupt:
            print '\nManually Stopped\n'
            return True
        else:
            return exc_type is None
        print '\nFinish running\n' 

def motifile(file, cwd, tmpdir):
    os.system('cp %s/%s %s/%s' % (cwd, file, tmpdir, file))
    text = ''
    f = open(file, 'r')
    for l in f.readlines():
        if not l.find('INCLUDE') == -1:
            a = l.split()
            #print a
            l = a[0] + ' '+cwd+'/'+a[1]
            if not l[-1] == '\n':
                l += '\n'
            text += l
        else:
            if not l[-1] == '\n':
                l += '\n'
            text += l
    f.close()
    f = open(file, 'w')
    f.write(text)
    f.close() #motify the tim file to enable those INCLUDE commands.
    
from ProgressBar import progressBar
def mcmc(Chain, runtime, mixingtime=1000):
    #mixingtime = 1000
    #runtime = 50000
    pb = progressBar(maxValue = runtime)
    cwd=os.getcwd()
    tmpdir = cwd+'/'+uniquename()
    if not tmpdir == None:
        if os.path.exists(tmpdir):
            os.chdir(tmpdir)
        else:
            os.mkdir(tmpdir)
            os.chdir(tmpdir)
    os.system('cp %s/%s %s/%s' % (cwd, parfile, tmpdir, parfile))
    motifile(toafile, cwd, tmpdir)
    dict = pickle.load(open('%s/bestpar.p' % cwd, 'r'))
    Omega, cosi, m2 = dict['BEST']
    MarkovChain = Chain
    pf = PARfile(parfile)
    pf.thawall()
    pf.write()
    pf.matrix(toafile)
    pf.freezeall()
    pf.write()
    p0 = probcal(Omega, cosi, m2, pf)
    p = p0
    best = (Omega, cosi, m2 )
    #dict = {'BEST':best}
    #pickle.dump(dict, open('bestpar.p', 'w'))
    n = count()
    m = count()
    while n.next() < mixingtime + runtime:
        npf = pf.randomnew()
        Ncosi = cosi + uniform(-0.05, 0.05)
        #Nm2 = m2 + normal(0, 0.02)
        NOmega = (Omega + normal(0, 5)) 
        #Nsini = float(str(npf.SINI[0]))
        #Ncosi = -1. * sqrt(1 - Nsini**2)
        Nm2 = float(str(npf.M2[0]))
        #print Ncosi, cosi,  Ncosi - cosi, Nm2 - m2
        p1 = probcal(NOmega, Ncosi, Nm2, npf)
        c = m.next()
        if p1 > p0:
            if c > mixingtime:
                MarkovChain.append((NOmega, Ncosi, Nm2 ))
                pb(c - mixingtime)
            pf = npf
            p0 = p1
            cosi, m2, Omega = Ncosi, Nm2, NOmega
            if p1 > p:
                best = (NOmega, Ncosi, Nm2)
                p = p1
                dict['BEST'] = best
                pickle.dump(dict, open('%s/bestpar.p' % cwd, 'w'))
        else:
            t = uniform(0,1,1)[0]
            if t < p1/p0:
                if c > mixingtime:
                    MarkovChain.append((NOmega, Ncosi, Nm2 ))
                    pb(c - mixingtime)
                pf = npf
                p0 = p1
                cosi, m2, Omega = Ncosi, Nm2, NOmega
            else:
                if c > mixingtime:
                    MarkovChain.append((Omega, cosi, m2))
    #print  MarkovChain
    print best
    print '%d new points generated.' %len(MarkovChain)
    #OldChain = pickle.load(open('MChain.p','r'))['Chain']
    #dict['Chain'] = OldChain + MarkovChain
    #dict['Chain'] = MarkovChain
    #os.rmdir(tmpdir)
    os.chdir(cwd)

        
if __name__ == '__main__':
    def run(s):
        seed(s) # assigning different initial seed for the random number generator in different threads.
        steps = 50000
        with MChain() as Chain:
            #print steps, tmpdir
            mcmc(Chain, steps )
        return Chain

    from multiprocessing import Pool
    p = Pool(3)
    p.map(run, range(24,27))
    #run(10)
    dict = pickle.load(open('bestpar.p', 'r'))
    best = dict['BEST']
    MarkovChain = pickle.load(open('MChain.p','r'))['Chain']
    Omega = [x[0] for x in MarkovChain]
    cosi = [x[1] for x in MarkovChain]
    m2 = [x[2] for x in MarkovChain]
    xlabel('$\cos i$')
    ylabel('$M_2$')
    plot(cosi, m2, 'b-')
    plot([best[1]], [best[2]], 'ro')
    show()

