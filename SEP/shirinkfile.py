import cPickle as pickle
import numpy as np
import random

#dic = pickle.load(open('bestpar.p', 'rb'))
#best = dic['BEST']
#plist = dic['parameters'] 
MChain = pickle.load(open('MChain.p','rb'))
#MChain = pickle.load(open('SmallMChain.p','rb'))
MarkovChain = MChain['Chain']
MCMCSize = len(MarkovChain)
print MCMCSize
#print MarkovChain

SmallerChain = random.sample(MarkovChain, 10000)
dct = {'Chain':SmallerChain}
pickle.dump(dct, open('TinyMChain.p', 'wb'), protocol=2)
