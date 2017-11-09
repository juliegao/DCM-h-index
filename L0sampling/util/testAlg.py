'''
Created on Aug 25, 2017

@author: Ruoyuan Gao
'''
import sys
from algorithms.f0 import F0_track
sys.path.append('/grad/users/rg522/workspace/L0sampling/L0sampling/')

from util.preprocess import SyntheticData
from util.config import NUM_DISTINCT, NUM_SITES_SYN,\
    FILE_ANS_SYN, DELTA, EPS, FILE_SHARED_SYN, NUM_ITEMS_SYN, SPARSITY, C_RAND
from applications.hindex import ExpHist, L0sampling, L0ThreshMonitor
from algorithms.L0_stream import L0
from operator import indexOf
from algorithms.countDown import CountTracking
from math import log, ceil
from random import randint, shuffle

def test_ANS_hindex():
    syn_data = SyntheticData()
    syn_data.loadFromFile(FILE_ANS_SYN)
    stream = syn_data.get_stream()
    
    for i in xrange(1,11): 
        endCountItems = i * NUM_DISTINCT / 10
        expHist = ExpHist(stream, EPS, endCountItems, 
            NUM_SITES_SYN, NUM_DISTINCT, DELTA)
        expHist.run()
        expHist.computeErr()
        print '#items, #msgs: %s %s' %(
            expHist.getCountItems(), expHist.getCountMsg())

def test_ExpCD():
    syn_data = SyntheticData()
    syn_data.loadFromFile(FILE_SHARED_SYN)
    stream = syn_data.get_stream()
    
    for _ in xrange(10):
        itemID = randint(0,NUM_DISTINCT-1)
        print 'item=', itemID
        for i in xrange(5000, 10000, 1000): 
            endCountItems = i * NUM_DISTINCT
            track = CountTracking(stream, EPS, endCountItems, 
                                  NUM_SITES_SYN, NUM_ITEMS_SYN, DELTA,
                                  [itemID, itemID])
            track.run()
            msg = track.getCountMsg()
            print 'i=%s, msgs=%s' %(i, msg)
            track.computeErr()
        
def test_L0():
    '''
    set k=100, N=10**3, 
    sample size S=25, 
    repetitions=5*(10**5) 
    output accuracy measured by 
    max_i |fi - f_star| / f_star
    
    given delta = 0.1
    delta = 2 ** (-s / 12)
    => sparsity = log(1/delta) * 12 =36
    
    '''
    n = 10**3
    r = 10**4
    S = 25 # = 1/eps^2
    f_star = 1.0*r/n
    
#     freq_list_bib = [0] * n
#     print 'bib draws'
#     for i in range(r*S):
#         if i % 10**5 == 0:
#             print 'iter=', i
#         bib_draw = random.randint(0, n-1)
#         freq_list_bib[bib_draw] += 1
#     acc_bib = (
#         abs(1.0*max(freq_list_bib) - f_star)
#         / f_star)
#     del freq_list_bib
#     print 'accuracy of bib=%.2f\n' % acc_bib
    
#     s = int(ceil(log(1/DELTA, 2) / (EPS**2)))  # 75
#     k = int(ceil(log(1/DELTA)))    # 3
    print 'L0 draws'
    countFail = 0
    stream = range(n)
    freq_list_l0 = [0] * n
    
#     # test uniformality
    seed_list = range(r*S)
    shuffle(seed_list)
    seed_i = 0
#     for i in range(r):
#         if i % 10**3 == 0:
#             print 'iter=', i
# #         l0_sampler = L0(n, 24, 12, seed)
#         for _ in range(S):
#             l0_sampler = L0(n, 4, 2, seed_list[seed_i])
#             for item in stream:
#                 l0_sampler.sample(item)
#             sample = l0_sampler.recover_lv(-1)
#             seed_i += 1
#             if sample:
#                 freq_list_l0[sample[0]] += 1
#                 break
#         if not sample:
#             countFail += 1
            
    
    # test sample size
    for i in range(r):
        if i % 10**3 == 0:
            print 'iter=', i
        l0_sampler = L0(n, 24, 12, seed_list[seed_i])
        for item in stream:
            l0_sampler.sample(item)
        samples = l0_sampler.recover_s()
        if samples:
            print '%s samples: %s' % ( len(samples), samples)
            for sample in samples:
                freq_list_l0[sample] += 1
        seed_i += 1

    max_freq = max(freq_list_l0)
    print 'most sampled item %s=freq_%s'%(
         indexOf(freq_list_l0, max_freq), max_freq)
    acc_l0 = ((max_freq - f_star) / (f_star))
    print 'accuracy of l0=%.2f' % acc_l0
    countFail = r*S - sum(freq_list_l0)
    print 'failure rate %s/%s=%s' % (countFail, r, 1.0*countFail/r/S)

def test_L0_hindex():
    syn_data = SyntheticData()
    syn_data.loadFromFile(FILE_SHARED_SYN)
    stream = syn_data.get_stream()
    del syn_data

# #     for i in xrange(100, 500, 100):
#     for i in xrange(50, 100, 10):
# #     for i in xrange(1):
#         endCountItems = i * NUM_DISTINCT
#         endCountItems = 500
#         l0_sampler = L0sampling(stream, EPS, endCountItems,
#                                 NUM_SITES_SYN, NUM_DISTINCT, DELTA)
#         l0_sampler.computeErr()
#         msg = l0_sampler.getCountMsg()
    
    sample_size = int(ceil((3.0*log(2/DELTA, 2) / (EPS**2))))
    numCDs = int(ceil(log(1.0*sample_size/(DELTA*(EPS**2)), 2)))
    boundRand = C_RAND* (1.0/(EPS**2) - 0.5/EPS)
    boundRand = int(boundRand * numCDs)
    detBound = int(2*NUM_SITES_SYN*log(1.0/(EPS**2),2))
    
    # fix k and eps, change n
    print 'numDistinct N=', NUM_DISTINCT
    print 'sample size: %s' %(int(ceil((3*log(2/DELTA, 2) / (EPS**2)))))

    for threshold in [10,100,1000, 5000, 10**4]:
#     for i in [10, 100, 500, 1000]:
        for i in [10,100,1000]:
            endCountItems = i * NUM_DISTINCT
            print '-----------'
            print 'numItems n=', endCountItems
            if int((EPS**2) * threshold / (C_RAND*NUM_SITES_SYN)) <1:
                print 'bound detCD = ', detBound
            else:
                print 'bound randCD =', int(boundRand)
            l0_sampler = L0ThreshMonitor(stream, EPS, endCountItems, 
                     NUM_SITES_SYN, NUM_DISTINCT, DELTA, threshold)
            l0_sampler.run()
            del l0_sampler
            print 'done.'
            

def test_f0_tracking():
    syn_data = SyntheticData()
    syn_data.loadFromFile(FILE_SHARED_SYN)
    stream = syn_data.get_stream()
    del syn_data
    for i in [10,100,1000,10**4,10**5,10**6, 10**7]:
        endCountItems = i
        print '-----------'
        print 'numItems n=', endCountItems
        f0_track = F0_track(stream, NUM_DISTINCT, 
                 endCountItems, NUM_SITES_SYN, DELTA, eps=0.1)
        f0_track.run()
        f0_track.computeErr()

        
def main():
#     test_ANS_hindex()
#     test_L0()
#     test_ExpCD()
#     test_L0_hindex()
    test_f0_tracking()
    pass


if __name__ == '__main__':
    main()