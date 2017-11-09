'''
@author: Ruoyuan
'''
from algorithms.L0_stream import L0, L0_single_level
from math import log, ceil
from random import randint, shuffle
from algorithms.countDown import CountTracking, RandCD, MultiRandCD
from util.config import SPARSITY, SIZE_WORD, C_RAND
import copy
from algorithms.hash_functions import my_hashxx
from algorithms.f0 import F0_track

class L0_CDM(object):
    '''
    classdocs
    '''
    def __init__(self, stream, eps, endCountItems, 
                 numSites, numDistinct, delta, numSamples=1):
        '''
        @param n: number of distinct items. 
        @param stream: each item is a tuple (siteID, itemID) 
        '''
        self.stream = stream
        self.k = numSites
        self.endCountItems = endCountItems
        self.n = numDistinct
        self.delta = delta
        self.eps = eps
        self.numSamples = numSamples
        self.numRepeats = int(ceil(log(1/delta, 2)))
        self.countMsg = 0
        self.countItems = 0
        self.lv = 0
        self.f0 = 0
        self.init_seeds()
        self.resetSamples(0)
        self.samples_all = [[-1,0] for _ in xrange(self.numSamples)]

        # fake f0 for testing purpose,
        # to remove
        self.itemSet = set()
    
    def resetSamples(self, sampleSize):
        self.samples = [[-1,0] for _ in xrange(sampleSize)]
    
    def update_lv(self, F0):
        '''
        update L0 sample level based on current estimate of F0
        '''
        self.lv = 2 + int(ceil(log(F0, 2)))
    
    def init_seeds(self):
        seeds = range(self.numRepeats * self.numSamples)
        shuffle(seeds)
        seed_i = 0
        self.seed_list = []
        for _ in range(self.numSamples):
            r_list = []
            for _ in range(self.numRepeats):
                r_list.append(seeds[seed_i])
                seed_i += 1
            self.seed_list.append(r_list)
            
    def individual_site_getSample(self, stream, sample_ind, sample_set):
        for seed in self.seed_list[sample_ind]:
            l0_sampler = L0_single_level(
                self.f0, self.lv, SPARSITY, SPARSITY/2, seed)
            for item in stream:
                l0_sampler.sample(item)
            sample = l0_sampler.recover()
            if sample and sample not in sample_set:
                return sample
        return None
    
    def individual_site_getSamples(self, stream, stream_pre):
        return set.difference(set(stream), set(stream_pre))
    
    def individual_site_sample(self, item):
        siteID, itemID = item
        for s in xrange(self.numSamples):
            for r in xrange(self.numRepeats):
                self.sites[siteID][s][r].sample_lv(itemID, self.lv)
        
    def individual_site_recover(self, sample_index, siteID):
        '''
        @return: if not None, return itemID
        '''
        for r in xrange(self.numRepeats):
            sample = self.sites_l0[siteID][sample_index][r].recover_lv(self.lv)
            if sample:
                return sample[0] # sampleID
        return None
    
    def center(self):
        print '>> begin L0 dcm'
        countFail = 0
        lastRoundTracked = False
        log_r = int(ceil(log(self.numRepeats,2)))
        log_n = int(ceil(log(self.n,2)))
        log_r = 1
        log_n = 1
        site_stream = [[] for _ in xrange(self.k)]
        site_stream_preRound = [[] for _ in xrange(self.k)]
        sample_succeed = False
        countF0_change=0  # to delete
        num_samples = 0
        sample_set_all = set()
        
        while self.countItems < self.endCountItems:
            item = self.stream[self.countItems]
            site_stream[item[0]].append(item[1])
            self.countItems += 1
            lastRoundTracked = False
        
            ###########################
            ######## F0 tracking ######
            # pass item to F0
            # use fake f0 tracker
            if self.trackF0(item[1]):
                countF0_change += 1
            ###########################
                if sample_succeed and num_samples == self.numSamples:
                    self.samples = self.trackCount(self.samples)
                    lastRoundTracked = True
#                 elif num_samples == self.numSamples:
#                     self.trackCount(self.samples_all)
                
                if countF0_change %10 == 0:
                    print 'f0=%s' % self.getF0()
                
                self.update_lv(self.f0)
                # update sample upon change of F0
                sites = range(self.k)
                num_samples = min(self.numSamples, self.f0)
                self.resetSamples(num_samples)
                sample_set = set()
                
#                 if countF0_change %10 == 0:
#                     print '>> recover %s samples' %(num_samples)
                    
                if num_samples == self.f0:
                    for site_i in sites:
                        samples = self.individual_site_getSamples(
                            site_stream[site_i], site_stream_preRound[site_i])
                        if not samples:
                            continue
                        for sample in samples:
                            if sample and sample not in sample_set_all:
                                self.countMsg += log_n  # report itemID
                                sample_set_all.add(sample)
                                if len(sample_set_all) >= num_samples:
                                    break
                        if len(sample_set_all) >= num_samples:
                            break
                        
                    count_sample = 0
                    for sample in sample_set_all:
                        self.samples[count_sample][0] = sample
                        self.samples_all[count_sample][0] = sample
                        count_sample += 1
                    sample_set = sample_set_all
                    site_stream_preRound = copy.deepcopy(site_stream)
                else:
                    shuffle(sites)
                    site_samples = [set() for _ in xrange(self.k)]
                    for site_i in sites:
                        self.countMsg += 1  # center asks site to send sample
                        for s in xrange(self.numSamples):
                            sample = self.individual_site_getSample(
                                site_stream[site_i], s, site_samples[site_i])
                            if sample:
                                site_samples[site_i].add(sample)
                            if sample and sample not in sample_set:
                                self.countMsg += log_r  # report hash instance
                                sample_set.add(sample)
                                self.samples[s][0] = sample
                                self.samples_all[s][0] = sample
                                break
                            else:
                                self.countMsg += 1  # report -1 to center
                
                if len(sample_set) < num_samples:
#                     print '>> # of samples recovered: %s' % len(sample_set)
                    countFail += 1
                    sample_succeed = False
                else:
                    if num_samples == self.numSamples:
                        self.countMsg += self.k * log_r # broadcast hash instance
                    sample_succeed = True
#                     if countF0_change %10 == 0:
#                         print '>> success.'
                
        if not lastRoundTracked:
            self.samples = self.trackCount(self.samples_all)
        print '# times failed: %s, fail/sampletimes: %s/%s=%.2f' %(
            countFail, countFail, countF0_change, 1.0*countFail/countF0_change)
               
    def center_tmp(self):
        print '>> begin L0 dcm'
        countFail = 0
        lastRoundTracked = False
        log_r = int(ceil(log(self.numRepeats,2)))
        
        countF0_change=0  # to delete
        
        while self.countItems < self.endCountItems:
            item = self.stream[self.countItems]
            self.countItems += 1
            self.individual_site_sample(item)
            lastRoundTracked = False
            ###########################
            ######## F0 tracking ######
            # pass item to F0
            # use fake f0 tracker
            if self.trackF0(item[1]):
                countF0_change += 1
            ###########################
                self.trackCount()
                lastRoundTracked = True
                
                if countF0_change %100 == 0:
                    print '>> f0 changed to', self.getF0()
                    print '>> track count of last f0'
                
                self.update_lv(self.f0)
                # update sample upon change of F0
                self.resetSamples()
                
                if countF0_change %100 == 0:
                    print '>> recover s samples at each site'
                
                sites = range(self.k)
                sample_set = set()
                for s in xrange(self.numSamples):
                    shuffle(sites)
                    for site_i in sites:
                        sample = self.individual_site_recover(s, site_i)
                        self.countMsg += log_r
                        if sample and sample not in sample_set:
                            sample_set.add(sample)
                            self.samples[s][0] = sample
                            break
                
                if len([sample[0]>-1 for sample in self.samples]) < self.numSamples:
                    countFail += 1
                else:
                    self.countMsg += self.k * log_r # broadcast hash instance
                
        if not lastRoundTracked:
            self.trackCount()
        print '# times failed: %s, fail/F0= %s' %(
            countFail, 1.0*countFail/countF0_change)
        
    
    def trackF0(self, itemID):
        self.itemSet.add(itemID)
        n_distinct = len(self.itemSet)
        change = False
        if n_distinct <= 1:
            change = n_distinct > self.f0
        else:
            change = int(log(n_distinct, 1.1)) > int(log(self.f0, 1.1))
                
        if change:
            self.f0 = n_distinct
            return True
        else:
            return False
        
#         # mimic 1.1 approximate of f0
#         err = uniform(0.9, 1.1)
#         estimate_f0 = n_distinct * err

    def trackCount(self, samples):
#         print 'track %s samples' % sum([i[0]>-1 for i in samples])
        upperBound = 1.0*self.f0 / (1-self.eps)
        for sample in samples:
            sampleID = sample[0]
            if sampleID != -1:
                track = CountTracking(self.stream, self.eps, 
                                      self.countItems, self.k, 
                                      upperBound, self.delta,
                                      [sampleID, sampleID])
                track.run()
                self.countMsg += track.getCountMsg()
                sample[1] = track.getCountEst()
#                 print 'sample count_est, true=', track.getCountEst(), track.getCountTrue()
        return samples
    
    def trackThreshold(self, samples):
        for sample in samples:
            sampleID = sample[0]
            if sampleID != -1:
                track = RandCD(self.stream, self.threshold, 
                               self.eps, self.endCountItems, 
                               self.k, [sampleID, sampleID])
                reachThresh = track.center()
                self.countMsg += track.getCountMsg()
                sample[1] = track.getCountEst()
#                 print 'sample count_est, true=', track.getCountEst(), track.getCountTrue()
        return samples
    
    
    def getSamples(self):
        return self.samples
    
    def getF0(self):
        return self.f0
    
    def getCountMsg(self):
        return self.countMsg    
    
    
class DistinctSampling():
    def __init__(self, stream, eps, threshold, endCountItems, 
                 numSites, numDistinct, delta, numSamples=1):
        '''
        @param n: number of distinct items. 
        @param stream: each item is a tuple (siteID, itemID) 
        '''
        self.stream = stream
        self.threshold = threshold
        self.threshCD = threshold*1.0/(1-eps)
        self.k = numSites
        self.endCountItems = endCountItems
        self.n = numDistinct
        self.eps = eps
        self.numSamples = numSamples
        self.countMsg = 0
        self.countItems = 0
        self.seeds = [randint(0,numDistinct) for _ in xrange(numSamples)]
        self.hash_function = my_hashxx
        self.sample_list = [0]*numSamples
        self.sample_dict = dict()   # sampleID: [reachThresh, CD instance]
        self.globalMax_list = [-1]*numSamples   # global max hash
        self.localMax_list = [[-1]*numSamples for _ in xrange(self.k)]
        self.numInstances = int(ceil(log(1.0*numSamples/(delta*(eps**2)), 2)))
#         boundRand = C_RAND* (1.0/(self.eps**2) - 0.5/self.eps)  
#         boundRand *= self.numInstances
#         print 'bound rand CD =', int(boundRand)
#         print 'bound detCD = ', int(2*self.k*log(1.0/(eps**2),2))
        
        # fake f0 for testing purpose,
        # to remove
        self.f0 = 0
        
    def individual_site(self, item):
        siteID, itemID = item
        sampleID = None
        localMax_list = self.localMax_list[siteID]
        for i in xrange(self.numSamples):
            hash_val = self.hash_function(itemID, self.seeds[i])
            if hash_val > localMax_list[i]:
                self.countMsg += 2    # site and center exchange max item
                localMax_list[i] = hash_val
                if hash_val > self.globalMax_list[i]:
                    self.globalMax_list[i] = hash_val
                    sampleID = itemID
                    
                    self.sample_dict.pop(self.sample_list[i], None)
                    self.sample_list[i] = itemID  # update sample for i-th instance
                else:
                    localMax_list[i] = self.globalMax_list[i]
        return sampleID
        
    def center(self):
        if self.raiseAlert():
            return True
        
        itemSet = set()
        while self.countItems < self.endCountItems:
            item = self.stream[self.countItems]
            itemID = item[1]
            itemSet.add(itemID)
            self.f0 = len(itemSet) 
            self.countItems += 1
            sampleID = self.individual_site(item)
            # if not sampled the first time appears,
            # won't be sampled in the future no matter what instance
            if sampleID:
#                 self.countMsg += (self.k-1) * SIZE_WORD * countInstance
                self.sample_dict[sampleID] = [0,
                                              MultiRandCD(self.stream, self.threshCD, 
                                                self.eps, self.numInstances, self.countItems, 
                                                self.k, [sampleID, sampleID])
                                              ]
            if itemID in self.sample_dict:
                if self.trackThreshold(itemID):
                    return True

        print '>> end of stream' 
        return self.raiseAlert(True)
    

    def trackThreshold(self, itemID):
        reachThresh = self.sample_dict[itemID][0]
        if not reachThresh:
            reachThresh = self.sample_dict[itemID][1].updateEndCountItems(self.countItems)
            if reachThresh:
                self.sample_dict[itemID][0] = 1
                if self.raiseAlert():
                    return True
        return False
    
        
    def getCountMsg(self):
        return self.countMsg
    
    def raiseAlert(self, isLastRound=False):
        countSampleReach = 0
        for reach,_ in self.sample_dict.values():
            countSampleReach += reach
        s = len(self.sample_dict)
        
        if (countSampleReach >= self.threshold or 
            1.0 * countSampleReach * self.f0 / s >= self.threshold):
            
            print '#samples reached threshold, f0: %s/%s, %s' %(
                countSampleReach, s, self.f0)
            
            msgCD = []
            for _,cd in self.sample_dict.values():
                msgCD.append(cd.getCountMsg())
            print 'total cd cost: ', sum(msgCD)
            if sum(msgCD):
                print 'per sample:'
                print msgCD
            print 'msgs before counting:', self.countMsg
            self.countMsg *= SIZE_WORD
            self.countMsg += sum(msgCD)
            return True
        
        if isLastRound:
            print '#samples reached threshold, f0: %s/%s, %s' %(
                countSampleReach, s, self.f0)
            msgCD = []
            for _,cd in self.sample_dict.values():
                msgCD.append(cd.getCountMsg())
            print 'total cd cost: ', sum(msgCD)
            if sum(msgCD):
                print 'per sample:'
                print msgCD
            print 'msgs before counting:', self.countMsg
            self.countMsg *= SIZE_WORD
            self.countMsg += sum(msgCD)
        return False
    
    
    
    
    