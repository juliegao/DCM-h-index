'''
Created on Nov 7, 2017

@author: rg522
'''
from algorithms.hash_functions import pairwiseHashF, pairwiseHashG, NUM_BINS,\
    PRIME
from math import ceil, log, floor
from random import randint
from algorithms.countDown import getThresholdList

class F0_threshMonitor():
    def __init__(self, stream, threshold, eps, endCountItems, numSites, seed):
        self.stream = stream
        self.endCountItems = endCountItems
        self.seed = seed
        self.eps = eps
        self.k = numSites
        self.threshold = threshold
        self.hash_function_f = pairwiseHashF
        self.hash_function_g = pairwiseHashG
        self.local_buff = [set() for _ in xrange(self.k)]
        self.global_buff = set()
        self.countMsg = 0
        self.msg_size = int(ceil(log(NUM_BINS,2)))
        self.countItems = 0
        
        self.init_params()
        
    def init_params(self):
        '''
        48/Te^2 <= T/2^t < 96/Te^2
        '''
        low = (self.eps**2) * self.threshold / 48.0
#         high = (self.eps**2) * self.threshold / 96.0
        low = ceil(log(low, 2))
#         high = floor(log(high, 2))
        self.t = int(low)
        self.reportsRequired = int(ceil(
            (1 - self.eps/2) * self.threshold / (2**self.t)
                                ))
        
    def checkLastBits(self, val):
        mask = (1<<self.t) - 1
        if val & mask:
            return True
        return False
        
    def individual_site(self, item):
        siteID, itemID = item
        site_buff = self.local_buff[siteID]
        hash_f = self.hash_function_f(itemID, self.seed)
        if self.checkLastBits(hash_f):
            hash_g = self.hash_function_g(itemID)
            if hash_g not in site_buff:
                self.countMsg += 1
                site_buff.add(hash_g)
                return hash_g
        return None
        
    def center(self):
        if self.raiseAlert():
            return True
        while self.countItems < self.endCountItems:
            item = self.stream[self.countItems]
            self.countItems += 1
            hash_g = self.individual_site(item)
            if hash_g and hash_g not in self.global_buff:
                self.global_buff.add(hash_g)
                if self.raiseAlert():
                    return True
        return False
                
    def raiseAlert(self):
        return (len(self.global_buff) >= self.reportsRequired)
                
    def getCountMsg(self):   
        return self.countMsg * self.msg_size    
        
    def updateEndCountItems(self, endCountItems):
        self.endCountItems = endCountItems
      
        
class MultiF0():
    def __init__(self, stream, threshold, eps, 
                 endCountItems, numSites, numInstances):
        self.stream = stream
        self.threshold = threshold
        self.eps = eps
        self.endCountItems = endCountItems
        self.k = numSites
        self.numInstances = numInstances
        self.instance_list = []  # f0 threshold monitor instances
        
        self.init_seeds()
        
    def init_seeds(self):
        self.seed_list = []
        countSeeds = 0
        while countSeeds < self.numInstances:
            a = randint(1,PRIME-1)
            b = randint(1,PRIME-1)
            seed = [a,b]
            if seed not in self.seed_list:
                self.seed_list.append(seed)
                countSeeds += 1
        
    def run(self):
        countTrue = 0
        self.countMsg = 0
        if not self.instance_list:
            for i in xrange(self.numInstances):
                self.instance_list.append(
                    F0_threshMonitor(self.stream, self.threshold, 
                                     self.eps, self.endCountItems,
                                     self.k, self.seed_list[i]))
        else:
            for i in xrange(self.numInstances):
                self.instance_list[i].updateEndCountItems(self.endCountItems)
                
        for i in xrange(self.numInstances):
            reachThresh = self.instance_list[i].center()
            self.countMsg += self.instance_list[i].getCountMsg()
            if reachThresh:
                countTrue += 1
                
        halfIns = ceil(self.numInstances/2.0)                  
        if countTrue >= halfIns:
            return True
        else:
            return False
        
    def updateEndCountItems(self, endCountItems):
        self.endCountItems = endCountItems
        return self.run()
    
        
class F0_track():
    def __init__(self, stream, numDistinct, 
                 endCountItems, numSites, delta, eps=0.1):
        self.stream = stream
        self.n = numDistinct
        self.endCountItems = endCountItems
        self.k = numSites
        self.eps = eps
        
        self.f0_true = 0
        self.f0_est = 0
        self.epsTM = eps / (1+eps)  # threshold monitor
        self.epsHist = eps / (1-eps)
        self.setThresholdList(1+self.epsHist)
        
        deltaTMinverse = log(self.n, 2) / (delta * self.epsHist)
        self.numInstances = int(ceil(log(deltaTMinverse, 2)))
        self.itemSet = set()
        self.multiF0_list = []
        
    def run(self):
        self.countMsg = 0
        i = 0
        if not self.multiF0_list:
            while i < len(self.thresholdList):
                t = self.thresholdList[i]
                if t > self.endCountItems:
                    return
                self.multiF0_list.append(
                    MultiF0(self.stream, t, self.epsTM, 
                            self.endCountItems, self.k, self.numInstances)
                                         )
                if self.multiF0_list[i].run():
                    self.f0_est = t
                self.countMsg += self.multiF0_list[i].getCountMsg()
                i += 1
        else:
            while i < len(self.thresholdList):
                t = self.thresholdList[i]
                if t > self.endCountItems:
                    return
                reachThresh = self.multiF0_list[i].updateEndCountItems(self.endCountItems)
                if reachThresh:
                    self.f0_est = t
                
        
    def setThresholdList(self, windowSize):
        self.thresholdList = getThresholdList(self.n, windowSize)
    
    def getF0(self):
        return self.f0_est
    
    def computeErr(self):
        self.f0_true = len(self.itemSet)
        err = abs(self.f0_est - self.f0_true) * 1.0 / self.f0_true
        print 'f0_true, f0_est, relative err: %s, %s, %s' %(
            self.f0_true, self.f0_est, err)
        print ''