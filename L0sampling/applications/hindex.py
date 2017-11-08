'''
Created on Aug 25, 2017

1. To have randCD better than DetCD, let
    1/eps^2 log( log(n, 2) / (delta * eps))/2
    only need half raise alert
    < k log 1/eps
2. Only use detCD when eps^2 t/(ck) > 1

with eps = 0.2, delta = 0.1,
then k > 62.5, 
set k =100, then t > 240,000

@author: Ruoyuan Gao
'''
from algorithms.countDown import RandCD, CountTracking, DetCD, getThresholdList
from math import log, ceil
from algorithms.L0_CDM import L0_CDM, DistinctSampling
from algorithms.f0 import F0_track

class Hindex(object):
    '''
    compute h-index:
    1. exponential histogram
    2. l0 sampling
    '''

    def __init__(self, stream, eps, endCountItems, 
                 numSites, numDistinct, delta):
        self.stream = stream
        self.eps = eps
        self.k = numSites
        self.n = numDistinct
        self.delta = delta
        self.endCountItems = endCountItems
        if self.endCountItems > len(stream):
            self.endCountItems = len(stream)
        
        self.countItems = 0
        self.countMsg = 0
        self.h_est = 0
        self.h_true = 0            

    def getCountMsg(self):
        return self.countMsg
    
    def computeTrueHindex_ANS(self):
        N = self.n

        freqCountList = [0] * (N+1)  # max of h is N
        for i in xrange(self.endCountItems):
            freq = self.stream[i][1]
            if freq >= N:
                freqCountList[N] += 1
            else:
                freqCountList[freq] += 1
        
        if freqCountList[N] >= N:
            self.h_true = N
        else:
            N = N - 1
            while N >= 0:
                freqCountList[N] += freqCountList[N+1]
                if freqCountList[N] >= N:
                    self.h_true = N
                    return
                N -= 1
                
    def computeTrueHindex_CS(self):
        print '>> compute True h-index'
        N = self.n
        
        itemFreqList = [0] * (N+1)
        for i in xrange(self.endCountItems):
            itemID = self.stream[i][1]
            itemFreqList[itemID] += 1
         
        freqCountList = [0] * (N+1)  # max of h is N
        for freq in itemFreqList:
            if freq >= N:
                freqCountList[N] += 1
            else:
                freqCountList[freq] += 1
        del itemFreqList
        
        if freqCountList[N] >= N:
            self.h_true = N
        else:
            N = N - 1
            while N >= 0:
                freqCountList[N] += freqCountList[N+1]
                if freqCountList[N] >= N:
                    self.h_true = N
                    return
                N -= 1
    
    def setThresholdList(self, windowSize):
        self.thresholdList = getThresholdList(self.n, windowSize)
        
    def computeErr(self):
        '''
        if estimate is within range:
        (1-eps)h_true <= h_est <= (1+eps)h_true
        '''
        if not self.h_true:
            print 'h=%s,h*=0' % self.h_est
            return
        
        ratio = 1.0 * self.h_est / self.h_true
        left = 1 - self.eps
        right = 1 + self.eps
        print 'h/h*, [1-eps,1+eps], in range?: %s/%s=%.2f, [%s,%s], %s' %(
            self.h_est, self.h_true, ratio,
            left, right, left <= ratio and ratio <= right)        
    
    
class ExpHist(Hindex):
    '''
    Exponential histogram with randomized countdown
    '''
    def __init__(self, stream, eps, endCountItems, 
                 numSites, numDistinct, delta):
        Hindex.__init__(self, stream, eps, endCountItems, 
                 numSites, numDistinct, delta)
        self.epsCD = eps / (1+eps) 
        self.epsHist = eps / (1-eps)
        self.setThresholdList(1+self.epsHist)
   
    def run(self):
        deltaCDinverse = log(self.n, 2) / (self.delta * self.epsHist)
        numCDs = ceil(log(deltaCDinverse, 2))
        if ((1/(self.epsCD**2) + self.k) * numCDs 
            < self.k * log(1/(self.epsCD**2),2)):
            self.useRandCD(numCDs)
        else:
            self.useDetCD()
        self.computeTrueHindex_ANS()
        
    def useRandCD(self, numCDs):
        i = 0
        while i < len(self.thresholdList):
            t = self.thresholdList[i]
            if t > self.endCountItems:
                return
            costCD = ( self.k * 2
                       + ( 1/(self.epsCD**2) + self.k ) 
                       * log(1.0*t/self.k , 1+self.epsHist)
                       * numCDs)
            if costCD >= t:
                self.reportAll(t)
                i += 1
            else: break
        
        print 't=%s, use rand' %self.thresholdList[i]
        while i < len(self.thresholdList):
            t = self.thresholdList[i]
            if t > self.endCountItems:
                return
        
            countTrue = 0
            self.countItems = 0
            for j in xrange(numCDs):
                randCD = RandCD(self.stream, t, self.epsCD, 
                            self.endCountItems, self.k,
                            [t])
                if randCD.center():
                    countTrue += 1
                self.countMsg += randCD.getCountMsg()
                if randCD.getCountItems() > self.countItems:
                    self.countItems = randCD.getCountItems()
                if randCD.get_useDetCD():
                    # used detCD
                    # no need to repeat
                    if randCD.center():
                        self.h_est = t
                    break
                if countTrue >= ceil(numCDs / 2):    
                    # half raises alert
                    self.h_est = t
                    break
                elif j+1 - countTrue > ceil(numCDs / 2):
                    # half failed
                    break
    #             print 'h=%s\nitems, msg, local: %s %s %s' %(
    #                 h, countItems, self.countMsg, localMsg)
            i += 1
    
    def useDetCD(self):
        i = 0
        while i < len(self.thresholdList):
            t = self.thresholdList[i]
            if t > self.endCountItems:
                return
            costCD = ( self.k * 2
                       + 2*self.k * log(1/self.epsCD,2)
                       * log(1.0*t/self.k, 1+self.epsHist))
            if t <= self.k or costCD >= t:
                self.reportAll(t)
                i += 1
            else: break
        
#         print 't=%s, use det' %self.thresholdList[i]
        while i < len(self.thresholdList):
            t = self.thresholdList[i]
            if t > self.endCountItems:
                return
            detCD = DetCD(self.stream, t, self.epsCD, 
                        self.endCountItems, self.k,
                        [t])
            if detCD.centerApprox():
                self.h_est = t
            self.countMsg += detCD.getCountMsg()
            if detCD.getCountItems() > self.countItems:
                self.countItems = detCD.getCountItems()
            i += 1

        
    def reportAll(self, h):
        while self.countItems < self.endCountItems:
            self.countItems += 1
            self.countMsg += 1
            self.computeTrueHindex_ANS()
            if self.h_true >= h:
                return
    
    def getCountItems(self):
        return self.countItems
        
class L0sampling(Hindex):
    '''
    L0 with randomized countdown
    '''
    def __init__(self, stream, eps, endCountItems, 
                 numSites, numDistinct, delta):
        Hindex.__init__(self, stream, eps, endCountItems, 
                 numSites, numDistinct, delta)
        self.epsHist = eps / (1-eps)
        self.sample_size = int(ceil((3*log(2/delta, 2) / (eps**2))))
        # see "streaming algorithms for measuring h-impact threorem 14"
    
    def run(self):
        l0 = L0_CDM(self.stream, self.eps, self.endCountItems, 
                    self.k, self.n, self.delta, self.sample_size)
        l0.center()
        samples = l0.getSamples()
#         print 'final samples:', samples
        self.countMsg = l0.getCountMsg()
        thresholdList = getThresholdList(self.n, 1+self.epsHist)
        for t in thresholdList:
            R = sum([s[1]>=t for s in samples])
            r = 1.0 * R * l0.getF0() / self.sample_size
            if r >= t * (1-self.epsHist):
                self.h_est = t
            else:
                break
        del thresholdList
        self.computeTrueHindex_CS()
        
class L0ThreshMonitor(Hindex):
    def __init__(self, stream, eps, endCountItems, 
                 numSites, numDistinct, delta, threshold):
        Hindex.__init__(self, stream, eps, endCountItems, 
                 numSites, numDistinct, delta)
        self.sample_size = int(ceil((3*log(2/delta, 2) / (eps**2))))
        self.threshold = threshold
        # see "streaming algorithms for measuring h-impact threorem 14"
    
    def run(self):
        print '>> start L0 with threshold counting, sample size: %s' %self.sample_size
        l0_dcm = DistinctSampling(self.stream, self.eps, 
                                  self.threshold, self.endCountItems, 
                                  self.k, self.n, self.delta, 
                                  self.sample_size)
        reachThresh = l0_dcm.center()
        countMsg = l0_dcm.getCountMsg()
#         samples = l0_dcm.getSamples()
        # f0 tracking algorithm
#         f0_track = F0_track(self.stream, self.endCountItems)
#         f0 = f0_track.getF0()
#         reachThresh = self.ifReachThresh(samples, f0)
        self.computeTrueHindex_CS()
        print 'h_true, threshold, h_est over thresh?: %s, %s, %s' %(
            self.h_true, self.threshold, reachThresh)
        print 'total bits excluding f0: %s' %countMsg
        
    
#     def ifReachThresh(self, samples, f0):
#         s = len(samples)
#         countThresh_sample = sum(samples)
#         countThresh_all = 1.0 * countThresh_sample * f0 / s
#         return countThresh_all >= self.threshold
#         
        
        
        
    
