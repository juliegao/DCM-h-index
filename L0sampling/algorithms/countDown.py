'''
Created on Mar 15, 2017

@author: Julie Gao

class CountDown: base class for several count down algorithms
class DetCD(CountDown): 
    Deterministic countdown
    Exact: O(k log T/k)
    Approximate: O(k log 1/eps)
class RandCD(CountDown): 
    Randomized countdown: O(1/eps^2)
'''
from random import randint
from math import ceil, log
from util.config import C_RAND

class CountDown(object):
    '''
    base class for several count down algorithms
    '''
    def __init__(self, stream, threshold, eps,
                 endCountItems, numSites, countCondition=[]):
        '''
        :param a stream assigned to diff sites.
            Each item is a tuple (siteID, itemID)
        :param threshold: stop when value_est reaches threshold T
        :param countCondition: only count items that are in the range
            e.g., [a,b]: count items a<= x <=b
            e.g., [a]: count items a<= x
            e.g., []: count any items
        '''
        self.stream = stream
        self.threshold = threshold
        self.eps = eps
        self.k = numSites
        self.countCondition = countCondition
        self.endCountItems = endCountItems
        if self.endCountItems > len(stream):
            self.endCountItems = len(stream)
         
        self.countMsg = 0
        self.countItems = 0
        self.sites = [0] * numSites
        
    def meetCountCondition(self, itemID):
        if not self.countCondition:
            return True
        if len(self.countCondition) == 2:
            return (itemID >= self.countCondition[0] 
                    and itemID <= self.countCondition[1])
        else:
            return itemID >= self.countCondition[0]
    
    def getCountMsg(self):
        # Question: include k or not?
#         return self.countMsg + self.k   # broadcast to stop
        return self.countMsg
    
    def getCountItems(self):
        return self.countItems
    
    def updateEndCountItems(self, endCountItems):
        self.endCountItems = endCountItems

 
class DetCD(CountDown):
    '''
    Deterministic countdown
    Exact: O(k log T/k)
    Approximate: O(k log 1/eps)
    '''
    def __init__(self, stream, threshold, eps,
                 endCountItems, numSites, countCondition=[]):
        CountDown.__init__(self,stream, threshold, eps,
                 endCountItems, numSites, countCondition)
        self.roundNum = 0
        self.localThresh = 0
        self.currCount = 0
        self.nextRound()
        
    def nextRound(self, tag=''):
        '''
        Start next round.
        Computes the local threshold for current round:
        T/(2^roundNum * k)
        '''
        if self.roundNum > 0:
            self.countMsg += self.k # broadcast to next round
        self.numReportsCurrRound = 0
        self.roundNum += 1
        self.localThresh = int( 
            self.threshold*1.0 / (2**(self.roundNum) * self.k)
            )
        if self.localThresh < 1:
            self.localThresh = 1
            
        for _ in xrange(2):
            # check remainders of each site
            # at most 2 iters will suffice
            for siteID in xrange(self.k):
                signal = self.individual_site([siteID, -1])
                if signal:
                    self.onSignal(tag)
        
    def raiseAlertExact(self):
        if self.currCount >= self.threshold:
#             self.countMsg += self.k # broadcast to end
            return True
        else:
            return False
    
    def raiseAlertApprox(self):
        boundUnreported = int(self.eps * self.threshold)
        if self.threshold - self.currCount <= boundUnreported:
#             self.countMsg += self.k # broadcast to end
            return True
        else:
            return False
        
    def individual_site(self, item):
        '''
        @return: True if reached local threshold
        '''
        siteID, itemID = item
        if itemID != -1:    # normal receiving items
            if not self.meetCountCondition(itemID):
                return False
            self.sites[siteID] += 1 # increment local count
            
        # if itemID =-1
        # check remainder directly
        if self.sites[siteID] >= self.localThresh:
            self.sites[siteID] -= self.localThresh
            return True
        else:
            return False

    
    def centerExact(self):
        if self.raiseAlertExact():
            return True
        while self.countItems < self.endCountItems:
            item = self.stream[self.countItems]
            self.countItems += 1
            signal = self.individual_site(item)
            if signal and self.onSignal('exact'):
                return True            
        return False

    def centerApprox(self):
        if self.raiseAlertApprox():
            return True
        while self.countItems < self.endCountItems:
            item = self.stream[self.countItems]
            self.countItems += 1
            signal = self.individual_site(item)
            if signal and self.onSignal():
                return True
        return False

    def onSignal(self, tag=''):
        self.countMsg += 1
        self.numReportsCurrRound += 1
        self.currCount += self.localThresh
        if tag == 'exact':
            if self.raiseAlertExact():
                return True
        elif self.raiseAlertApprox():
            return True
        if (self.localThresh > 1 and 
            self.numReportsCurrRound >= self.k):
            self.nextRound(tag)
            
class RandCD(CountDown):
    '''
    Randomized countdown: O(1/eps^2)
    set constant c=96
    '''
    def __init__(self, stream, threshold, eps,
                 endCountItems, numSites, countCondition=[]):
        CountDown.__init__(self,stream, threshold, eps,
                 endCountItems, numSites, countCondition)
        self.countReports = 0
        c = C_RAND   # or 48?
        self.localThresh = int(
            (self.eps**2) * self.threshold / (c*self.k)
            )   # e^2*T/(ck)
        self.reportsRequired = int(
            ceil(c/self.eps**2 - c/(2*self.eps))
            )   # c/e^2 - c/(2e)
        
        self.detCD = None
        
    def raiseAlert(self):
        if self.countMsg >= self.reportsRequired:
#             self.countMsg += self.k # broadcast to end
            return True
        else:
            return False
    
    def individual_site(self, item):
        siteID, itemID = item
        if self.meetCountCondition(itemID):
            self.sites[siteID] += 1 # increment local count
            if self.sites[siteID] >= self.localThresh:
                self.sites[siteID] = 0  # reset local count
                return randint(1,self.k) == 1  # p=1/k
        return False
            
    def center(self):
        if self.localThresh < 1:    # use detCD
            self.useDetCD = True
            if not self.detCD:  # new tracking
                self.detCD = DetCD(self.stream, self.threshold, self.eps,
                     self.endCountItems, self.k, self.countCondition)
            else:   # continue tracking
                self.detCD.updateEndCountItems(self.endCountItems)
            alert = self.detCD.centerApprox()
            self.countMsg = self.detCD.getCountMsg()
            self.countItems = self.detCD.getCountItems()
            return alert
        
        # else: use randomized
        self.useDetCD = False   
        if self.raiseAlert():
            return True 
        while self.countItems < self.endCountItems:
            item = self.stream[self.countItems]
            self.countItems += 1
            signal = self.individual_site(item)
            if signal:
                self.countMsg += 1
                if self.raiseAlert():
                    return True
        return False
    
    def get_useDetCD(self):
        return self.useDetCD
    

class MultiRandCD():
    def __init__(self, stream, threshold, eps, numCDs, 
                 endCountItems, numSites, countCondition):
        self.stream = stream
        self.threshold = threshold
        self.eps = eps
        self.endCountItems = endCountItems
        self.k = numSites
        self.countCondition = countCondition
        self.numInstances = numCDs
        self.countMsg = 0
        self.instance_list = []
        
    def run(self):
        countTrue = 0
        self.countMsg = 0
        if not self.instance_list:
            for _ in xrange(self.numInstances):
                self.instance_list.append(RandCD(self.stream, self.threshold, 
                               self.eps, self.endCountItems, 
                               self.k, self.countCondition))
        else:
            for i in xrange(self.numInstances):
                self.instance_list[i].updateEndCountItems(self.endCountItems)
                
        for i in xrange(self.numInstances):
            reachThresh = self.instance_list[i].center()
            self.countMsg += self.instance_list[i].getCountMsg()
            if self.instance_list[i].get_useDetCD():
                return reachThresh
            if reachThresh:
                countTrue += 1             
    
        halfCDs = ceil(self.numInstances/2.0)             
        if countTrue >= halfCDs:
            return True
        else:
            return False
        
    
    def getCountMsg(self):
        return self.countMsg
    
    def updateEndCountItems(self, endCountItems):
        self.endCountItems = endCountItems
        return self.run()
    
class CountTracking(CountDown):
    def __init__(self, stream, eps, endCountItems, 
                 numSites, numDistinct, delta,
                 countCondition):
        CountDown.__init__(self, stream, 0, eps,
                 endCountItems, numSites, countCondition)
        self.n = numDistinct
        self.delta = delta
        
        self.value_est = 0  # count to be tracked
        self.value_true = 0
        self.epsCD = eps / (1+eps) 
        self.epsHist = eps / (1-eps)
        self.setThresholdList(1+self.epsHist)    
        del self.sites, self.threshold
   
    def run(self):
        deltaCDinverse = log(self.n, 2) / (self.delta * self.epsHist)
        numCDs = int(ceil(log(deltaCDinverse, 2)))
        if ((1/(self.epsCD**2) + self.k) * numCDs 
            < self.k * log(1.0/(self.epsCD**2),2)):
            self.useRandCD(numCDs)
        else:
            self.useDetCD()
        
        
    def useRandCD(self, numCDs):
        i = 0
        while i < len(self.thresholdList):
            t = self.thresholdList[i]
            if t > self.endCountItems:
                return
            costCD = ( self.k * 2
                       + ( 1.0/(self.epsCD**2) + self.k ) 
                       * log(1.0*t/self.k , 1+self.epsHist)
                       * numCDs)
            if costCD >= t:
                self.reportAll(t)
                i += 1
            else: break
        
#         print 't=%s, use rand' %self.thresholdList[i]
        while i < len(self.thresholdList):
            t = self.thresholdList[i]
            if t > self.endCountItems:
                return
        
            localMsg = 0    # to delete
            countTrue = 0
            self.countItems = 0
            for j in xrange(numCDs):
                randCD = RandCD(self.stream, t, self.epsCD, 
                            self.endCountItems, self.k,
                            self.countCondition)
                if randCD.center():
                    countTrue += 1
                self.countMsg += randCD.getCountMsg()
                localMsg += randCD.getCountMsg()
                if randCD.getCountItems() > self.countItems:
                    self.countItems = randCD.getCountItems()
                if randCD.get_useDetCD():
                    # used detCD
                    # no need to repeat
                    if randCD.center():
                        self.value_est = t
                    break
                if countTrue >= ceil(numCDs / 2.0):    
                    # half raises alert
                    self.value_est = t
                    break
                elif j+1 - countTrue > ceil(numCDs / 2.0):
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
                        self.countCondition)
            if detCD.centerApprox():
                self.value_est = t
            self.countMsg += detCD.getCountMsg()
            if detCD.getCountItems() > self.countItems:
                self.countItems = detCD.getCountItems()
            i += 1
    
    def getCountEst(self):
        return self.value_est
    
    def getCountTrue(self):
        count_true = 0
        for i in range(self.endCountItems):
            itemID = self.stream[i][1]
            if self.meetCountCondition(itemID):
                count_true += 1
        return count_true
    
    def computeErr(self):
        '''
        1-eps <= count_estimate / count_true <= 1+eps
        '''
        self.value_true = self.getCountTrue()
        if not self.value_true:
            print 'est=%s, true=0' % self.value_est
        
        ratio = self.value_est * 1.0 / self.value_true
        left = 1 - self.eps
        right = 1 + self.eps
        print 'h/h*, [1-eps,1+eps], in range?: %s/%s=%.2f, [%s,%s], %s' %(
            self.value_est, self.value_true, ratio,
            left, right, left <= ratio and ratio <= right) 
      
    def setThresholdList(self, windowSize):
        self.thresholdList = getThresholdList(self.n, windowSize)
    
    def reportAll(self, t):
        while self.countItems < self.endCountItems:
            itemID = self.stream[self.countItems][1]
            self.countItems += 1
            if self.meetCountCondition(itemID):
                self.value_est += 1
                self.countMsg += 1
                if self.value_est >= t:
                    return

def getThresholdList(n, windowSize):
    '''
    @param n: upper bound threshold
    '''
    thresholdList = [1]
    i = 1
    while True:
        threshold = int(windowSize ** i)
        if threshold > n:
            break
        elif threshold > thresholdList[-1]:
            thresholdList.append(threshold)
        i += 1
    if n not in thresholdList:
        thresholdList.append(n)
    return thresholdList