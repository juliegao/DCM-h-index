'''
Created on Nov 7, 2017

@author: rg522
'''

class F0_threshMonitor(object):
    '''
    classdocs
    '''


    def __init__(self, params):
        '''
        Constructor
        '''


class F0_track():
    def __init__(self, stream, endCountItems):
        self.stream = stream
        self.endCountItems = endCountItems
        self.f0 = 0
    
    def fake(self):
        itemSet = set()
        while self.countItems < self.endCountItems:
            itemID = self.stream[self.countItems][1]
            itemSet.add(itemID)
        self.f0 = len(itemSet)
    
    def getF0(self):
        return self.f0