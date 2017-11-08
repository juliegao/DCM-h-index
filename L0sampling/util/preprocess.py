'''
Created on Aug 25, 2017
- generate synthetic data for testing
- convert real world data to DCM setting
@author: Julie Gao
'''
from random import randint
from util.config import NUM_ITEMS_SYN, \
    NUM_SITES_SYN, NUM_DISTINCT, FILE_SHARED_SYN, FILE_ANS_SYN, RANGE_OF_ITEM


class SyntheticData():
    '''
    Each item is a tuple (siteID, itemID),
    arranged based on random assignment to each site
    '''
    def __init__(self):
        self.stream = []    # stream: a list of items.
        
    def genCashShared(self, filename=FILE_SHARED_SYN):
        '''
        Items are shared
        - two sites can observe update of the same item
        - the count of that item is the summation of all counts
          
        assuming each update is of count 1,
        and items arrive one by one
        '''
        
        with open(filename,'wb') as fo:
            print ">> saving synthetic data to file:"
            print "  ", filename
            strbuff = ''
            j = 0
            print '>> generating synthetic data: cash register shared'
            for i in xrange(NUM_ITEMS_SYN):
                siteID = randint(0, NUM_SITES_SYN-1)
                itemID = randint(0, NUM_DISTINCT-1)
#                 self.stream.append((siteID, itemID))
                j += 1
                strbuff += '%s %s,'%(siteID, itemID)
                if j %(10**6) == 0:
                    fo.write(strbuff)
                    j = 0
                    strbuff = ''
                    print 'i=',i
            if strbuff: # write remains in buffer
                fo.write(strbuff)
            print "done.\n"
    
    def genAggNonShared(self):
        '''
        Items are not shared
        - two sites cannot observe update of the same item
        - each item is a tuple (siteID, itemCount) for a unique item.
          
        items arrive one by one
        '''
        print '>> generating synthetic data: aggregate nonshared'
        for _ in xrange(NUM_DISTINCT):
            siteID = randint(0, NUM_SITES_SYN-1)
            itemCount = randint(0, RANGE_OF_ITEM-1)
            self.stream.append((siteID, itemCount))
        print "done.\n"
        self.saveToFile(FILE_ANS_SYN)
        
        
    def saveToFile(self):
        pass
    
    
    def loadFromFile(self, filename):
        '''
        @return: a list of items.
        '''
        self.stream = []
        with open(filename, 'rb') as fr:
            print ">> loading synthetic data from file:"
            print "  ", filename
            
            item_list = fr.read().split(',')
            for item in item_list:
                if not item: continue
                siteID, itemID = item.split()
                self.stream.append((int(siteID), int(itemID)))
            print "done.\n"
    
    def get_stream(self):
        return self.stream

if __name__ == '__main__':
    syn_data = SyntheticData()
    syn_data.genCashShared()
#     syn_data.genAggNonShared()
#     syn_data.loadFromFile(FILE_SHARED_SYN)
#     stream = syn_data.get_stream()

#     print "stream sample"
#     for i in xrange(10):
#         print "site %s: %s" %(stream[i][0], stream[i][1])

    pass