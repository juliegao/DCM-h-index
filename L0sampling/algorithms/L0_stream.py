#!/usr/bin/env python
# encoding: utf-8

"""
heavily modified based on
https://github.com/venantius/droplet/blob/master/droplet/samplers.py

Implements L0 sampling on streaming setting
"""

import math

from algorithms.hash_functions import my_hashxx

class OneSparseRecovery():
    def __init__(self):
        self.phi = 0
        self.iota = 0
        self.tau = 0

    def __eq__(self, other):
        if (self.phi == other.phi and \
            self.iota == other.iota and \
            self.tau == other.tau):
            return True
        else:
            return False

    def __hash__(self):
        return hash((self.phi, self.iota, self.tau))

    def update(self, index, weight=1):
        '''
        update finger print upon item update
        '''
        self.phi += weight
        self.iota += weight * index
        self.tau += weight * (index ** 2)

    def is_one_sparse(self):
        return self.iota ** 2 == self.phi * self.tau 


class SSparseRecovery(object):
    """
    An s-sparse recovery comprised of an array of 1-sparse estimators
    """
    def __init__(self, s, num_rows, seed, hash_function):
        '''
        2d array, with log(s/delta_r) rows and 2s columns.
        Each cell contains an instance of 1-sparse recovery
        '''
        self.s = s
        self.num_cols = 2 * s
        self.num_rows = num_rows
        self.array = [[OneSparseRecovery() 
                       for _ in xrange(self.num_cols)] 
                      for _ in xrange(self.num_rows)]
        self.hash_function = hash_function
        self.seed = seed
    
    def is_s_sparse(self):
        '''
        at most s non-zero entries
        '''
        num_zeros = 0
        size = self.num_cols * self.num_rows
        for row in self.array:
            num_zeros += [item.phi for item in row].count(0)
        num_nonzeros = size - num_zeros
        return (num_nonzeros >= 1  and num_nonzeros <= self.s)

    def update(self, i, value):
        """
        update item value at index i
        """
        for row in xrange(self.num_rows):
            seed = str(self.seed) + str(row)
#             seed = row
            col = self.hash_function(i, seed) % self.num_cols
            self.array[row][col].update(i, value)   # update 1-sparse fingerprint

    def recover(self):
        """
        @return v_set: a nonzero vector from the sampled vectors
        """
        v_set = set()
        for row in self.array:
            for item in row:
                if item.phi != 0 and item.is_one_sparse():
                    v_set.add(item)
        return v_set
    
class L0(object):
    """
    A naive implementation of an L0-Sampling data structure, as described in 
    Cormode and Firmani's 2013 paper, "On Unifying the Space of L0-Sampling 
    Algorithms"
    N refers to the n of the input space (e.g. an unsigned 64-bit int in the
    case of most cookie ID spaces)
    k refers to the number of hash functions used in the s-sparse recovery data
    structure.
    s refers to the s of the s-sparse recovery data structure.
    
    In theory, one generally should hold k >= s/2, but in practice C&F 
    note that "it suffices to use small values of k, for instance k=7, to 
    ensure that the failure rate holds steady, independent of the number of 
    repetitions made." Additional notes on this can be found in the 
    accompanying README.
    Also of note: "When time is important, using s<=12 and k<=6 ensures fast 
    computation. On the other hand, by selecting bigger values for both s and 
    k, the process becomes slower than the FIS variant."
    """
    def __init__(self, n, s, k, seed, hash_function = my_hashxx):
        self.n = n
        self.s = s
        self.k = k
        self.hash_function = hash_function
        num_levels = int((math.log(n, 2)))
        self.levels = [SSparseRecovery(s, k, seed, hash_function) 
                       for _ in xrange(num_levels)]
#         # trick to speed up
#         # keep the last levels only
#         self.levels = [0 for _ in xrange(num_levels-1)]
#         self.levels.append(SSparseRecovery(s, k, seed,
#             hash_function))
        self.seed = seed
        self.compute_hash_range()

    def compute_hash_range(self):
        self.min_h=2**32
        self.max_h=0
        for i in range(self.n):
            h = my_hashxx(i, self.seed)
            if h < self.min_h:
                self.min_h = h
            elif h > self.max_h:
                self.max_h = h
    
    def normalize_hash(self, hash_value):
        '''
        normalize hash value to range [0,1]
        '''
        return 1.0*(hash_value - self.min_h) / (self.max_h - self.min_h)
    
    def recover_lv(self, lv):
        """
        recover_lv a nonzero vector from one of the L0 Sampler's levels.
        """
        v_set = set()
        level = self.levels[lv]
        if level.is_s_sparse():
            v_set = level.recover()
        if v_set:
            return self.select(v_set)
        else:
            return None
            
    def recover_greedy(self):
        v_set = set()   
        lv = len(self.levels) - 1
        while lv >= 0:
            level = self.levels[lv]
            if level.is_s_sparse():
                v = level.recover()
                if v:
                    item = self.select_min(v)
                    v_set.add(item[0])
                    self.sample(item[0], -item[1])
                    continue
            lv -= 1
        return v_set
        

    def recover_s(self, sample_size=0):
        """
        recover_lv s=sample_size items
        """
        if not sample_size:
            sample_size = self.n
        samples = []
        v_set = self.recover_greedy()
        while len(samples) < sample_size:
            if v_set:
                item = v_set.pop()
                samples.append(item)
            else:
                break
        return samples

    def select(self, vector_set):
        """
        each item in vector list is a tuple (item, rowID)
        @return: item with hash=0.
        """
        indexes = sorted(vector_set, key=lambda x: 
                self.hash_function(x.iota / x.phi, self.seed))
         
        item = indexes[0]
        i = item.iota / item.phi
        hash_v = self.hash_function(i, self.seed)
#         print 'recovered item=%s, hash=%s, normalize_hash=%s'%(
#             i, hash_v, self.normalize_hash(hash_v))
        if self.normalize_hash(hash_v) == 0:
#             print 'recovered item=%s, hash=%s, normalize_hash=%s'%(
#             i, hash_v, self.normalize_hash(hash_v))
            return (i, item.phi)
        else:
            return None

    def select_min(self, vector_set):
        '''
        @return: one with the lowest hash value.
        '''
        indexes = sorted(vector_set, key=lambda x: 
                self.hash_function(x.iota / x.phi, self.seed))
        item = indexes[0]
        i = item.iota / item.phi
        return (i, item.phi)
        
    def sample(self, i, value=1):
        '''
        if n^3 2^(-j) >= h(i), then a(j)_i = a_i
        by default, each update is +1 
        '''
        hash_v = self.hash_function(i, self.seed)
        for j, level in enumerate(self.levels):
#             if (self.n**3) * (2**-(j+1)) >= hash_v:
#             if (self.max_h - self.min_h) * (2.0 ** -(j+1)) >= hash_v:
            if 2.0**-(j+1) >= self.normalize_hash(hash_v):
                level.update(i, value)

        # trick to speed up:
        # knowing logn,
        # used only second to the last level
        
#             print 'item=%s, hash=%s' % (i, hash_v)
#             print 'hash % n=', hash_v % self.n
        

#         if (self.n * (2 ** -(len(self.levels))) # >=1
#             >= hash_v % self.n + 1):
#         if (self.max_h - self.min_h) * (2.0 ** -(len(self.levels))) >= hash_v:
#             self.levels[-1].update(i, value)
#             print 'hash_v=%s,n=%s'% (hash_v, self.n)
#             print 'n*s^-(lv+1)=%s, hash_v%%n+1=%s' %(
#             (self.n * (2 ** -(len(self.levels)))), (hash_v % self.n + 1))

                
class L0_single_level():
    def __init__(self, n, top_level, s, k, seed, hash_function = my_hashxx):
        self.n = n
        self.hash_function = hash_function
        self.level = SSparseRecovery(s, k, seed, hash_function)
        self.seed = seed
        self.top_level = top_level
        self.compute_hash_range()

    def compute_hash_range(self):
        self.min_h=2**32
        self.max_h=0
        for i in range(self.n):
            h = my_hashxx(i, self.seed)
            if h < self.min_h:
                self.min_h = h
            elif h > self.max_h:
                self.max_h = h
                
    def normalize_hash(self, hash_value):
        '''
        normalize hash value to range [0,1]
        '''
        return 1.0*(hash_value - self.min_h) / (self.max_h - self.min_h)
    
    def sample(self, i, value=1):
        '''
        @param lv: sample at levels >= lv   
        if R 2^(-j) >= h(i), then a(j)_i = a_i
        by default, each update is +1
        '''
        hash_v = self.hash_function(i, self.seed)
        if 2.0**-(self.top_level+1) >= self.normalize_hash(hash_v):
            self.level.update(i, value)
            
    def recover(self):
        """
        recover_lv a nonzero vector from one of the L0 Sampler's levels.
        """
        if self.level.is_s_sparse():
            v_set = self.level.recover()
            if v_set:
                return self.select(v_set)
        return None
        
        
    def select(self, v_set):
        indexes = sorted(v_set, key=lambda x: 
                self.hash_function(x.iota / x.phi, self.seed))
        item = indexes[0]
        i = item.iota / item.phi
        hash_v = self.hash_function(i, self.seed)
        if self.normalize_hash(hash_v) == 0:
            return i
        else:
            return None
        
    def recover_all(self):
        if self.level.is_s_sparse():
            v_set = self.level.recover()
            if v_set:
                return [v.iota/v.phi for v in v_set]
        return None