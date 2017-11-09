'''
https://github.com/venantius/droplet/blob/master/droplet/hash_functions.py
'''
from util.config import NUM_DISTINCT, EPS

#!/usr/bin/env python
# encoding: utf-8

"""
A set of hash functions known to have good uniformity properties
"""
PRIME = 2147483647  # large prime 2^31-1
NUM_BINS = 6* ((96.0/(EPS**2))**2)

from pyhashxx import hashxx
# from algorithms.pymmh3 import hash


def my_hashxx(x, seed=0):
    x = str(x)
    seed = int(seed)
    return hashxx(x, seed = seed)


def pairwiseHashF(x, seed):
    '''
    h(x) = (ax + b) % p % n
    p is 2^31-1
    @param seed: [a,b] 
    '''
    a,b = seed
    result = 1.0*(a*x + b)
    result = result % PRIME
    result = int(result % NUM_DISTINCT)
    return result
    
def pairwiseHashG(x, seed):
    '''
    h(x) = (ax + b) % p % [6*(96/e^2)^2]
    p is 2^31-1
    @param seed: [a,b] 
    '''
    a,b = seed
    result = 1.0*(a*x + b)
    result = result % PRIME
    result = int(result % NUM_BINS)
    return result
        
    
    
    