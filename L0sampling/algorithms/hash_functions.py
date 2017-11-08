'''
https://github.com/venantius/droplet/blob/master/droplet/hash_functions.py
'''

#!/usr/bin/env python
# encoding: utf-8

"""
A set of hash functions known to have good uniformity properties
"""

from pyhashxx import hashxx
# from algorithms.pymmh3 import hash

# def murmurhash3_32(item, seed = 0):
#     if type(item) is not str: 
#         item = str(item)
#     if type(seed) is not int:
#         seed = int(seed)
#     return hashxx(item, seed = seed)

def my_hashxx(item, seed=0):
    item = str(item)
    seed = int(seed)
    return hashxx(item, seed = seed)

