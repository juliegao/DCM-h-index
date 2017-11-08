'''
Created on Aug 25, 2017

@author: Julie Gao
'''
import sys
import os
from math import ceil, log

#===============================================================================
# parameters
#===============================================================================
EPS = 0.3
C_RAND = 48.0   # constant in randomized CD bundle
EPS_L0 = 0.02
                # pairwise hash s=2/eps=100
DELTA = 0.2     # repetitions 1/eps*log(1/delta)=300

# to make e^2 T /ck > 1 and use randCD,
# e=0.3, s=1/e^2=12, delta=0.2, repetition=log(s/delta), k>s= 15
# then T> 48*k/e^2=8000

#------------------------------------------------------------------------------ 
# synthetic data
#------------------------------------------------------------------------------ 
# 1. For aggregate non shared (ANS)
# NUM_SITES_SYN = 100
# NUM_DISTINCT  = 10**6
# NUM_ITEMS_SYN = 2*NUM_DISTINCT
# RANGE_OF_ITEM = 2*NUM_DISTINCT

# 2. For cash register shared
NUM_SITES_SYN = 15
# NUM_ITEMS_SYN = (10**5)*5   # for testing
# NUM_DISTINCT = 10**5    # for testing
# NUM_ITEMS_SYN = 10**8   # for testing
NUM_DISTINCT = 10**4    # for testing
NUM_ITEMS_SYN = 10**7   # for testing
RANGE_OF_ITEM = NUM_DISTINCT

SPARSITY = 24
SIZE_WORD = int(ceil(log(NUM_DISTINCT,2)))

# NUM_ITEMS_SYN = 5*(10**7)
# NUM_DISTINCT = 5*(10**6) #inverse distri paper

#===============================================================================
# directories
#===============================================================================
# HOME_DIR = "C:/Users/rg522/workspace/research/L0sampling/"
HOME_DIR = "/grad/users/rg522/workspace/L0sampling"
DATA_DIR = os.path.join(HOME_DIR, 'data')
SYNTHETIC_DIR = os.path.join(DATA_DIR, 'synthetic')
REAL_DIR = os.path.join(DATA_DIR, 'real')

# cash register model, shared, file format:
# siteID1 itemID1,siteID2 itemID2
FILE_SHARED_SYN = os.path.join(
    SYNTHETIC_DIR, "kNn_%s.%s.%s" % (NUM_SITES_SYN,NUM_DISTINCT,NUM_ITEMS_SYN))
FILE_ANS_SYN = os.path.join(
    SYNTHETIC_DIR, "ANS_kNr_%s.%s.%s" % (NUM_SITES_SYN,NUM_DISTINCT,RANGE_OF_ITEM))

