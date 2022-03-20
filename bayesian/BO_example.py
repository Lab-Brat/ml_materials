import numpy as np
#import cPickle as pickle
import scipy
import combo
import os
import urllib
import time

start = time.time()

def load_data():
    A =  np.asarray(np.loadtxt('SiSi_14.data'))
    print (A.shape)
    X = A[:,0:14]
    t  = A[:,14]
    return X, t

log=open('bo_log_02','w')

# Load the data. 
# X is the N x d dimensional matrix. Each row of X denotes the d-dimensional feature vector of search candidate. 
# t is the N-dimensional vector that represents the corresponding negative energy of search candidates. 
# ( It is of course unknown in practice. )
X, t = load_data()
 
# Normalize the mean and standard deviation along the each column of X to 0 and 1, respectively
#X = combo.misc.centering( X )

# Declare the class for calling the simulator. 
# In this tutorial, we simply refer to the value of t. 
# If you want to apply combo to other problems, you have to customize this class. 
class simulator:
    def __init__( self ):
        _, self.t = load_data()
    
    def __call__( self, action ):
        structure_current = X[action,:]
        log.write(str(self.t[action]))
        log.write('\t')
        log.write(str(structure_current))
        log.write('\n')
        log.flush()
        return self.t[action]

# Design of policy

# Declaring the policy by 
policy = combo.search.discrete.policy(test_X=X)
# test_X is the set of candidates which is represented by numpy.array.
# Each row vector represents the feature vector of the corresponding candidate

# set the seed parameter 
policy.set_seed( 5 )

# If you want to perform the initial random search before starting the Bayesian optimization, 
# the random sampling is performed by 

res = policy.random_search(max_num_probes=20, simulator=simulator())
# Input: 
# max_num_probes: number of random search 
# simulator = simulator
# output: combo.search.discreate.results (class)


# single query Bayesian search
# The single query version of COMBO is performed by 

res = policy.bayes_search(max_num_probes=1000, simulator=simulator(), score='EI', interval=20, num_rand_basis=5000)

# Input
# max_num_probes: number of searching by Bayesian optimization
# simulator: the class of simulator which is defined above
# score: the type of aquision funciton. TS, EI and PI are available
# interval: the timing for learning the hyper parameter. 
#               In this case, the hyper parameter is learned at each 20 steps
#               If you set the negative value to interval, the hyper parameter learning is not performed 
#               If you set zero to interval, the hyper parameter learning is performed only at the first step
# num_rand_basis: the number of basis function. If you choose 0,  ordinary Gaussian process runs


# The result of searching is summarized in the class combo.search.discrete.results.history()
# res.fx: observed negative energy at each step
# res.chosed_actions: history of choosed actions
# fbest, best_action= res.export_all_sequence_best_fx(): current best fx and current best action 
#                                                                                                   that has been observed until each step
# res.total_num_search: total number of search
#print 'f(x)='
#print res.fx[0:res.total_num_search]
#best_fx, best_action = res.export_all_sequence_best_fx()
#print 'current best'
#print best_fx
#print 'current best action='
#print best_action
#print 'history of chosed actions='
#print res.chosed_actions[0:res.total_num_search]


# save the results
res.save('test.npz')
np.savetxt('00_X.out', X, fmt='%d')

del res

end = time.time()
print("Execution Time: ", end - start)

