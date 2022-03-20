import mdts
import numpy as np
import itertools



### Thermal conductivity for (Si:Ge=1:1) alloy with 16 atoms in the interfacial structure
# DB_16=np.load('./DB_16.npy', allow_pickle=True).item()
# print(max(DB_16.values()))
# print(len(DB_16))

### Thermal conductivity for (Si:Si=1:1) alloy with 14 atoms in the interfacial structure
DB_14=np.load('./DB_14.npy', allow_pickle=True).item()
print(max(DB_14.values()))
print(len(DB_14))


### The get_reward simulates the structure: takes a python list represents the structure and return the reward
def get_reward(struct):
    s = ''.join(str(x) for x in struct)
    if s in DB_14.keys():
        cond = DB_14[s]
    else:
        print ("error")
    return cond

 
myTree=mdts.Tree(no_positions=14,   # number of positions in each structure. For example, 16 atoms.
                 atom_types=[0,1],  # atom types. For example, types of atoms: 0 for Si and 1 for Ge.
                 atom_const=[7,7],  # number of each atom type in the structure. For example, 8 atoms Si and 8 atoms Ge.
                 get_reward=get_reward,           # the experiment simulation function.
                 positions_order=list(range(14)), # define the order to assign atoms to the positions in the structure.
                 max_flag=False,     # if True the algorithm searches for maximum reward, else for minimum.
                 expand_children=2, # number of children to expand at each node, 1 means expand one child at a time.
                 play_out=5,        # number of play outs et each node.
                 play_out_selection="best", # when performing multiple playouts, best or mean is returned.
                 space=None,                # numpy ndarray representing the candidates space. 
                 candidate_pool_size=100,
                 ucb="mean",           # taking either average or best ucb score for MC tree search.
                 use_combo=False,      # weather to use Bayesian optimisation or not in combination with MC tree search.
                 combo_play_out=20,    # total number of candidates to be examind by COMBO.
                 combo_init_random=5,  # the initial random selection for Bayesian optimisation.
                 combo_step=5,         # the interval for Bayesian optimisation to perfrom hyperparameter optimization.
                 combo_lvl=5)          # the level of the tree at which start to apply Bayesian optimisation.
 
 
def output_func():
    # ### Optimal reward
    print (res.optimal_fx)

    # ### List of optimal candidates
    print (res.optimal_candidate)

    # ### List of tuples with the candidates examined and their rewards
    print (res.checked_candidates_DS)

    # ### Number of examined candidates
    print (res.checked_candidates_size)

    ### List of chosen candidates in order
    print (res.checked_candidates)

    ### List of simulated candidates rewards in order
    print (res.fx)

    ### List of current best reward
    print (res.best_fx)

    ### Maximum depth reached
    print (res.max_depth_reached)

    ### Number of nodes constructed
    print (res.no_nodes)

    ### Average visit per node
    print (res.avg_node_visit)


# def gen_res(itrs):
#     res_list = []
#     for i in range(1, itrs+1):
#         res_str = 'res'+str(i)
#         res_list.append(res_str)

#     return res_list


### Start the search for certain number of candidates and returns an object of type Result contains the result of the search
res=myTree.search(display=True,no_candidates=1010)
output_func()
res.save('results10.npz')


# for i in range(10):
    # res=myTree.search(display=True,no_candidates=20)

    ### Save results
    # res.save('results'+ str(i) +'.npz')
    # del(res)






# ### Load results
# new_res=mdts.Result()
# new_res.load('results.npz')
# print (new_res.optimal_fx)


# ### The tree can be saved
# mdts.save_tree(myTree,'Tree_file')
# del myTree

# ### Load the tree and continue the search
# myNewTree = mdts.load_tree('Tree_file')
# newTree_res=myNewTree.search(display=True, no_candidates=600)
# print (newTree_res.checked_candidates_size)
