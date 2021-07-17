import numpy as np
from BipartiteMatching import matching
from sampler import gibbs_sampler

def create_bipartite_graph(partial_order):
    graph = {}    
    for i in range(partial_order.shape[0]):        
        graph[i] = np.flatnonzero(partial_order[i, :] == 1)
        graph[i] = np.setdiff1d(graph[i], np.array([i]))
    return graph

def compute_minimal_chain_cover(partial_order):
    size = partial_order.shape[0]
    graph = create_bipartite_graph(partial_order)
    bipartite_matching = matching(graph)[0]
    chain_cover = []
    match_items = [it for it in bipartite_matching.values()]
    top_items = np.setdiff1d(np.arange(size), match_items)
    for it in top_items:
        curr_chain = [it]
        
        curr_it = it
        while True:
            if curr_it in bipartite_matching:
                curr_it = bipartite_matching[curr_it]
                curr_chain.append(curr_it)
            else:
                break
        chain_cover.append(np.array(curr_chain))
    return chain_cover

def chain_to_states(partial_order, chain):
    # Assumes that "chain" is actually a chain
    state_list = [[] for i in range(len(chain) + 1)]
    for it in chain:
        st = chain[partial_order[chain, it] == 1]
        state_list[len(st)] = st
    return state_list

def run_estimator(partial_order,
                  num_iter=1000,
                  p0=0.1):
    """
    Estimates the size of an ordinal knowledge space using a variation of the 
    subset simulation algorithm introduced by Au and Beck (2001).

    Parameters
    ----------
    partial_order : array-like, shape = [n_items, n_items]
        partial_order[i, j] = 1 if i <= j; 0 otherwise
    
    num_iter : int, optional
        Number of samples sets to generate at each level

    p0: float in (0, 1], optional 
        Target conditional probability estimate

    Returns
    -------
    num_states : float
        Estimated number of states in the knowledge space

    Ei_dict : dictionary
        Conditional probability estimate for each level

    References
    ----------
    Au, S.-K. and Beck, J. L. (2001). Estimation of small failure probabilities 
        in high dimensions by subset simulation. Probabilistic Engineering 
        Mechanics, 16(4):263 277.
    """    
    cover = compute_minimal_chain_cover(partial_order)
    size = partial_order.shape[0]
    p0_thresh = int(num_iter*p0)
    chain_state_dict = {}
    chain_size = {}
    chain_ordered_list = {}
    num_chains = len(cover)
    chain = []
    chain_index = [1]

    for i in range(num_chains):
        chain_size[i] = len(cover[i])
        chain_state_dict[i] = chain_to_states(partial_order, cover[i])
        chain_ordered_list[i] = []
        for j in range(chain_size[i]):
            # Add 1 since indices will be used in Fortran arrays
            chain_ordered_list[i].append(np.setdiff1d(
                chain_state_dict[i][j + 1], chain_state_dict[i][j])[0] + 1)
        if i > 0:
            chain_index.append(len(chain) + 1)
        chain.extend(chain_ordered_list[i])
    chain = np.array(chain, dtype='int')
    chain_index = np.array(chain_index, dtype='int')

    Ei = size - 1
    Eii = Ei    
    curr_set = np.random.randint(0, 2, size) # Random starting set
    Ei_dict = {}
    card_list = []
    card_diff = np.zeros(1, dtype='int')
    last_iter = False
    while True:
        print("Sampling with Ei = " + str(Ei) + "...")
        Ei_list = [] 
        diff_list = []

        for i in range(num_iter):
            card_diff[0] = int(0)
            # After running gibbs_sampler, curr_set is updated with the newly
            # sampled set from E_i, while card_diff contains the distance from
            # curr_set to the smallest knowledge state containing it
            gibbs_sampler(curr_set, partial_order, chain, chain_index, Ei,
                          card_diff, num_chains, size)           
            Ei_list.append(curr_set.copy())
            diff_list.append(card_diff[0])

        # Update E_ii value using p0 threshold
        Eii = np.sort(diff_list)[p0_thresh]
        if Eii == Ei:
            Eii = Ei - 1
        # Find all sets in E_ii (i.e., the sets with a distance at most Eii from
        # the nearest knowledge state)
        seed_list = np.array(Ei_list)[np.flatnonzero(diff_list <= Eii), :]
        # Compute the proportion of sets in E_ii
        Ei_dict[Eii] = len(seed_list)/num_iter

        if Ei == 1:
            break
        # Update family of sets and randomly choose one of the sets in E_ii
        Ei = Eii                   
        curr_set = seed_list[np.random.randint(len(seed_list))].copy()

    # Compute estimated probability of sampling a knowledge state 
    prob_prod = 1.0
    for prob in list(Ei_dict.values()):
        prob_prod *= prob
    # Compute size of the Cartesian product of the chains
    card_prod = 1.0
    for i in range(len(cover)):
        card_prod *= (len(cover[i]) + 1)
    # Estimate of the knowledge state size
    num_states = prob_prod * card_prod
    
    return num_states, Ei_dict
