#!/usr/bin/env python
"""Fast implementation of UniFrac for use with very large datasets"""

from random import shuffle
from numpy import ones, ma, where
from numpy.random import permutation
from cogent.maths.unifrac.fast_tree import *
# not imported by import *
from cogent.maths.unifrac.fast_tree import _weighted_unifrac, _branch_correct 
from cogent.parse.tree import DndParser
from cogent.cluster.metric_scaling import *
from cogent.core.tree import PhyloNode, TreeError
from cogent.cluster.UPGMA import UPGMA_cluster
from cogent.phylo.nj import nj
from StringIO import StringIO

__author__ = "Rob Knight and Micah Hamady"
__copyright__ = "Copyright 2007-2011, The Cogent Project"
__credits__ = ["Rob Knight", "Micah Hamady", "Daniel McDonald", 
    "Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.6.0.dev"
__maintainer__ = "Rob Knight, Micah Hamady"
__email__ = "rob@spot.colorado.edu, hamady@colorado.edu"
__status__ = "Prototype"

UNIFRAC_DIST_MATRIX = "distance_matrix"
UNIFRAC_PCOA = "pcoa"
UNIFRAC_CLUST_ENVS = "cluster_envs"
UNIFRAC_NJ_ENVS = "nj_envs"
UNIFRAC_DIST_VECTOR = "distance_vector"
UNIFRAC_VALID_MODES = set([UNIFRAC_DIST_MATRIX, UNIFRAC_PCOA, 
    UNIFRAC_CLUST_ENVS,UNIFRAC_DIST_VECTOR, UNIFRAC_NJ_ENVS])
UNIFRAC_DEFAULT_MODES = set([UNIFRAC_DIST_MATRIX, UNIFRAC_PCOA,
    UNIFRAC_CLUST_ENVS])

TEST_ON_PAIRWISE = "Pairwise"
TEST_ON_TREE = "Tree"
TEST_ON_ENVS = "Envs"

def identity(x): return range(x)

def num_comps(num_envs):
    """ Calc number of comparisons for Bonferroni correction """
    return (num_envs * (num_envs-1.0)) / 2.0

def mcarlo_sig(real_val, sim_vals, num_comps, tail='low'):
    """ Calc significance value for Monte Carlo simulation
   
    real_val: real value of test
    sim_vals: numpy array of values from simulation runs
    num_envs: number of environments (for correction)
    tail: which side of distro to check, high=larger, low=smaller

    returns (raw pval, corrected pval)
    """
    if num_comps < 1:
        raise ValueError, "num_comps must be > 0"
    sim_vals = array(sim_vals) 
    out_count = 0
    pop_size = float(len(sim_vals))
    if tail == 'low':
        out_count = sum(sim_vals < real_val) 
    elif tail == 'high':
        out_count = sum(sim_vals > real_val) 
    else:
        raise ValueError, "tail must be 'low' or 'high'"
    raw_pval = out_count/pop_size
    cor_pval = raw_pval * num_comps
    # reset to 1.0 if corrected if > 1.0 
    if cor_pval > 1.0:
        cor_pval = 1.0
    elif raw_pval == 0.0:
        cor_pval = "<=%.1e" % (1.0/pop_size)
    return (raw_pval, cor_pval)

def fast_unifrac_file(tree_in, envs_in, weighted=False, metric=unifrac, is_symmetric=True, modes=UNIFRAC_DEFAULT_MODES):
    """Takes tree and envs file and returns fast_unifrac() results

    typical results: distance matrix, UPGMA cluster and PCoA

    tree_in: open file object or list of lines with tree in Newick format
    envs_in: open file object of list of lines with tab delimited sample mapping file
    see fast_unifrac() for further description
    """
    tree = DndParser(tree_in, UniFracTreeNode)
    envs = count_envs(envs_in)
    return fast_unifrac(tree, envs, weighted, metric, is_symmetric=is_symmetric, modes=modes)

def fast_unifrac_permutations_file(tree_in, envs_in, weighted=False, 
    num_iters=1000, verbose=False, test_on=TEST_ON_PAIRWISE):
    """ Wrapper to read tree and envs from files. """
    result = []
    t = DndParser(tree_in, UniFracTreeNode)
    envs = count_envs(envs_in)
    if test_on == TEST_ON_PAIRWISE:
        # calculate real values
        results = fast_unifrac(t, envs,weighted=weighted, metric=unifrac,
            is_symmetric=True, modes=[UNIFRAC_DIST_MATRIX])
        real_env_mat, unique_envs = results[UNIFRAC_DIST_MATRIX]
        num_uenvs = real_env_mat.shape[0] 
        cur_num_comps = num_comps(num_uenvs)
        for i in range(num_uenvs): 
            first_env = unique_envs[i]
            for j in range(i+1, num_uenvs): 
                second_env = unique_envs[j]
                real = real_env_mat[i][j]
                sim = fast_unifrac_permutations(t, envs, weighted, num_iters, 
                    first_env=first_env, second_env=second_env)
                raw_pval, cor_pval = mcarlo_sig(real, sim, cur_num_comps, 
                    tail='high')
                result.append((first_env, second_env, raw_pval, cor_pval))
                if verbose:
                    print "env %s vs %s" % (first_env, second_env)
                    print raw_pval, cor_pval, num_uenvs, cur_num_comps, 'high'
    # calculate single p-value for whole tree
    elif test_on == TEST_ON_TREE:
        # will be using env_unique_fraction
        real_ufracs, sim_ufracs = fast_unifrac_whole_tree(t, envs, num_iters)
        raw_pval, cor_pval = mcarlo_sig(sum(real_ufracs), [sum(x) for x in sim_ufracs], 1, tail='high')
        result.append(('whole tree', raw_pval, cor_pval))
    # calculate one p-value per env 
    elif test_on == TEST_ON_ENVS:
        unique_envs, num_uenvs = get_unique_envs(envs)
        real_ufracs, sim_ufracs = fast_unifrac_whole_tree(t, envs, num_iters)
        sim_m = array(sim_ufracs) 
        # for each env, cal paval
        for i in range(len(real_ufracs)):
            raw_pval, cor_pval = mcarlo_sig(real_ufracs[i], sim_m[:,i], 1, 
                tail='high')
            result.append((unique_envs[i], raw_pval, cor_pval))
    else:
        raise ValueError, "Invalid test_on value: %s" % str(test_on)
    return result

def fast_p_test_file(tree_in, envs_in, num_iters=1000, verbose=False, 
    test_on=TEST_ON_PAIRWISE):
    """ Wrapper to read tree and envs from files. """
    result = []
    t = DndParser(tree_in, UniFracTreeNode)
    envs = count_envs(envs_in)
    unique_envs, num_uenvs = get_unique_envs(envs)
    # calculate real, sim vals and p-vals for each pair of envs in tree 
    if test_on == TEST_ON_PAIRWISE:
        cur_num_comps = num_comps(num_uenvs) 
        for i in range(num_uenvs): 
            first_env = unique_envs[i]
            for j in range(i+1, num_uenvs): 
                second_env = unique_envs[j]
                real = fast_p_test(t, envs, num_iters=1, first_env=first_env, 
                    second_env=second_env, permutation_f=identity)[0]
                sim = fast_p_test(t, envs, num_iters, first_env=first_env, 
                    second_env=second_env)
                raw_pval, cor_pval = mcarlo_sig(real, sim, cur_num_comps, 
                    tail='low')
                result.append((first_env, second_env, raw_pval, cor_pval))
                if verbose:
                    print "P Test: env %s vs %s" % (first_env, second_env)
                    print raw_pval, cor_pval, num_uenvs, cur_num_comps, 'low'
    # calculate real, sim vals and p-vals for whole tree
    elif test_on == TEST_ON_TREE:
        real = fast_p_test(t, envs, num_iters=1, permutation_f=identity)[0]
        sim = fast_p_test(t, envs, num_iters)
        raw_pval, cor_pval = mcarlo_sig(real, sim, 1, tail='low')
        result.append(('Whole Tree', raw_pval))
    else:
        raise ValueError, "Invalid test_on value: %s" % str(test_on)

    return result

def _fast_unifrac_setup(t, envs, make_subtree=True):
    """Setup shared by fast_unifrac and by significance tests."""
    if make_subtree:
        t2 = t.copy()
        wanted = set(envs.keys())
        all_tips = list(t2.tips())
        # can't delete nodes while itering, thus list()
        for tip in all_tips: 
            if tip.Name not in wanted:
                curr_node = tip.Parent
                did_remove = curr_node.removeNode(tip)
                if not did_remove:
                    raise RuntimeError('failed to remove tip in tree')
                # and travel up tree 
                while curr_node.istip():
                    new_node = curr_node.Parent
                    # better to ask forgiveness that permission
                    try:
                        new_node.removeNode(curr_node)
                    except AttributeError:
                        if curr_node.isroot():
                            break
                        else:
                            raise
                    curr_node = new_node
        t2.prune()
        t = t2
    #index tree
    node_index, nodes = index_tree(t)
    #get good nodes, defined as those that are in the env file.
    good_nodes=dict([(i.Name,envs[i.Name]) for i in t.tips() if i.Name in envs])
    envs = good_nodes
    count_array, unique_envs, env_to_index, node_to_index = index_envs(envs, node_index)
    env_names = sorted(unique_envs)
    #Note: envs get sorted at the step above
    branch_lengths = get_branch_lengths(node_index)
    if not envs:
        raise ValueError, "No valid samples/environments found. Check whether tree tips match otus/taxa present in samples/environments"
    return envs, count_array, unique_envs, env_to_index, node_to_index, env_names, branch_lengths, nodes, t

def fast_unifrac_whole_tree(t, envs, num_iters, permutation_f=permutation):
    """Performs UniFrac permutations on whole tree """
    sim_ufracs = []
    envs, count_array, unique_envs, env_to_index, node_to_index, env_names, \
        branch_lengths, nodes, t = _fast_unifrac_setup(t, envs)
    
    bound_indices = bind_to_array(nodes, count_array)
    orig_count_array = count_array.copy()

    # calculate real values 
    bool_descendants(bound_indices)
    real_bl_sums, real_bl_ufracs = env_unique_fraction(branch_lengths, 
        count_array)
    tip_indices = [n._leaf_index for n in t.tips()]
    for i in range(num_iters):
        permute_selected_rows(tip_indices, orig_count_array, count_array, 
            permutation_f)
        bool_descendants(bound_indices)
        cur_bl_sums, cur_bl_ufracs = env_unique_fraction(branch_lengths, 
            count_array)
        sim_ufracs.append(cur_bl_ufracs)
    return real_bl_ufracs, sim_ufracs 

def PD_whole_tree(t, envs):
    """Run PD on t and envs for each env.

    Note: this is specific for PD per se, use PD_generic_whole_tree if you
    want to calculate a related metric.
    """
    envs, count_array, unique_envs, env_to_index, node_to_index, env_names, \
        branch_lengths, nodes, t = _fast_unifrac_setup(t, envs)
    count_array = count_array.astype(bool)
    bound_indices = bind_to_array(nodes, count_array)
    #initialize result
    bool_descendants(bound_indices)
    result = (branch_lengths * count_array.T).sum(1)
    return unique_envs, result

def PD_generic_whole_tree(t, envs, metric=PD):
    """Run PD on t and envs for each env.

    Note: this is specific for PD per se, use PD_generic_whole_tree if you
    want to calculate a related metric.
    """
    envs, count_array, unique_envs, env_to_index, node_to_index, env_names, \
        branch_lengths, nodes, t = _fast_unifrac_setup(t, envs)
    count_array = count_array.astype(bool)
    bound_indices = bind_to_array(nodes, count_array)
    #initialize result
    bool_descendants(bound_indices)
    result = PD_vector(branch_lengths, count_array,metric)
    return unique_envs, result

def fast_unifrac_permutations(t, envs, weighted, num_iters, first_env, 
    second_env, permutation_f=permutation, unifrac_f=_weighted_unifrac):
    """Performs UniFrac permutations between specified pair of environments.
    
    NOTE: this function just gives you the result of the permutations, need to 
    compare to real values from doing a single unifrac.
    """
    result = []
    envs, count_array, unique_envs, env_to_index, node_to_index, env_names, branch_lengths, nodes, t = _fast_unifrac_setup(t, envs)

    first_index,second_index = env_to_index[first_env], env_to_index[second_env]
    count_array = count_array[:,[first_index,second_index]] #ditch rest of array
    bound_indices = bind_to_array(nodes, count_array)
    orig_count_array = count_array.copy()
    first_col, second_col = count_array[:,0], count_array[:,1]
    tip_indices = [n._leaf_index for n in t.tips()]

    #figure out whether doing weighted or unweighted analysis: for weighted,
    #need to figure out root-to-tip distances, but can skip this step if
    #doing unweighted analysis.
    if weighted:
        tip_ds = branch_lengths.copy()[:,newaxis]
        bindings = bind_to_parent_array(t, tip_ds)
        tip_distances(tip_ds, bindings, tip_indices)
        if weighted == 'correct':
            bl_correct = True
        else:
            bl_correct = False
        first_sum, second_sum = [sum(take(count_array[:,i], tip_indices)) for i in range(2)]
        for i in range(num_iters):
            permute_selected_rows(tip_indices, orig_count_array, count_array, permutation_f)
            sum_descendants(bound_indices)
            curr = unifrac_f(branch_lengths, first_col, second_col, first_sum, second_sum)
            if bl_correct:
                curr /= _branch_correct(tip_ds, first_col, second_col, first_sum, second_sum)
            result.append(curr)
    else:
        for i in range(num_iters):
            permute_selected_rows(tip_indices, orig_count_array, count_array, permutation_f)

            # commenting out these (??)
            #first = count_array.copy()
            #permute_selected_rows(tip_indices, orig_count_array, count_array, permutation_f)
            #second = count_array.copy()

            bool_descendants(bound_indices)
            curr = unifrac(branch_lengths, first_col, second_col)
            result.append(curr)
    return result

def fast_p_test(t, envs, num_iters, first_env=None, second_env=None, 
    permutation_f=permutation):
    """Performs Andy Martin's p test between specified pair of environments.

    t: tree 
    envs: envs 
    first_env: name of first env, or None if doing whole tree
    second_env: name of second env, or None if doing whole tree

    NOTE: this function just gives you the result of the permutations, need to 
    compare to real Fitch parsimony values. Sleazy way to get the real values 
    is to set num_iters to 1, permutation_f to identity."""
    result = []
    envs, count_array, unique_envs, env_to_index, node_to_index, env_names, branch_lengths, nodes, t = _fast_unifrac_setup(t, envs)

    # check if doing pairwise
    if not (first_env is None or second_env is None):
        first_ix, second_ix = env_to_index[first_env], env_to_index[second_env]
        count_array = count_array[:,[first_ix,second_ix]] #ditch rest of array
    # else if both envs are none, we're doing whole tree 
    elif not (first_env is None and second_env is None):
        raise ValueError, "Both envs must either have a value or be None."

    bound_indices = bind_to_array(nodes, count_array)
    orig_count_array = count_array.copy()
    tip_indices = [n._leaf_index for n in t.tips()]
    for i in range(num_iters):
        count_array *= 0
        permute_selected_rows(tip_indices, orig_count_array, count_array, 
            permutation_f=permutation_f)
        curr = fitch_descendants(bound_indices)
        result.append(curr)
    return result

def shared_branch_length(t, envs, env_count=1):
    """Returns the shared branch length env_count combinations of envs
    
    t: phylogenetic tree relating the sequences.
    envs: dict of {sequence:{env:count}} showing environmental abundance.
    env_count: number of envs that must be within the subtree

    Returns {(env1,env2,...env_count):shared_branch_length}
    """

    envs, count_array, unique_envs, env_to_index, node_to_index, env_names, branch_lengths, nodes, t = _fast_unifrac_setup(t, envs)

    

    if len(unique_envs) < env_count:
        raise ValueError, "Not enough environments for env_count"

    index_to_env = dict([(i,e) for i,e in enumerate(unique_envs)])

    bound_indices = bind_to_array(nodes, count_array)
    bool_descendants(bound_indices)

    # determine what taxa meet the required number of environments
    count_array = where(count_array > 0, 1, 0)
    counts = count_array.sum(axis=1)
    taxa_to_investigate = (counts == env_count).nonzero()[0]

    # determine what environments corrispond to what taxa
    envs_to_investigate = {} 
    for row_index in taxa_to_investigate:
        taxa_envs = count_array[row_index]
        row_envs = tuple([index_to_env[i] for i,v in enumerate(taxa_envs) if v])
        try: 
            envs_to_investigate[row_envs].append(row_index)
        except KeyError:
            envs_to_investigate[row_envs] = [row_index]
            

    # compute shared branch length for each environments
    result = {}
    for envs_tuple, taxa_indices in envs_to_investigate.items():
        valid_rows = zeros(len(count_array))
        for i in taxa_indices:
            valid_rows[i] = 1.0
        result[envs_tuple] = sum(branch_lengths * valid_rows)

    return result

def shared_branch_length_to_root(t, envs):
    """Returns the shared branch length for a single env from tips to root
    
    t: phylogenetic tree relating sequences
    envs: dict of {sequence:{env:count}} showing environmental abundance

    Returns {env:shared_branch_length}
    """
    working_t = t.copy()
    result = {}

    # decorate nodes with environment information
    for n in working_t.postorder():
        # for tip, grab and set env information
        if n.isTip():  
            curr_envs = envs.get(n.Name, None)
            if curr_envs is None:
                n.Envs = set([])
            else:
                n.Envs = set(curr_envs.keys())
        # for internal node, collect descending env information
        else:  
            n.Envs = set([])  # should only visit each internal node once
            for c in n.Children:
                n.Envs.update(c.Envs)

    # collect branch length for each environment
    for n in working_t.preorder(include_self=False):
        if not hasattr(n, 'Length') or n.Length is None:
            continue

        for e in n.Envs:
            if e not in result:
                result[e] = 0.0
            result[e] += n.Length

    return result

def fast_unifrac(t, envs, weighted=False, metric=unifrac, is_symmetric=True, 
    modes=UNIFRAC_DEFAULT_MODES, weighted_unifrac_f=_weighted_unifrac,make_subtree=True):
    """ Run fast unifrac.
    
    t: phylogenetic tree relating the sequences.  pycogent phylonode object
    envs: dict of {sequence:{env:count}} showing environmental abundance.
    weighted: if True, performs the weighted UniFrac procedure.
    metric: distance metric to use.  currently you must use unifrac only
        if weighted=True.
        see fast_tree.py for metrics (e.g.: G, unnormalized_G, unifrac, etc.)
    modes: tasks to perform on running unifrac.  see fast_unifrac.py
        default is to get a unifrac distance matrix, pcoa on that matrix, and a
        cluster of the environments
    is_symmetric: if the desired distance matrix is symmetric 
        (dist(sampleA, sampleB) == dist(sampleB, sampleA)), then set this True
        to prevent calculating the same number twice

    using default modes, returns a dictionary with the following (key:value) pairs:

    'distance_matrix': a tuple with a numpy array of pairwise distances between 
    samples and a list of names describing the order of samples in the array 
    
    'cluster_envs': cogent.core.PhyloNode object containing results of running 
    UPGMA on the distance matrix.
    
    'pcoa': a cogent.util.Table object with the results of running Principal 
    Coordines Analysis on the distance matrix.

    Usage examples: (these assume the example files exist)
    from cogent.parse.tree import DndParser
    from cogent.maths.unifrac.fast_unifrac import count_envs, fast_unifrac
    from cogent.maths.unifrac.fast_tree import UniFracTreeNode, count_envs, G
    tree_in = open('Crump_et_al_example.tree')
    envs_in = open('Crump_et_al_example_env_file.txt')
    tree = DndParser(tree_in, UniFracTreeNode)
    envs = count_envs(envs_in)
    unifrac_result = fast_unifrac(tree, envs)
    G_result = fast_unifrac(tree, envs, metric=G, is_symmetric=False)

    WARNING: PCoA on asymmetric matrices (e.g.: G metric) is meaningless
    because these are not distance matrices.
    """
    modes = set(modes)  #allow list, etc. of modes to be passed in.
    if not modes or modes - UNIFRAC_VALID_MODES:
        raise ValueError, "Invalid run modes: %s, valid: %s" % (str(modes),str(UNIFRAC_VALID_MODES))

    envs, count_array, unique_envs, env_to_index, node_to_index, env_names, branch_lengths, nodes, t = _fast_unifrac_setup(t, envs, make_subtree)
    bound_indices = bind_to_array(nodes, count_array)
    #initialize result
    result = {}
    
    #figure out whether doing weighted or unweighted analysis: for weighted,
    #need to figure out root-to-tip distances, but can skip this step if
    #doing unweighted analysis.
    if weighted:
        tip_indices = [n._leaf_index for n in t.tips()]
        sum_descendants(bound_indices)
        tip_ds = branch_lengths.copy()[:,newaxis]
        bindings = bind_to_parent_array(t, tip_ds)
        tip_distances(tip_ds, bindings, tip_indices)
        if weighted == 'correct':
            bl_correct = True
        else:
            bl_correct = False
        u = weighted_unifrac_matrix(branch_lengths, count_array, tip_indices, \
            bl_correct=bl_correct, tip_distances=tip_ds, \
            unifrac_f=weighted_unifrac_f)
        #figure out if we need the vector
        if UNIFRAC_DIST_VECTOR in modes:
            result[UNIFRAC_DIST_VECTOR] = (weighted_unifrac_vector(
                branch_lengths, count_array, tip_indices, 
                bl_correct=bl_correct, tip_distances=tip_ds, 
                unifrac_f=weighted_unifrac_f), env_names)
    else:
        bool_descendants(bound_indices)
        u = unifrac_matrix(branch_lengths, count_array, metric=metric, is_symmetric=is_symmetric)
        if UNIFRAC_DIST_VECTOR in modes:
            result[UNIFRAC_DIST_VECTOR] = (unifrac_vector(branch_lengths, 
                count_array), env_names)
    
    #check if we have to do the matrix calculations, which are expensive
    if modes - set([UNIFRAC_DIST_VECTOR]):
        result.update(unifrac_tasks_from_matrix(u, env_names, modes=modes))
    return result

def fast_unifrac_one_sample(one_sample_name, t, envs, weighted=False, 
    metric=unifrac, weighted_unifrac_f=_weighted_unifrac, make_subtree=False):
    """ performs fast unifrac between a specified sample and all other samples
    
    one_sample_name: unifrac will be calculated bet this and each other sample
    t: phylogenetic tree relating the sequences.  pycogent phylonode object
    envs: dict of {sequence:{env:count}} showing environmental abundance.
    weighted: if True, performs the weighted UniFrac procedure.  If "correct",
        performs weighted unifrac with branch length correction.  If false
        (default), performs supplied metric (default: unweighted unifrac)
    metric: distance metric to use, unless weighted=True
        see fast_tree.py for metrics (e.g.: G, unnormalized_G, unifrac, etc.)

    returns distances, environments. e.g. when one_sample_name = 'B': 
    array([ 0.623, 0., 0.4705]), ['A', 'B', 'C']
    returns a ValueError if we have no data on one_sample_name
    """
    envs, count_array, unique_envs, env_to_index, node_to_index, env_names, \
        branch_lengths, nodes, t = _fast_unifrac_setup(t, envs, make_subtree)
    bound_indices = bind_to_array(nodes, count_array)
    result = {}
    try:
        one_sample_idx = env_names.index(one_sample_name)
    except ValueError:
        raise ValueError('one_sample_name not found, ensure that there are'+\
            ' tree tips and corresponding envs counts for sample ' +\
            one_sample_name )

    if weighted:
        tip_indices = [n._leaf_index for n in t.tips()]
        sum_descendants(bound_indices)
        tip_ds = branch_lengths.copy()[:,newaxis]
        bindings = bind_to_parent_array(t, tip_ds)
        tip_distances(tip_ds, bindings, tip_indices)
        if weighted == 'correct':
            bl_correct = True
        else:
            bl_correct = False
        u = weighted_one_sample(one_sample_idx, branch_lengths, count_array,
            tip_indices, bl_correct=bl_correct, tip_distances=tip_ds, \
            unifrac_f=weighted_unifrac_f)
    else: # unweighted
        bool_descendants(bound_indices)
        u = unifrac_one_sample(one_sample_idx, branch_lengths, count_array, 
            metric=metric)

    return (u, env_names)
        
def unifrac_tasks_from_matrix(u, env_names, modes=UNIFRAC_DEFAULT_MODES):
    """Returns the UniFrac matrix, PCoA, and/or cluster from the matrix."""
    result = {}

    if UNIFRAC_DIST_MATRIX in modes:
        result[UNIFRAC_DIST_MATRIX] = (u, env_names)

    if UNIFRAC_PCOA in modes:
        point_matrix, eigvals = principal_coordinates_analysis(u)
        result[UNIFRAC_PCOA] =  output_pca(point_matrix, eigvals, env_names)

    if UNIFRAC_CLUST_ENVS in modes:
        nodes = map(PhyloNode, env_names)
        BIG = 1e305
        U = u.copy()
        for i in range(len(U)):
            U[i,i] = BIG
        c = UPGMA_cluster(U, nodes, BIG)
        result[UNIFRAC_CLUST_ENVS] = c

    if UNIFRAC_NJ_ENVS in modes:
        c = nj(dists_to_nj(u, env_names))
        result[UNIFRAC_NJ_ENVS] = c

    return result

def unifrac_recursive(tree, envs, ref_tree):#, metric=weighted):
    """Performs UniFrac recursively over a tree.

    Specifically, for each node in the tree, performs UniFrac clustering.
    Then compares the UniFrac tree to a reference tree of the same taxa using
    the tip-to-tip distances and the subset distances. Assumption is that if
    the two trees match, the node represents a group in which evolution has
    mirrored the evolution of the reference tree.

    tree: contains the tree on which UniFrac will be performed recursively.
    envs: environments for UniFrac clustering (these envs should match the
          taxon labels in the ref_tree)
    ref_tree: reference tree that the clustering is supposed to match.
    metric: metric for UniFrac clustering.

    Typically, will want to estimate significance by comparing the actual
    values from ref_tree to values obtained with one or more shuffled versions
    of ref_tree (can make these with permute_tip_labels).
    """
    lengths, dists, sets = [], [], []
    for node in tree.traverse(self_before=True, self_after=False):
        try:
            result = fast_unifrac(node, envs, weighted=False, modes=set([UNIFRAC_CLUST_ENVS]))
            curr_tree = result[UNIFRAC_CLUST_ENVS]
        except AttributeError:
            #hit a zero branch length
            continue
        if curr_tree is None:
            #hit single node?
            continue
        try:
            l = len(curr_tree.tips())
            d = curr_tree.compareByTipDistances(ref_tree)
            s = curr_tree.compareBySubsets(ref_tree, True)
            #want to calculate all values before appending so we can bail out
            #if any of the calculations fails: this ensures that the lists
            #remain synchronized.
            lengths.append(l)
            dists.append(d)
            sets.append(s)
        except ValueError:
            #no common taxa
            continue
    return lengths, dists, sets

def dists_to_nj(matrix, labels):
    """Wraps matrix and labels together for format NJ requires."""
    result = {}
    for outer, row in zip(labels, matrix):
        for inner, i in zip(labels, row):
            result[(outer, inner)] = i
    return result

def shuffle_tipnames(t):
    """Returns copy of tree t with tip names shuffled."""
    result = t.copy()
    names = [i.Name for i in t.tips()]
    shuffle(names)
    for name, tip in zip(names, result.tips()):
        tip.Name = name
    return result

def weight_equally(tree_list, envs_list):
    """Returns equal weights for all trees."""
    num_trees = len(tree_list)
    return ones(num_trees)

def weight_by_num_tips(tree_list, envs_list):
    """Weights each tree by the number of tips it contains."""
    return array([len(list(t.tips())) for t in tree_list])

def weight_by_branch_length(tree_list, envs_list):
    """Weights each tree by the sum of its branch length."""
    return array([sum(filter(None, [i.Length for i in \
        t.traverse(self_before=True,self_after=False) if hasattr(i, 'Length')])) for t in tree_list])

def weight_by_num_seqs(tree_list, envs_list):
    """Weights each tree by the number of seqs it contains."""
    return array(map(sum_env_dict, envs_list))

def get_all_env_names(envs):
    """Returns set of all env names from envs."""
    result = set()
    for e in envs.values():
        result.update(e.keys())
    return result

def consolidate_skipping_missing_matrices(matrices, env_names, weights, 
    all_env_names):
    """Consolidates matrices, skipping any that are missing envs"""
    weight_sum = 0
    result = zeros((len(all_env_names),len(all_env_names)), float)
    for m, e, w in zip(matrices, env_names, weights):
        if e == all_env_names: #note -- assumes sorted
            result += m * w
            weight_sum += w
    #readjust weights for missing matrices
    result /= weight_sum
    return result

def consolidate_missing_zero(matrices, env_names, weights, all_env_names):
    """Consolidates matrices, setting missing values to 0 distance"""
    result = zeros((len(all_env_names),len(all_env_names)), float)
    for m, e, w in zip(matrices, env_names, weights):
        result += reshape_by_name(m, e, all_env_names, 0) * w
    return result

def consolidate_missing_one(matrices, env_names, weights, all_env_names):
    """Consolidates matrices, setting missing values to 1 distance"""
    result = zeros((len(all_env_names),len(all_env_names)), float)
    for m, e, w in zip(matrices, env_names, weights):
        result += reshape_by_name(m, e, all_env_names, 1) * w
    return result

def consolidate_skipping_missing_values(matrices, env_names, weights, 
    all_env_names):
    """Consolidates matrices, skipping only values from missing envs"""
    result = []
    for m, e, w in zip(matrices, env_names, weights):
        reshaped = reshape_by_name(m, e, all_env_names, masked=True)
        reshaped *= w
        result.append(reshaped)

    data = array([i.data for i in result], float)
    masks = array([i.mask for i in result], bool)
    masked_result = ma.array(data, mask=masks)
    #figure out mask of weights so we can figure out per-element weighting
    masked_weights = ma.array(zeros(data.shape), mask=masks) + \
        array(weights,float).reshape((len(weights),1,1))
    return masked_result.sum(0)/masked_weights.sum(0)

def reshape_by_name(m, old_names, new_names, default_off_diag=0,default_diag=0, 
    masked=False):
    """Reshape matrix m mapping slots from old names to new names.  """
    num_names = len(new_names)
    result = zeros((num_names,num_names), float) + default_off_diag
    for i in range(num_names):
        result[i,i] = default_diag
    pairs = {}
    for i, n in enumerate(old_names):
        if n in new_names:
            pairs[i] = new_names.index(n)
    for i, row in enumerate(m):
        new_i = pairs[i]
        for j, val in enumerate(row):
            new_j = pairs[j]
            result[new_i, new_j] = val
    if masked:
        mask = ones((num_names, num_names), float)
        for i in pairs.values():
            for j in pairs.values():
                mask[i,j] = 0
        result = ma.array(result, mask=mask)
    return result

def meta_unifrac(tree_list, envs_list, weighting_f, 
    consolidation_f=consolidate_skipping_missing_values, 
    modes=UNIFRAC_DEFAULT_MODES, **unifrac_params):
    """Perform metagenomic UniFrac on a list of trees and envs.

    tree_list: list of tree objects
    env_list: list of sample x env count arrays
    weighting_f: f(trees, envs) -> weights
    consolidation_f: f(matrix_list, name_list, weight_list, all_env_names) -> matrix
    unifrac_params: parameters that will be passed to unifrac to build the matrix

    Notes:
    - tree list and env list must be same length and consist of matched pairs,
      i.e. each tree must have a corresponding envs array.
    """
    all_weights = weighting_f(tree_list, envs_list)
    all_env_names = set()
    for e in envs_list:
        all_env_names.update(get_all_env_names(e))
    all_env_names = sorted(all_env_names)
    matrices = []
    env_names = []
    weights = []
    #need to be robust to failure to build UniFrac tree, and to combine the
    #matrices that survive in a reasonable way that handles missing data
    for t, e, w in zip(tree_list, envs_list, all_weights):
        #try:
            u, en = fast_unifrac(t, e, modes=[UNIFRAC_DIST_MATRIX], **unifrac_params)[UNIFRAC_DIST_MATRIX]
            matrices.append(u)
            env_names.append(en)
            weights.append(w)
        #except ValueError:
        #    pass
    #normalize weights so sum to 1
    weights = array(weights, float)/sum(weights)
    final_matrix = consolidation_f(matrices, env_names, weights, all_env_names)
    return unifrac_tasks_from_matrix(final_matrix, all_env_names, modes)

# Crump example tree 
CRUMP_TREE = """((((((((((((((((AF141409:0.25346,(((AF141563:0.00000,AF141568:0.00000):0.03773,AF141410:0.06022)Pseudomonadaceae:0.08853,((AF141521:0.12012,AF141439:0.03364):0.01884,(((((((((((AF141494:0.00000,AF141503:0.00000):0.00000,AF141553:0.00000):0.00000,AF141555:0.00000):0.00000,AF141527:0.00000):0.00000,AF141477:0.00513):0.00000,AF141554:0.00000):0.00000,AF141438:0.00000):0.00000,AF141493:0.00000):0.00000,AF141489:0.00000):0.00000,AF141437:0.00513):0.00000,AF141490:0.00771):0.07637):0.00573):0.00896):0.00327,((AF141546:0.03666,AF141510:0.04080):0.03273,AF141586:0.05358):0.00492):0.02041,((AF141577:0.04055,((AF141558:0.00000,AF141569:0.00771):0.01436,AF141430:0.01282):0.01144)Legionellaceae:0.06630,(AF141501:0.01323,AF141528:0.00000)SUP05:0.06854):0.02220):0.00651,AF141498:0.12862):0.00081,AF141428:0.08582):0.00570,((((((((((AF141464:0.00485,AF141422:0.00256):0.00000,AF141473:0.00770)Alcaligenaceae:0.06178,((((((AF141492:0.01026,AF141453:0.00513):0.00000,AF141443:0.01027):0.00323,(AF141441:0.00162,(AF141404:0.00000,AF141398:0.00000):0.00161):0.02267):0.01703,((AF141465:0.00513,AF141491:0.00256):0.00000,AF141405:0.00256):0.01214)Polynucleobacter:0.02121,AF141431:0.03263)Ralstoniaceae:0.01786,AF141456:0.01328):0.01630):0.00325,((((((((AF141462:0.01064,AF141414:0.00513):0.00000,AF141457:0.00000):0.00161,AF141413:0.00383):0.00646,AF141459:0.08285):0.02827,((((((AF141468:0.00256,AF141392:0.00256):0.00324,AF141394:0.00000):0.00323,AF141388:0.05424):0.01136,AF141581:0.00998):0.00323,(AF141444:0.00000,AF141446:0.00000):0.01673):0.02266,AF141478:0.02523):0.00162):0.02266,AF141595:0.05097):0.02108,((AF141393:0.00513,AF141458:0.00256):0.02621,AF141469:0.03147):0.02285):0.01541,AF141538:0.10137)Comamonadaceae:0.06378):0.01134,AF141449:0.06970)Burkholderiales:0.02283,(AF141600:0.07983,AF141482:0.04023):0.00573):0.01304,(((AF141486:0.00000,AF141403:0.00256):0.00000,AF141525:0.01542):0.01546,AF141461:0.01645)Methylophilales:0.06427):0.00326,(AF141551:0.02426,AF141507:0.04615)Neisseriales:0.09487)Betaproteobacteria:0.07947,AF141517:0.09902):0.05527,((AF141532:0.01292,AF141542:0.00257):0.00000,AF141424:0.00262)Ellin307/WD2124:0.12728):0.01546)Gamma_beta_proteobacteria:0.02854,((((((AF141557:0.00513,AF141480:0.00256)Roseobacter:0.13561,((AF141529:0.01265,AF141434:0.03671):0.02966,(((AF141526:0.00256,AF141495:0.00000):0.00000,AF141556:0.00513):0.00000,AF141548:0.00514):0.00164)Rhodobacter:0.04051)Rhodobacterales:0.22265,(AF141421:0.04514,AF141432:0.01539)Beijerinckiaceae:0.06527):0.00244,(((AF141544:0.00000,AF141545:0.00000):0.11111,AF141450:0.02767):0.10513,AF141396:0.07575):0.00653):0.00487,((AF141530:0.00000,AF141531:0.00000):0.02863,AF141448:0.04585)Sphingomonadales:0.13392):0.04010,(((((((((((AF141598:0.00323,AF141479:0.00836):0.00646,AF141593:0.01564):0.01291,AF141588:0.01033)Pelagibacter:0.00646,AF141583:0.00257):0.00000,AF141539:0.00000):0.01778,(AF141601:0.00162,AF141580:0.00769):0.02595):0.06022,(((AF141582:0.00000,AF141585:0.00000):0.02425,AF141590:0.02328):0.02274,AF141395:0.04863):0.00498)SAR11:0.04348,(AF141447:0.05616,(AF141594:0.00256,AF141584:0.00000):0.06659):0.02347):0.00663,AF141567:0.22932):0.01333,AF141589:0.05577):0.01247,AF141435:0.09271)Consistiales:0.05054)Alphaproteobacteria:0.08879):0.07221,((((AF141472:0.16321,(AF141537:0.01302,AF141496:0.00838)Desulfobulbaceae:0.12237):0.01639,AF141454:0.14646):0.01228,AF141504:0.15471):0.01396,AF141505:0.13406)'Deltaproteobacteria':0.00581)Proteobacteria:0.00332,AF141560:0.18908):0.00658,((((((AF141549:0.10860,((((((((AF141499:0.00258,AF141547:0.00258):0.00000,AF141559:0.00258):0.07190,AF141518:0.07560):0.04246,AF141474:0.07724)Cytophaga:0.01871,((((AF141515:0.02581,(AF141519:0.02952,AF141452:0.04379):0.01140):0.00484,AF141451:0.02101):0.01471,AF141466:0.06719)Sporocytophaga:0.02370,AF141524:0.04740):0.01239):0.01304,((AF141500:0.00256,AF141552:0.00000):0.00000,AF141502:0.00000):0.02059):0.04357,AF141407:0.14560):0.03134,AF141488:0.10821)Flavobacteriales:0.01232):0.02462,AF141543:0.10581):0.00743,AF141397:0.07576):0.04542,((AF141436:0.01862,AF141497:0.00325):0.00162,AF141460:0.00797)Flexibacteraceae:0.13515):0.00494,(AF141550:0.09051,AF141418:0.03264)Saprospiraceae:0.16732):0.01398,AF141514:0.25839)Bacteroidetes:0.13111):0.00907,((((AF141516:0.06810,AF141562:0.03493)'"Planctomycetacia"':0.09418,(AF141399:0.00528,AF141417:0.00262)'"Gemmatae"':0.21158)Planctomycetes:0.09116,((AF141391:0.02575,(AF141508:0.16401,(((AF141475:0.07537,AF141487:0.02921):0.02882,AF141513:0.08164):0.03217,AF141541:0.05472):0.02157)'Verrucomicrobiae(1)':0.00906):0.03990,((AF141455:0.00512,AF141387:0.00256):0.01132,(AF141408:0.00767,AF141406:0.00257):0.00815)'Opitutae(4)':0.18896)Verrucomicrobia:0.10524):0.01404,AF141536:0.11511):0.00906):0.00332,((((AF141463:0.01550,AF141476:0.01081)'Agrococcusetal.':0.09011,(((((AF141592:0.00664,AF141426:0.01026):0.02268,(AF141402:0.01525,(AF141484:0.00418,AF141445:0.00257):0.00486):0.00811):0.00486,AF141587:0.00170):0.02354,(AF141442:0.00674,AF141389:0.00000):0.01457):0.00489,AF141411:0.06787)Cellulomonadaceae:0.10377)Actinobacteridae:0.09674,((((AF141471:0.01091,AF141467:0.01136):0.01798,((AF141400:0.00767,AF141481:0.00000):0.05116,(AF141401:0.00512,AF141433:0.00256):0.03594):0.01066):0.01065,(((AF141522:0.00000,AF141423:0.01279):0.00000,AF141427:0.00000):0.00000,AF141420:0.00257):0.04350):0.09155,(((AF141440:0.00421,AF141520:0.00256):0.00161,AF141597:0.00256):0.00000,AF141485:0.00256)'BD2-10group':0.11726)Acidimicrobidae:0.04260)Actinobacteria:0.05311,((((((((AF141574:0.00000,AF141591:0.00256):0.00000,AF141579:0.00256):0.00000,AF141572:0.00769):0.00000,AF141565:0.00256):0.00000,AF141571:0.00513):0.00000,AF141575:0.00771):0.00324,((AF141596:0.00000,AF141564:0.00513):0.00162,(AF141599:0.00256,AF141578:0.01026):0.00162):0.00809)Prochlorales:0.18823,((((((AF141523:0.01285,AF141425:0.00258):0.00000,AF141429:0.00260):0.02350,AF141470:0.10407):0.00486,((AF141534:0.00809,AF141533:0.04316):0.00647,(AF141412:0.00000,AF141419:0.00257):0.00162):0.01456):0.04144,(((((AF141561:0.00372,AF141570:0.00000):0.00000,AF141576:0.00000):0.00000,AF141506:0.00743):0.00000,AF141573:0.00743):0.05520,AF141566:0.02577):0.01571)'Euglenaetal.chloroplasts':0.08248,AF141512:0.08318)Chloroplasts:0.03457)Cyanobacteria:0.16029):0.00661):0.01395,(AF141416:0.20329,AF141511:0.19093)'"Anaerolines"':0.09929):0.00986,AF141540:0.22389):0.03152,(AF141509:0.23520,(AF141390:0.25622,AF141483:0.15819)OP11-5:0.02369):0.08699):0.00913,(AF141415:0.00647,AF141535:0.00260)OP10:0.19926)Bacteria; 
"""

# Crump example envs
CRUMP_ENVS = """AF141399	R_FL	1
AF141411	R_FL	1
AF141408	R_FL	2
AF141403	R_FL	1
AF141410	R_FL	1
AF141398	R_FL	2
AF141391	R_FL	1
AF141389	R_FL	1
AF141395	R_FL	1
AF141401	R_FL	1
AF141390	R_FL	1
AF141393	R_FL	1
AF141396	R_FL	1
AF141402	R_FL	1
AF141407	R_FL	1
AF141387	R_FL	1
AF141394	R_FL	2
AF141409	R_FL	1
AF141400	R_FL	1
AF141397	R_FL	1
AF141405	R_FL	1
AF141388	R_FL	1
AF141424	R_PA	1
AF141421	R_PA	1
AF141433	R_PA	1
AF141428	R_PA	1
AF141432	R_PA	1
AF141426	R_PA	1
AF141430	R_PA	1
AF141413	R_PA	1
AF141419	R_PA	1
AF141423	R_PA	2
AF141429	R_PA	2
AF141422	R_PA	1
AF141431	R_PA	1
AF141415	R_PA	1
AF141418	R_PA	1
AF141416	R_PA	1
AF141420	R_PA	1
AF141417	R_PA	1
AF141434	R_PA	1
AF141412	R_PA	1
AF141414	R_PA	1
AF141463	E_FL	1
AF141493	E_FL	14
AF141459	E_FL	1
AF141461	E_FL	1
AF141447	E_FL	1
AF141479	E_FL	1
AF141449	E_FL	1
AF141465	E_FL	2
AF141435	E_FL	1
AF141457	E_FL	3
AF141468	E_FL	1
AF141487	E_FL	1
AF141472	E_FL	1
AF141466	E_FL	1
AF141444	E_FL	4
AF141473	E_FL	3
AF141439	E_FL	1
AF141436	E_FL	1
AF141455	E_FL	1
AF141443	E_FL	3
AF141483	E_FL	1
AF141476	E_FL	1
AF141441	E_FL	1
AF141440	E_FL	2
AF141474	E_FL	1
AF141486	E_FL	1
AF141481	E_FL	1
AF141480	E_FL	1
AF141470	E_FL	1
AF141458	E_FL	1
AF141460	E_FL	1
AF141478	E_FL	1
AF141450	E_FL	1
AF141471	E_FL	2
AF141442	E_FL	1
AF141454	E_FL	1
AF141488	E_FL	1
AF141451	E_FL	1
AF141456	E_FL	1
AF141452	E_FL	1
AF141482	E_FL	1
AF141448	E_FL	1
AF141484	E_FL	3
AF141475	E_FL	1
AF141548	E_PA	4
AF141520	E_PA	1
AF141508	E_PA	1
AF141523	E_PA	1
AF141547	E_PA	3
AF141530	E_PA	2
AF141510	E_PA	1
AF141513	E_PA	1
AF141524	E_PA	1
AF141516	E_PA	1
AF141543	E_PA	1
AF141503	E_PA	11
AF141498	E_PA	1
AF141532	E_PA	2
AF141549	E_PA	1
AF141501	E_PA	2
AF141525	E_PA	1
AF141534	E_PA	1
AF141506	E_PA	1
AF141519	E_PA	1
AF141550	E_PA	1
AF141496	E_PA	1
AF141535	E_PA	3
AF141512	E_PA	1
AF141514	E_PA	1
AF141537	E_PA	1
AF141533	E_PA	2
AF141511	E_PA	1
AF141539	E_PA	1
AF141545	E_PA	2
AF141500	E_PA	3
AF141505	E_PA	1
AF141497	E_PA	1
AF141541	E_PA	1
AF141518	E_PA	1
AF141551	E_PA	1
AF141504	E_PA	1
AF141546	E_PA	1
AF141515	E_PA	1
AF141529	E_PA	1
AF141507	E_PA	1
AF141540	E_PA	1
AF141521	E_PA	1
AF141538	E_PA	1
AF141522	E_PA	1
AF141509	E_PA	1
AF141517	E_PA	1
AF141536	E_PA	1
AF141557	O_UN	1
AF141569	O_UN	2
AF141561	O_UN	4
AF141578	O_UN	2
AF141567	O_UN	1
AF141566	O_UN	1
AF141568	O_UN	2
AF141579	O_UN	6
AF141562	O_UN	1
AF141559	O_UN	1
AF141560	O_UN	1
AF141577	O_UN	1
AF141583	O_FL	2
AF141595	O_FL	1
AF141599	O_FL	2
AF141582	O_FL	2
AF141600	O_FL	1
AF141594	O_FL	2
AF141592	O_FL	1
AF141597	O_FL	1
AF141590	O_FL	3
AF141581	O_FL	1
AF141587	O_FL	1
AF141588	O_FL	1
AF141586	O_FL	1
AF141591	O_FL	1
AF141589	O_FL	1
AF141593	O_FL	1
"""

def _load_tree(tree_str, envs_str, title):
    """Load tree and return """
    tree_in = StringIO(CRUMP_TREE)
    envs_in = StringIO(CRUMP_ENVS)
    print """\n(((((Running: %s""" % title
    return  tree_in, envs_in

def _display_results(out):
    print "Results:", out, "\n)))))"

if __name__ == '__main__':
    """
    Examples below using Crump tree and envs
    """
    from StringIO import StringIO

    tree_in, envs_in = _load_tree(CRUMP_TREE, CRUMP_ENVS, "unifrac - pairwise example (unweighted)")
    out = fast_unifrac_permutations_file(tree_in, envs_in, weighted=False, num_iters=100, verbose=True)
    _display_results(out)

    tree_in, envs_in = _load_tree(CRUMP_TREE, CRUMP_ENVS, "unifrac - pairwise example (weighted, normalized)")
    out = fast_unifrac_permutations_file(tree_in, envs_in, weighted='correct', num_iters=100, verbose=True)
    _display_results(out)

    tree_in, envs_in = _load_tree(CRUMP_TREE, CRUMP_ENVS, "unifrac - pairwise example (weighted)")
    out = fast_unifrac_permutations_file(tree_in, envs_in, weighted=True, num_iters=100, verbose=True)
    _display_results(out)

    tree_in, envs_in = _load_tree(CRUMP_TREE, CRUMP_ENVS, "unifrac - whole tree example")
    out = fast_unifrac_permutations_file(tree_in, envs_in, weighted=True, num_iters=100, verbose=True, test_on=TEST_ON_TREE) 
    _display_results(out)

    tree_in, envs_in = _load_tree(CRUMP_TREE, CRUMP_ENVS, "unifrac - each env example")
    out = fast_unifrac_permutations_file(tree_in, envs_in, weighted=True, num_iters=100, verbose=True, test_on=TEST_ON_ENVS) 
    _display_results(out)
    
    tree_in, envs_in = _load_tree(CRUMP_TREE, CRUMP_ENVS, "p test - pairwise example")
    out = fast_p_test_file(tree_in, envs_in, num_iters=10, verbose=True, test_on=TEST_ON_PAIRWISE)
    _display_results(out)

    tree_in, envs_in = _load_tree(CRUMP_TREE, CRUMP_ENVS, "p test - whole tree example")
    out = fast_p_test_file(tree_in, envs_in, num_iters=10, verbose=True, test_on=TEST_ON_TREE)
    _display_results(out)

    print "Done examples."
