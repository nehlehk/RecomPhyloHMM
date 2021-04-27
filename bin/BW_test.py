import numpy as np
import numpy.linalg as la
import pandas as pd
import matplotlib.pyplot as plt
import hmmlearn.base
import hmmlearn._utils
from dendropy import Tree, DnaCharacterMatrix
import dendropy
import xml.etree.ElementTree as ET
from sklearn.metrics import mean_squared_error
import math
import argparse
import csv
import operator



class GTR_model:
    def __init__(self, rates, pi):
        self.rates = rates
        self.pi = pi
    #     ========================================================================
    def get_pi(self):
        return self.pi
    #     ========================================================================
    def p_matrix(self , br_length):
        p = np.zeros((4, 4))

        mu = 0
        freq = np.zeros((4, 4))
        q = np.zeros((4, 4))
        sqrtPi = np.zeros((4, 4))
        sqrtPiInv = np.zeros((4, 4))
        exchang = np.zeros((4, 4))
        s = np.zeros((4, 4))
        fun = np.zeros(1)
        a, b, c, d, e = self.rates
        f = 1

        freq = np.diag(self.pi)
        sqrtPi = np.diag(np.sqrt(self.pi))
        sqrtPiInv = np.diag(1.0 / np.sqrt(self.pi))
        mu = 1 / (2 * ((a * self.pi[0] * self.pi[1]) + (b * self.pi[0] * self.pi[2]) + (c * self.pi[0] * self.pi[3]) + (d * self.pi[1] * self.pi[2]) + (
                e * self.pi[1] * self.pi[3]) + (self.pi[2] * self.pi[3])))
        exchang[0][1] = exchang[1][0] = a
        exchang[0][2] = exchang[2][0] = b
        exchang[0][3] = exchang[3][0] = c
        exchang[1][2] = exchang[2][1] = d
        exchang[1][3] = exchang[3][1] = e
        exchang[2][3] = exchang[3][2] = f


        q = np.multiply(np.dot(exchang, freq), mu)

        for i in range(4):
            q[i][i] = -sum(q[i][0:4])


        s = np.dot(sqrtPi, np.dot(q, sqrtPiInv))


        eigval, eigvec = la.eig(s)
        eigvec_inv = la.inv(eigvec)

        left = np.dot(sqrtPiInv, eigvec)
        right = np.dot(eigvec_inv, sqrtPi)

        p = np.dot(left, np.dot(np.diag(np.exp(eigval * br_length)), right))


        return p
# **********************************************************************************************************************
def get_DNA_fromAlignment(alignment):

    alignment_len = alignment.sequence_size
    tips = len(alignment)
    column = []
    for l in range(alignment_len):
        col = ""
        for t in range(tips):
            col += str(alignment[t][l])
        column.append(col)

    return column
# **********************************************************************************************************************
def set_index(tree, dna):
    sequence_count = len(dna)

    for node in tree.postorder_node_iter():
        node.index = -1
        node.annotations.add_bound_attribute("index")

    s = sequence_count
    for node in tree.postorder_node_iter():
        if not node.is_leaf():
            node.index = s
            node.label = str(node.index)
            s += 1
        else:
            for idx, name in enumerate(dna):
                if str(name) == str(node.taxon):
                    node.index = idx
                    node.label = str(node.index)
                    break
# **********************************************************************************************************************
def set_tips_partial(tree, alignment):
    partial = np.zeros(((alignment_len, tips_num, 4)))
    for tip in range(tips_num):
      for site in range(alignment_len):
        dna = column[site]
        i = give_index(str(dna[tip]))
        partial[site][tip][i] = 1
    return partial
# **********************************************************************************************************************
def give_index(c):
    if c == "A":
        return 0
    elif c == "C":
        return 1
    elif c == "G":
        return 2
    elif c == "T":
        return 3
# **********************************************************************************************************************
def make_hmm_input_mixture(tree, alignment, tip_partial, model):
    sitell, partial = computelikelihood_mixture(tree, alignment, tip_partial, model)
    children = tree.seed_node.child_nodes()
    children_count = len(children)
    x = np.zeros((alignment.sequence_size, children_count * 4))
    for id, child in enumerate(children):
        # print(child.index)
        x[:, (id * 4):((id + 1) * 4)] = partial[:, child.index, :]
    return x
# **********************************************************************************************************************
def computelikelihood_mixture(tree, alignment ,tip_partial, model):
    alignment_len = alignment.sequence_size
    tips = len(dna)
    partial = np.zeros(((alignment_len,(2 * tips) -1, 4)))
    partial[:,0:tips,:] = tip_partial
    persite_ll = np.zeros(alignment_len)


    column = get_DNA_fromAlignment(alignment)

    uniqueCol = list(set(column))
    for u in range(len(uniqueCol)):
      indexes = []
      indexes = [id for id, x in enumerate(column) if x == uniqueCol[u]]
      site = indexes[0]
      for node in tree.postorder_node_iter():
          if not node.is_leaf():
              children = node.child_nodes()
              partial[site,node.index] = np.dot(model.p_matrix(children[0].edge_length), partial[site,children[0].index])
              for i in range(1, len(children)):
                  partial[site, node.index] *= np.dot(model.p_matrix(children[i].edge_length),partial[site,children[i].index])
      partial[indexes,:,:] = partial[site,:,:]
      p = np.dot(partial[site,tree.seed_node.index] , model.get_pi())
      persite_ll[indexes] = np.log(p)


    return persite_ll, partial
# **********************************************************************************************************************
def tree_evolver_rerooted(tree ,node ,nu):
    co_recom = nu/2
    if (node.edge_length is None):
       node.edge.length = 0
    node.edge.length = node.edge.length + co_recom
    recombination_tree = tree.as_string(schema="newick")

    return recombination_tree
# **********************************************************************************************************************
class phyloLL_HMM(hmmlearn.base._BaseHMM):
    def __init__(self, n_components, trees, model ,child_order,recom_child_order,n_iter=10, tol=1e-2):
        super().__init__(n_components)
        self.trees = trees
        self.model = model
        self.child_order = child_order
        self.recom_child_order = recom_child_order
        self.n_iter = n_iter
        self.tol = tol

    def _init(self, X, lengths):

        """Initializes model parameters prior to fitting.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Feature matrix of individual samples.
            n_samples == number of alignment sites
            n_features == 12: 4 site partials for each of 3 neighbour nodes

        lengths : array-like of integers, shape (n_sequences, )
            Lengths of the individual sequences in ``X``. The sum of
            these should be ``n_samples``.
        """
        init = 1. / self.n_components
        if 's' in self.init_params or not hasattr(self, "startprob_"):
            self.startprob_ = np.full(self.n_components, init)
        if 't' in self.init_params or not hasattr(self, "transmat_"):
            self.transmat_ = np.full((self.n_components, self.n_components), init)
        # n_fit_scalars_per_param = self._get_n_fit_scalars_per_param()
        # n_fit_scalars = sum(n_fit_scalars_per_param[p] for p in self.params)
        # if X.size < n_fit_scalars:
        #     _log.warning("Fitting a model with {} free scalar parameters with "
        #                  "only {} data points will result in a degenerate "
        #                  "solution.".format(n_fit_scalars, X.size))

    #     ==========================================================================
    def _check(self):
        """Validates model parameters prior to fitting.

        Raises
        ------

        ValueError
            If any of the parameters are invalid, e.g. if :attr:`startprob_`
            don't sum to 1.
        """
        self.startprob_ = np.asarray(self.startprob_)
        if len(self.startprob_) != self.n_components:
            raise ValueError("startprob_ must have length n_components")
        if not np.allclose(self.startprob_.sum(), 1.0):
            raise ValueError("startprob_ must sum to 1.0 (got {:.4f})"
                             .format(self.startprob_.sum()))

        self.transmat_ = np.asarray(self.transmat_)
        if self.transmat_.shape != (self.n_components, self.n_components):
            raise ValueError(
                "transmat_ must have shape (n_components, n_components)")
        if not np.allclose(self.transmat_.sum(axis=1), 1.0):
            raise ValueError("rows of transmat_ must sum to 1.0 (got {})"
                             .format(self.transmat_.sum(axis=1)))

    #     ==========================================================================
    def _compute_log_likelihood(self, X):
        """Computes per-component log probability under the model.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Feature matrix of individual samples.

        Returns
        -------
            Log probability of each sample in ``X`` for each of the
            model states.
        """

        return compute_logprob_phylo(X, self.trees, self.model, self.child_order, self.recom_child_order)

    def _initialize_sufficient_statistics(self):
        """Initializes sufficient statistics required for M-step.

        The method is *pure*, meaning that it doesn't change the state of
        the instance.  For extensibility computed statistics are stored
        in a dictionary.

        Returns
        -------
        nobs : int
            Number of samples in the data.

        start : array, shape (n_components, )
            An array where the i-th element corresponds to the posterior
            probability of the first sample being generated by the i-th
            state.

        trans : array, shape (n_components, n_components)
            An array where the (i, j)-th element corresponds to the
            posterior probability of transitioning between the i-th to j-th
            states.
        """
        n, dim = X.shape
        stats = {'nobs': n,
                 'start': np.asarray(self.startprob_),
                 'trans': np.asarray(self.transmat_)}
        return stats
# **********************************************************************************************************************
def compute_logprob_phylo(X, recom_trees, model,child_order,recom_child_order):
    n, dim = X.shape
    result = np.zeros((n, len(recom_trees)))
    for tree_id, item in enumerate(recom_trees):
        state_tree = dendropy.Tree.get(data=item, schema="newick")
        children = state_tree.seed_node.child_nodes()
        for site_id, partial in enumerate(X):
            order = child_order.index(recom_child_order[tree_id * len(children)])
            p = np.zeros(4)
            p = np.dot(model.p_matrix(children[0].edge_length), partial[order * 4:(order + 1) * 4])
            for i in range(1, len(children)):
                order = child_order.index(recom_child_order[(tree_id* len(children)) + i])
                p *= np.dot(model.p_matrix(children[i].edge_length), partial[order * 4:(order + 1) * 4])
            site_l = np.dot(p, model.get_pi())
            result[site_id, tree_id] = np.log(site_l)
    print(result.shape)
    print(result[0])
    return result
# **********************************************************************************************************************




if __name__ == "__main__":

    tree_path = '/home/nehleh/PhyloCode/Data/num_9/num_9_RAxML_bestTree.tree'
    tree = Tree.get_from_path(tree_path, 'newick')
    alignment = dendropy.DnaCharacterMatrix.get(file=open("/home/nehleh/PhyloCode/Data/num_9/num_9_wholegenome_9.fasta"),schema="fasta")
    recomLog = '/home/nehleh/PhyloCode/Data/num_9/num_9_BaciSim_Log.txt'




    tree = Tree.get_from_path(tree_path, 'newick')
    pi = [0.2184, 0.2606, 0.3265, 0.1946]
    rates = [0.975070, 4.088451, 0.991465, 0.640018, 3.840919]
    nu = 0.03
    mixtureProb = 0.9
    p_start = np.array([0.95, 0.05])
    p_trans = np.array([[0.997, 0.003],
                        [0.003, 0.997]])



    GTR_sample = GTR_model(rates, pi)
    column = get_DNA_fromAlignment(alignment)
    dna = column[0]
    tips_num = len(alignment)
    alignment_len = alignment.sequence_size
    nodes_number = len(tree.nodes())
    set_index(tree, alignment)



    tipdata = set_tips_partial(tree, alignment)
    recombination_trees = []
    child_order = []
    filter_fn = lambda n: hasattr(n, 'index') and n.index == 10
    target_node = tree.find_node(filter_fn=filter_fn)

    # ----------- Step 1 : Make input for hmm ------------------------------------------------------
    # --------------  Stetp 1.1 : re-root the tree based on the target node where the target node is each internal node of the tree.

    tree.reroot_at_node(target_node, update_bipartitions=False, suppress_unifurcations=True)
    recombination_trees.append(tree.as_string(schema="newick"))

    # --------------  Step 1.2: Calculate X based on this re-rooted tree
    X = make_hmm_input_mixture(tree, alignment, tipdata, GTR_sample)

    # ----------- Step 2: make recombination trees -----------------------------------------------
    temptree = {}
    recom_child_order = []
    for id, child in enumerate(target_node.child_node_iter()):  # to keep the order of clonal tree children
        recom_child_order.append(child.index)

    for id, child in enumerate(target_node.child_node_iter()):
        # print(id, child.index, child.taxon)
        temptree["tree{}".format(id)] = Tree.get_from_path(tree_path, 'newick')
        set_index(temptree["tree{}".format(id)], alignment)

        filter_fn = lambda n: hasattr(n, 'index') and n.index == target_node.index
        target_node_temp = temptree["tree{}".format(id)].find_node(filter_fn=filter_fn)
        temptree["tree{}".format(id)].reroot_at_node(target_node_temp, update_bipartitions=False,
                                                     suppress_unifurcations=True)

        kids = temptree["tree{}".format(id)].seed_node.child_nodes()  # to keep the order of recombination trees children
        # print(kids)
        for kid in kids:
            recom_child_order.append(kid.index)

        filter_fn = lambda n: hasattr(n, 'index') and n.index == child.index
        recombined_node = temptree["tree{}".format(id)].find_node(filter_fn=filter_fn)
        recombination_trees.append(tree_evolver_rerooted(temptree["tree{}".format(id)], recombined_node, nu))
        # hmm_trees.append(tree_evolver_rerooted(temptree["tree{}".format(id)], recombined_node, nu))
        child_order.append(recombined_node.index)



    model = phyloLL_HMM(n_components=2, trees=recombination_trees[0:2], model=GTR_sample, child_order=child_order,recom_child_order=recom_child_order)
    model.startprob_ = p_start
    model.transmat_ = p_trans

    posterior = model.predict_proba(X)
    hiddenStates = model.predict(X)
    score = model.score(X)

    print(hiddenStates[2000:2010])
    print("Score:",score)

    newmodel =  phyloLL_HMM(n_components=2, trees=recombination_trees[0:2], model=GTR_sample, child_order=child_order,recom_child_order=recom_child_order,n_iter = 100 )
    newmodel.startprob_ = p_start
    newmodel.transmat_ = p_trans
    print(newmodel)
    newmodel.fit(X)


    print(newmodel.startprob_)
    print(newmodel.transmat_)
    print(newmodel.score(X))

    print(newmodel.monitor_.converged)








