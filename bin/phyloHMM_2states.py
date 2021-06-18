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
import argparse
import csv
import operator
import itertools




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
class phyloLL_HMM(hmmlearn.base._BaseHMM):
    def __init__(self, n_components, trees, model ,child_order,recom_child_order):
        super().__init__(n_components)
        self.trees = trees
        self.model = model
        self.child_order = child_order
        self.recom_child_order = recom_child_order

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
        n_fit_scalars_per_param = self._get_n_fit_scalars_per_param()
        n_fit_scalars = sum(n_fit_scalars_per_param[p] for p in self.params)
        if X.size < n_fit_scalars:
            _log.warning("Fitting a model with {} free scalar parameters with "
                         "only {} data points will result in a degenerate "
                         "solution.".format(n_fit_scalars, X.size))

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
def set_tips_partial(tree, alignment):
    partial = np.zeros(((alignment_len, tips_num, 4)))
    for tip in range(tips_num):
      for site in range(alignment_len):
        dna = column[site]
        i = give_index(str(dna[tip]))
        partial[site][tip][i] = 1
    return partial
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
def update_mixture_partial(alignment,tree,node,tipdata,posterior):
  column = get_DNA_fromAlignment(alignment)
  for site in range(alignment_len):
    dna = column[site]
    my_number = give_index(dna[node.index])
    rho = give_rho(node,posterior,site)
    for i in range(4):
      if i == my_number:
        tipdata[site,node.index,i] = 1
      else:
        tipdata[site,node.index,i] = rho

  return tipdata
# **********************************************************************************************************************
def give_rho(node,recom_prob,site):
  parent = node.parent_node
  if parent == tree.seed_node:
    myindex = parent.index -1
  else:
    myindex = parent.index
  rho = recom_prob[site][1]
  return rho
# **********************************************************************************************************************
def tree_evolver_rerooted(tree ,node ,nu):
    co_recom = nu/2
    if (node.edge_length is None):
       node.edge.length = 0
    node.edge.length = node.edge.length + co_recom
    recombination_tree = tree.as_string(schema="newick")

    return recombination_tree
# **********************************************************************************************************************
def give_taxon(tree,index):
    for node in tree.postorder_node_iter():
        if int(node.index) == index:
            return int(node.taxon.label)
# **********************************************************************************************************************
def compute_logprob_phylo_fast(X, recom_trees, model):
    n, dim = X.shape
    result = np.zeros((n, len(recom_trees)))
    for tree_id, item in enumerate(recom_trees):
        state_tree = dendropy.Tree.get(data=item, schema="newick")
        children = state_tree.seed_node.child_nodes()

        uniqueX = np.unique(X, axis=0)
        for uniq_id, partial in enumerate(uniqueX):
            indexes = []
            indexes = [id for id, x in enumerate(X) if (x == uniqueX[uniq_id]).all()]
            p = np.zeros(4)
            p = np.dot(model.p_matrix(children[0].edge_length), partial[0:4])
            for i in range(1, len(children)):
                p *= np.dot(model.p_matrix(children[i].edge_length), partial[i * 4:(i + 1) * 4])
            site_l = np.dot(p, model.get_pi())
            result[indexes, tree_id] = np.log(site_l)
    return result
# **********************************************************************************************************************
def compute_logprob_phylo(X,recom_trees,model,child_order,X_child_order):
    n, dim = X.shape
    result = np.zeros((n, len(recom_trees)))
    for tree_id, item in enumerate(recom_trees):
        state_tree = dendropy.Tree.get(data=item, schema="newick")
        children = state_tree.seed_node.child_nodes()
        # print(children)
        for site_id, partial in enumerate(X):
            # order = child_order.index(X_child_order[tree_id * len(children)])
            order = X_child_order.index(child_order[0])
            p = np.zeros(4)
            p = np.dot(model.p_matrix(children[0].edge_length), partial[order * 4:(order + 1) * 4])
            for i in range(1, len(children)):
                # order = child_order.index(X_child_order[(tree_id* len(children)) + i])
                order = X_child_order.index(child_order[i])
                p *= np.dot(model.p_matrix(children[i].edge_length), partial[order * 4:(order + 1) * 4])
            site_l = np.dot(p, model.get_pi())
            result[site_id, tree_id] = np.log(site_l)
    return result
# **********************************************************************************************************************
def recom_maker(r_tree,index,nu):
    filter_fn = lambda n: hasattr(n, 'index') and n.index == index
    recombined_node = r_tree.find_node(filter_fn=filter_fn)
    return tree_evolver_rerooted(r_tree, recombined_node, nu)
# **********************************************************************************************************************
def phylohmm(tree,alignment,nu,p_start,p_trans,threshold):
    mytree = []
    posterior = []
    hiddenStates = []
    score = []
    tipdata = set_tips_partial(tree, alignment)
    r_node = []
    t_node = []

    print(tree.as_ascii_plot(show_internal_node_labels=True))

    for id_tree, target_node in enumerate(tree.postorder_internal_node_iter(exclude_seed_node=True)):
      # if target_node.index < 15:
        print(target_node.index)
        recombination_trees = []
        recombination_nodes = []
        child_order = []
        mytree.append(Tree.get_from_path(tree_path, 'newick'))
        set_index(mytree[id_tree], alignment)

        # ----------- Step 1 : Make input for hmm ------------------------------------------------------
        # --------------  Stetp 1.1 : re-root the tree based on the target node where the target node is each internal node of the tree.

        mytree[id_tree].reroot_at_node(target_node, update_bipartitions=False, suppress_unifurcations=True)
        recombination_trees.append(mytree[id_tree].as_string(schema="newick"))

        # --------------  Step 1.2: Calculate X based on this re-rooted tree

        X = make_hmm_input_mixture(mytree[id_tree], alignment, tipdata, GTR_sample)
        # to keep the order of clonal tree children
        X_child_order = []
        for id, child in enumerate(target_node.child_node_iter()):
            X_child_order.append(child.index)

        # ----------- Step 2: make recombination trees -----------------------------------------------

        temptree = Tree.get_from_path(tree_path, 'newick')
        set_index(temptree, alignment)

        filter_fn = lambda n: hasattr(n, 'index') and n.index == target_node.index
        target_node_temp = temptree.find_node(filter_fn=filter_fn)
        temptree.reroot_at_node(target_node_temp, update_bipartitions=False,suppress_unifurcations=True)
        kids = temptree.seed_node.child_nodes()

        for k1, kid1 in enumerate(kids):
            child_order.append(kid1.index) #keep the order of children after reroot
          # if (kid1.index < tips_num) or (target_node.index < kid1.index):  # to avoid having repeated internal nodes
            temptree = Tree.get_from_path(tree_path, 'newick')
            set_index(temptree, alignment)
            filter_fn = lambda n: hasattr(n, 'index') and n.index == target_node.index
            target_node_temp = temptree.find_node(filter_fn=filter_fn)
            temptree.reroot_at_node(target_node_temp, update_bipartitions=False,suppress_unifurcations=True)
            recombination_trees.append(recom_maker(temptree, kid1.index, nu))
            recombination_nodes.append(kid1)

        # ----------- Step 3: Call phyloHMM -----------------------------------------------------
        for h in range(1,len(recombination_trees)):
            model = phyloLL_HMM(n_components=2,trees=[recombination_trees[0],recombination_trees[h]],model=GTR_sample,child_order=child_order,recom_child_order=X_child_order)
            model.startprob_ = p_start
            model.transmat_ = p_trans

            p = model.predict_proba(X)
            hidden = model.predict(X)

            posterior.append(p)
            hiddenStates.append(hidden)
            score.append(model.score(X))

            r_node.append(recombination_nodes[h-1].index)
            t_node.append(target_node.index)

            tree_updatePartial = Tree.get_from_path(tree_path, 'newick')
            set_index(tree_updatePartial, alignment)
            filter_fn = lambda n: hasattr(n, 'index') and n.index == X_child_order[h-1]
            update_child = tree_updatePartial.find_node(filter_fn=filter_fn)
            # print("update_child:",update_child)
            if update_child.is_leaf():
                update_mixture_partial(alignment, tree_updatePartial, update_child, tipdata, p)


    recom_prob = pd.DataFrame({'recom_nodes':r_node, 'target_node':t_node,  'posterior':posterior})
    return tipdata,posterior,hiddenStates,score,recom_prob,r_node,t_node
# **********************************************************************************************************************
def give_descendents(tree,node_index,result):
  if node_index >= tips_num:
    internal_recom_node = tree.find_node_with_label(str(node_index))
    children = internal_recom_node.child_nodes()
    for n in range(len(children)):
      r_node= int(children[n].index)
      if r_node >= tips_num:
        give_descendents(tree,r_node,result)
      else:
        result.add(give_taxon(tree,r_node))
  return result
# **********************************************************************************************************************
def recom_resultFig_dm(recom_prob,mixtureProb):
    output = np.zeros((alignment_len, nodes_number))
    for i in range(len(recom_prob)):
        if (recom_prob['recom_nodes'][i] < tips_num):
            for j in range(alignment_len):
                if (recom_prob['posterior'][i][j][1] >= mixtureProb):
                    output[j, recom_prob['recom_nodes'][i]] = 1
        else:
            # for j in range(alignment_len):
            #     if (recom_prob['posterior'][i][j][1] >= mixtureProb):
            #         output[j, recom_prob['target_node'][i]] = 1
            for j in range(i+1,len(recom_prob)):
                if (recom_prob['recom_nodes'][i] == recom_prob['target_node'][j]) and (recom_prob['recom_nodes'][j] == recom_prob['target_node'][i]) :
                  for k in range(alignment_len):
                      if ((recom_prob['posterior'][i][k][1] >= mixtureProb) and (recom_prob['posterior'][j][k][1] >= mixtureProb)):
                          output[k, recom_prob['target_node'][i]] = 1
                        # if (recom_prob['posterior'][i][k] < recom_prob['posterior'][j][k]):
                        #   recom_prob['posterior'][i][k] = recom_prob['posterior'][j][k]
                        # if (recom_prob['posterior'][i][k] >= mixtureProb):
                        #     output[k, recom_prob['target_node'][i]] = 1

    fig = plt.figure(figsize=(tips_num + 9, tips_num / 2))
    color = ['red', 'green', 'purple', 'blue', 'black']
    clonaltree = Tree.get_from_path(tree_path, 'newick')
    set_index(clonaltree, alignment)
    for i in range(nodes_number):
        ax = fig.add_subplot(nodes_number, 1, i + 1)
        if i >= tips_num:
            desc = set()
            d = give_descendents(clonaltree, i, desc)
            ax.plot(output[:, i], label=str(i) + ' is mrca:' + str(d), color=color[i % 5])
        else:
            ax.plot(output[:, i], label=give_taxon(clonaltree, i), color=color[i % 5])
        ax.legend(bbox_to_anchor=(0.045, 1.5), prop={'size': 10})
        ax.set_frame_on(False)
        ax.axis('off')

    ax.axis('on')
    ax.set_yticklabels([])
    plt.savefig("PhyloHMM_Recombination_two.jpeg")
    # plt.show()

    return output
# **********************************************************************************************************************
def internal_plot(c_tree,posterior,hiddenStates,score,r_node,t_node):
    for i in range(len(posterior)):
        poster = posterior[i]
        hidden = hiddenStates[i]
        sc = score[i]

        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(2, 1, 1)
        ax.plot(hidden)
        ax.set_title("Node"+str(t_node[i]))
        ax.set_ylabel("Clonal - Recombination State")
        ax2 = fig.add_subplot(2, 1, 2)
        ax2.plot(poster[:, 0], label="ClonalFrame")
        if r_node[i] < tips_num:
            label = give_taxon(c_tree, r_node[i])
        else:
            label = r_node[i]
        ax2.plot(poster[:, 1], label="Node"+str(label))
        ax2.set_ylabel("posterior probability for each state")
        ax2.legend(loc=1, bbox_to_anchor=(1.13, 1.1))
        plt.savefig("posterior node"+str(t_node[i])+"_node"+str(label)+".jpeg")
        # plt.show()
# **********************************************************************************************************************
def make_beast_xml_partial(tipdata,tree,xml_path):
    my_tipdata = tipdata.transpose(1, 0, 2)
    my_xml = ET.parse(xml_path)
    root = my_xml.getroot()
    data = root.find("data")

    for i in range(my_tipdata.shape[0]):
        x = ''
        c = ET.Element("sequence")
        c.set("taxon" , str(give_taxon(tree,i)))
        c.set("uncertain" , "true")
        for j in range(my_tipdata.shape[1]):
          x = x + str(repr(my_tipdata[i,j,:]))[7:-2] + ';'
        c.text = '\n' + x +'\n'
        data.insert(i,c)
        c.tail = "\n"

    my_xml.write('RecomPartial_two.xml' ,encoding="utf-8", xml_declaration=True)
# **********************************************************************************************************************
def make_beast_xml_original(tree,xml_path):
    my_xml = ET.parse(xml_path)
    root = my_xml.getroot()
    data = root.find("data")

    for i in range(tips_num):
        x = ''
        c = ET.Element("sequence")
        c.set("taxon" , str(give_taxon(tree,i)))
        c.text = '\n' + str(alignment[i]) + '\n'
        data.insert(i,c)
        c.tail = "\n"

    my_xml.write('originalSeq.xml' ,encoding="utf-8", xml_declaration=True)
# **********************************************************************************************************************
def give_taxon_index(tree,taxa):
    node_mapping = np.zeros((nodes_number,2))
    i = 0
    for node in tree.postorder_node_iter():
      if node.is_leaf():
        node_mapping[i][0] = int(str(node.taxon.label))
      else:
        node_mapping[i][0] = int(str(node.index))
      node_mapping[i][1] = int(str(node.index))
      i = i+1

    for i in range(len(node_mapping)):
      if int(node_mapping[i][0]) == int(taxa):
        return node_mapping[i][1]
# **********************************************************************************************************************
def nodes_separation(nodes):
  nodes = str(nodes)
  nodes = nodes[1:-1]
  mynodes = nodes.split(",")
  return mynodes
# **********************************************************************************************************************
def real_recombination(recomLog):
    realData = np.zeros((alignment_len, nodes_number))
    df = pd.read_csv(recomLog,sep='\t', engine='python')
    # print(df)
    recom = df.loc[df['status'] != 'clonal']
    recom = recom.reset_index(drop=True)
    for i in range(len(recom)):
        s = recom['start'][i]
        t = recom['end'][i]
        # print(recom['nodes'][i])
        nodes = nodes_separation(recom['nodes'][i])
        # print(nodes)
        for i in range(len(nodes)):
            mynode = int(give_taxon_index(tree, nodes[i]))
            realData[s:t, mynode] = 1

    return realData
# **********************************************************************************************************************
def comparison_plot(RealData,predictionData):
    fig = plt.figure(figsize=(15, 25))

    for i in range(tips_num):
        ax = fig.add_subplot(tips_num, 1, i + 1)
        ax.plot(RealData[:, i], label="RealData")
        ax.plot(predictionData[:, i], label="predictionData ")
        ax.set_ylabel("posterior probability for each state")
        ax.legend(loc='best')
        plt.savefig("taxa" + str(i) + ".jpeg")
# **********************************************************************************************************************
def write_rmse_phylohmm(nu,mixtureProb,rmse_real_phylohmm):
    with open('rmse_phylohmm_two.csv', mode='w') as rmse_file:
        rmse_writer = csv.writer(rmse_file, delimiter=';', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        rmse_writer.writerow([nu,mixtureProb,rmse_real_phylohmm])
# **********************************************************************************************************************
def set_label(tree):
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            node.label = node.taxon.label
# **********************************************************************************************************************
def phyloHMM_Log(c_tree,output):
    nodes = []
    starts = []
    ends = []
    recomlens = []

    for i in range(output.shape[1]):
        non_zeros = [[i for i, value in it] for key, it in itertools.groupby(enumerate(output[:, i]), key=operator.itemgetter(1)) if key != 0]
        for j in range(len(non_zeros)):
            if i < tips_num:
                n = give_taxon(c_tree, i)
            else:
                n = i
            nodes.append(n)
            starts.append(non_zeros[j][0])
            ends.append(non_zeros[j][len(non_zeros[j]) - 1])
            recomlens.append(non_zeros[j][len(non_zeros[j]) - 1] - non_zeros[j][0])

    all_data = {'nodes': nodes, 'start': starts, 'end': ends, 'len': recomlens}
    df = pd.DataFrame(all_data)
    df = df.sort_values(by=['nodes'], ascending=[True])
    df.to_csv('Recom_phyloHMM_Log_two.txt', sep='\t', header=True , index = False)

    return df
# **********************************************************************************************************************
if __name__ == "__main__":

    tree_path = '/home/nehleh/Desktop/sisters/new_BaciSim/num_1_RAxML_bestTree.tree'
    # recomLog = '/home/nehleh/work/results/num_5/num_5_BaciSim_Log.txt'
    genomefile = '/home/nehleh/Desktop/sisters/new_BaciSim/num_1_wholegenome_1.fasta'
    # xml_path = '/home/nehleh/PhyloCode/RecomPhyloHMM/bin/GTR_template.xml'

    parser = argparse.ArgumentParser(description='''You did not specify any parameters.''')
    # parser.add_argument('-t', "--treeFile", type=str, required= True, help='tree')
    # parser.add_argument('-a', "--alignmentFile", type=str, required= True , help='fasta file')
    # parser.add_argument('-l', "--recomlogFile", type=str, help='recombination log file')
    parser.add_argument('-nu', "--nuHmm", type=float,default=0.03,help='nuHmm')
    parser.add_argument('-p', "--mixtureProb", type=float, default=0.9, help='mixtureProb')
    parser.add_argument('-f', "--frequencies", type=list, default= [0.2184,0.2606,0.3265,0.1946],help='frequencies')
    parser.add_argument('-r', "--rates", type=list, default= [0.975070 ,4.088451 ,0.991465 ,0.640018 ,3.840919 ], help='rates')
    parser.add_argument('-s', "--startProb", type=list, default= [0.85, 0.05, 0.05, 0.05],help='frequencies')
    parser.add_argument('-m', "--transmat", type=list, default= [[0.997, 0.001, 0.001, 0.001],[0.001, 0.997, 0.001, 0.001],[0.001, 0.001, 0.997, 0.001],[0.001, 0.001, 0.001, 0.997]], help='rates')
    # parser.add_argument('-x', "--xmlFile", type=str, help='xmlFile')
    parser.add_argument('-st', "--status", type=int, default=1, help='1 for simulated data, 0 is for real dataset')
    args = parser.parse_args()

    # tree_path = args.treeFile
    # genomefile = args.alignmentFile
    # recomLog = args.recomlogFile
    # xml_path = args.xmlFile
    pi = args.frequencies
    rates = args.rates
    nu = args.nuHmm
    p_start = args.startProb
    p_trans = args.transmat
    mixtureProb = args.mixtureProb
    status = args.status



    # ============================================ setting parameters ==================================================
    tree = Tree.get_from_path(tree_path, 'newick')
    alignment = dendropy.DnaCharacterMatrix.get(file=open(genomefile), schema="fasta")
    set_index(tree, alignment)

    GTR_sample = GTR_model(rates, pi)
    column = get_DNA_fromAlignment(alignment)
    dna = column[0]
    tips_num = len(alignment)
    alignment_len = alignment.sequence_size
    nodes_number = len(tree.nodes())

    p_start = np.array([0.99, 0.01])
    p_trans = np.array([[0.9999, 0.0001],
                        [0.0001, 0.9999]])


    # ============================================= main operation =====================================================

    tipdata,posterior,hiddenStates,score,recom_prob,r_node,t_node = phylohmm(tree, alignment, nu , p_start , p_trans, mixtureProb)

    c_tree = Tree.get_from_path(tree_path, 'newick')
    set_index(c_tree, alignment)

    # print(recom_prob)
    internal_plot(c_tree,posterior, hiddenStates, score,r_node,t_node)
    phyloHMMData = recom_resultFig_dm(recom_prob, mixtureProb)
    phyloHMM_log = phyloHMM_Log(c_tree, phyloHMMData)

    # ======================================= providing xml files for beast ============================================
    # make_beast_xml_partial(tipdata,c_tree,xml_path)
    # make_beast_xml_original(c_tree,xml_path)

    # ============================================    only for simulated data  =========================================
    # if status:
    #     realData = real_recombination(recomLog)
    #     rmse_real_phyloHMM= mean_squared_error(realData,phyloHMMData,squared=False)
    #     write_rmse_phylohmm(nu,mixtureProb,rmse_real_phyloHMM)