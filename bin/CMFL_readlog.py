import pandas as pd
import numpy as np
from dendropy import Tree, DnaCharacterMatrix
import dendropy
import matplotlib.pyplot as plt


# **********************************************************************************************************************
def give_descendents_CFML(tree,node_label,result):
    if "NODE" in str(node_label):
        internal_recom_node = tree.find_node_with_label(node_label)
        children = internal_recom_node.child_nodes()
        for n in range(len(children)):
          r_node= children[n].label
          if "NODE" in str(r_node):
            give_descendents_CFML(tree,r_node,result)
          else:
            result.add(r_node)
    return result
# **********************************************************************************************************************
def CFML_resultFig(tree,CFMLData):
    fig = plt.figure(figsize=(tips_num + 9, tips_num / 2))
    color = ['red', 'green', 'purple', 'blue', 'black']
    taxa = CFMLData.shape[1]
    for i in range(taxa):
        ax = fig.add_subplot(taxa, 1, i + 1)
        if i >= tips_num:
            node_label = str('NODE '+ str(i+1))
            desc = set()
            d = give_descendents_CFML(tree, node_label, desc)
            ax.plot(CFMLData[:, i], label=str(i+1) + ' is mrca:' + str(d), color=color[i % 5])
        else:
            ax.plot(CFMLData[:, i], label=i, color=color[i % 5])
        ax.legend(bbox_to_anchor=(0.045, 1.5), prop={'size': 10})
        # ax.plot(CFMLData[:, i],label=i, color=color[i % 5])
        # ax.legend(bbox_to_anchor=(0.04, 1.33) ,prop={'size':10} )
        ax.set_frame_on(False)
        ax.axis('off')
    ax.axis('on')
    ax.set_yticklabels([])
    # plt.show()
    plt.savefig("CFML_Recombination.jpeg")
# **********************************************************************************************************************
def set_label(tree):
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            node.label = node.taxon.label
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
def give_taxon_index(tree,taxa):
    node_mapping = np.zeros((tips_num,2))
    i = 0
    for node in tree.postorder_node_iter():
      if node.is_leaf():
        node_mapping[i][0] = int(str(node.taxon.label))
        node_mapping[i][1] = int(str(node.index))
        i = i+1

    for i in range(len(node_mapping)):
      if node_mapping[i][0] == taxa:
        return node_mapping[i][1]
# **********************************************************************************************************************
def CFML_recombination(CFML_recomLog):
    CFMLData = np.zeros((alignment_len, nodes_number))
    df = pd.read_csv(CFML_recomLog, sep='\t', engine='python')
    # print(df)
    for i in range(len(df)):
        s = df['Beg'][i]
        e = df['End'][i]
        node = df['Node'][i]
        if "NODE_" in str(node):
            node = node[5:]
        CFMLData[s:e,int(node)] = 1

    return CFMLData
# **********************************************************************************************************************


tree_path = '/home/nehleh/Documents/num_2/num_1_CFML.labelled_tree.newick'
tree = Tree.get_from_path(tree_path, 'newick')
alignment = dendropy.DnaCharacterMatrix.get(file=open("/home/nehleh/Documents/num_2/num_1_wholegenome.fasta"),schema="fasta")

nodes_number = len(tree.nodes())
tips_num = len(alignment)
alignment_len = alignment.sequence_size

set_label(tree)

# desc = set()
# d = give_descendents_CFML(tree, 'NODE 9', desc)
# print(d)

CFMLData= CFML_recombination('/home/nehleh/Documents/num_2/num_1_CFML.importation_status.txt')

CFML_resultFig(tree,CFMLData)

# print(CFMLData[4200:4205])