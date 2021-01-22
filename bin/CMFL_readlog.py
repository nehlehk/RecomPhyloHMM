import pandas as pd
import numpy as np
from dendropy import Tree, DnaCharacterMatrix
import dendropy

alignment_len = 5000
tips_num = 10


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
    CFMLData = np.zeros((alignment_len, tips_num))
    df = pd.read_csv(CFML_recomLog,sep='\t', engine='python')
    for i in range(len(df)):
        s = df['Beg'][i]
        e = df['End'][i]
        node = df['Node'][i]
        CFMLData[s:e,int(give_taxon_index(tree, node))] = 1

    return CFMLData
# **********************************************************************************************************************


tree_path = '/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/BaciSim/2/RAxML_bestTree.tree'
tree = Tree.get_from_path(tree_path, 'newick')
alignment = dendropy.DnaCharacterMatrix.get(file=open("/home/nehleh/Documents/0_Research/PhD/Data/simulationdata/BaciSim/2/wholegenome.fasta"),schema="fasta")
set_index(tree, alignment)

CFMLData= CFML_recombination('/home/nehleh/CF_ML.importation_status.txt')

print(CFMLData[4700:4706])