import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from dendropy import Tree
import dendropy


# **********************************************************************************************************************
def set_label(tree):
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            node.label = node.taxon.label
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''You did not specify any parameters.''')
    parser.add_argument('-t', "--treeFile", type=str, required=True, help='tree')
    parser.add_argument('-a', "--alignmentFile", type=str, required= True , help='fasta file')
    parser.add_argument('-c', "--cfmlFile", type=str, help='cfmlFile')
    parser.add_argument('-ct', "--cfmltreefile", type=str, help='cfmltreefile')
    args = parser.parse_args()

    tree_path = args.treeFile
    cfml_path = args.cfmlFile
    cfml_tree = args.cfmltreefile
    genomefile = args.alignmentFile


    tree = Tree.get_from_path(tree_path, 'newick')
    cfml_tree = Tree.get_from_path(cfml_tree, 'newick')
    set_label(cfml_tree)
    alignment = dendropy.DnaCharacterMatrix.get(file=open(genomefile), schema="fasta")
    nodes_number = len(tree.nodes())
    tips_num = len(alignment)
    alignment_len = alignment.sequence_size

    CFMLData = CFML_recombination(cfml_path)
    CFML_resultFig(cfml_tree, CFMLData)