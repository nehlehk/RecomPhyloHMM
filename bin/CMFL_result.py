import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from dendropy import Tree
import dendropy
from sklearn.metrics import mean_squared_error
import csv


# **********************************************************************************************************************
def set_label(tree):
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            node.label = node.taxon.label
# **********************************************************************************************************************
def CFML_recombination(CFML_recomLog):
    CFMLData = np.zeros((alignment_len, tips_num))
    # CFMLData = np.zeros((alignment_len, nodes_number))
    df = pd.read_csv(CFML_recomLog, sep='\t', engine='python')
    # print(df)
    for i in range(len(df)):
        s = df['Beg'][i]
        e = df['End'][i]
        node = df['Node'][i]
        # if "NODE_" in str(node):
        #     node = node[5:]
        mynode = int(give_taxon_index(tree, node))
        CFMLData[s:e, mynode] = 1
        # CFMLData[s:e,int(node)] = 1

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
    # taxa = CFMLData.shape[1]
    # for i in range(taxa):
    for i in range(tips_num):
        # ax = fig.add_subplot(taxa, 1, i + 1)
        ax = fig.add_subplot(tips_num, 1, i + 1)
        # if i >= tips_num:
        #     node_label = str('NODE '+ str(i+1))
        #     desc = set()
        #     d = give_descendents_CFML(tree, node_label, desc)
        #     ax.plot(CFMLData[:, i], label=str(i+1) + ' is mrca:' + str(d), color=color[i % 5])
        # else:
        # if i < tips_num:
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
def real_recombination(recomLog):
    realData = np.zeros((alignment_len, tips_num))
    # realData = np.zeros((alignment_len, nodes_number))
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
def nodes_separation(nodes):
  nodes = str(nodes)
  nodes = nodes[1:-1]
  mynodes = nodes.split(",")
  return mynodes
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
def write_rmse_CFML(rmse_real_CFML):
    with open('rmse_CFML.csv', mode='w') as rmse_file:
        rmse_writer = csv.writer(rmse_file, delimiter=';', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        rmse_writer.writerow([rmse_real_CFML])
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''You did not specify any parameters.''')
    parser.add_argument('-t', "--treeFile", type=str, required=True, help='tree')
    parser.add_argument('-a', "--alignmentFile", type=str, required= True , help='fasta file')
    parser.add_argument('-c', "--cfmlFile", type=str, help='cfmlFile')
    parser.add_argument('-ct', "--cfmltreefile", type=str, help='cfmltreefile')
    parser.add_argument('-l', "--recomlogFile", type=str, help='recombination log file')
    args = parser.parse_args()

    tree_path = args.treeFile
    cfml_path = args.cfmlFile
    cfml_tree = args.cfmltreefile
    genomefile = args.alignmentFile
    recomLog = args.recomlogFile


    tree = Tree.get_from_path(tree_path, 'newick')
    cfml_tree = Tree.get_from_path(cfml_tree, 'newick')
    set_label(cfml_tree)
    alignment = dendropy.DnaCharacterMatrix.get(file=open(genomefile), schema="fasta")
    nodes_number = len(tree.nodes())
    tips_num = len(alignment)
    alignment_len = alignment.sequence_size
    set_index(tree, alignment)

    CFMLData = CFML_recombination(cfml_path)
    realData = real_recombination(recomLog)
    CFML_resultFig(cfml_tree, CFMLData)

    rmse_real_CFML = mean_squared_error(realData, CFMLData, squared=False)

    write_rmse_CFML(rmse_real_CFML)