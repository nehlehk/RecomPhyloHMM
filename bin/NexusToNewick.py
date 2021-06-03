import dendropy
import argparse

parser = argparse.ArgumentParser(description='''You did not specify any parameters.''')
parser.add_argument('-t', "--treeFile", type=str,  help='tree')
parser.add_argument('-o', "--outputtree", type=str,  help='tree')

args = parser.parse_args()
tree_path = args.treeFile
outputtree = args.outputtree


# tree_path = '/home/nehleh/PhyloCode/Result/10nodes-20sample-len10000/short_philobacteria/num_1/num_1_beastSeqTree.nexus'
# outputtree = '/home/nehleh/PhyloCode/Result/10nodes-20sample-len10000/short_philobacteria/num_1/TEST_beastSeqTree.newick'

tree = dendropy.Tree.get_from_path(tree_path, 'nexus')
tree.write(path = outputtree, schema="newick")