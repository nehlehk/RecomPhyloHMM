from dendropy import Tree, DnaCharacterMatrix
import dendropy
import xml.etree.ElementTree as ET
import argparse



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
def give_taxon(tree,index):
    for node in tree.postorder_node_iter():
        if int(node.index) == index:
            return int(node.taxon.label)
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
    parser.add_argument('-t', "--treeFile", type=str, required= True, help='tree')
    parser.add_argument('-a', "--alignmentFile", type=str, required= True , help='fasta file')
    parser.add_argument('-x', "--xmlFile", type=str, help='xmlFile')
    args = parser.parse_args()


    tree_path = args.treeFile
    genomefile = args.alignmentFile
    xml_path = args.xmlFile
    alignment = dendropy.DnaCharacterMatrix.get(file=open(genomefile), schema="fasta")
    tips_num = len(alignment)
    tree = Tree.get_from_path(tree_path, 'newick')
    set_index(tree, alignment)

    make_beast_xml_original(tree, xml_path)