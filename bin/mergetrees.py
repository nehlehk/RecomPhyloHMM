import numpy as np
import dendropy
from dendropy.simulate import treesim
import random
from dendropy import Tree
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
import argparse




def regraft_recomnodes(tree,attached_node,attached_id,prunenode,recomnode,ancestor):
  relocated_nodes = ancestor[attached_id-1]  # relocated node is the adjacent node of recombinant node
  tree.prune_subtree(prunenode) # the original recombinant node was removed to reinsert in the other side
  attached_node.remove_child(relocated_nodes) # relocated node was removed to reinsert in to right side
  newborn = dendropy.datamodel.treemodel.Node()  # newborn in the new mrca of recombinant node and its sister
  newborn.add_child(prunenode)
  prunenode.edge_length = recomnode.edge_length
  if attached_node.edge_length is None:
    attached_node.edge_length = 0
  newborn.edge_length = tree.max_distance_from_root() - (recomnode.edge_length + attached_node.edge_length)
  relocated_nodes.edge_length = relocated_nodes.edge_length - newborn.edge_length
  newborn.add_child(relocated_nodes)
  attached_node.add_child(newborn)
# ----------------------------------------------------------------------------------------------------------------------
def find_attached_node(tree,node,prunenode):
  if ((node.edge_length) < tree.max_distance_from_root()) :
    ancestor = []
    for id,tmp_node in enumerate(prunenode.ancestor_iter()):
        ancestor.append(tmp_node)
        if node.edge_length < tmp_node.distance_from_tip() :
            attached_node = tmp_node
            attached_id = id
            break
  return attached_node,attached_id,ancestor
# ----------------------------------------------------------------------------------------------------------------------
def common_nodes(clonaltree,overlap_nodes):
  kids = []
  nodes = []
  s_overlap_nodes = np.sort(overlap_nodes)
  for i in range(len(s_overlap_nodes)):
    nodes.append(s_overlap_nodes[i])
    if s_overlap_nodes[i] >= tips_number:
      desc = set()
      kids.append(give_descendents(clonaltree,s_overlap_nodes[i],desc))
    else:
      kids.append({(s_overlap_nodes[i])})

  subsets = []
  for i in range(len(kids)):
    for j in range(i+1,len(kids)):
      if (kids[i].issubset(kids[j])):
        # print(i,j)
        subsets.append([nodes[i], nodes[j]])

  return subsets
# ----------------------------------------------------------------------------------------------------------------------
def set_index(tree):
    for node in tree.postorder_node_iter():
      node.index = -1
      node.annotations.add_bound_attribute("index")

    s = len(tree.leaf_nodes())
    for node in tree.postorder_node_iter():
      if not node.is_leaf():
          node.index = s
          node.label = str(node.index)
          s += 1
      else:
          node.index = node.taxon.label
          node.label = str(node.index)
# ----------------------------------------------------------------------------------------------------------------------
def give_equivalent_node(recomtree):
  e = []
  for edge in recomtree.postorder_edge_iter():
    if edge.length is None:
      edge.length = 0
    e.append(edge.length)

  m = np.max(e)
  for node in recomtree.postorder_node_iter():
    if node.edge_length == m and node.is_leaf():
      return node.label,m
    elif node.edge_length == m and node.is_internal():
      return node.label,node.distance_from_tip() + m
# ----------------------------------------------------------------------------------------------------------------------
def give_older_recom(s_equ):
  count = Counter(a[0] for a in s_equ)
  max = np.max([a[2] for a in s_equ])
  [a for a in s_equ if (count[a[0]] == 1) or (a[2]==max)]
  return [a for a in s_equ if (count[a[0]] == 1) or (a[2]==max)]
# ----------------------------------------------------------------------------------------------------------------------
def new_mrca(tree1,tree2,primenode,deletenode):
  a = []
  desc = set()
  d = give_descendents(tree1,primenode,desc)
  rlist = set()
  removelist = give_descendents(tree1,deletenode,rlist)
  if len(removelist)> 0 :
    a = [x for x in d if x not in removelist]
  else:
    a = list(d)
  final_node = my_mrca(tree2,a)
  return final_node
# ----------------------------------------------------------------------------------------------------------------------
def give_descendents(tree,node_index,result):
  if node_index >= tips_number:
    internal_recom_node = tree.find_node_with_label(str(node_index))
    children = internal_recom_node.child_nodes()
    for n in range(len(children)):
      r_node= int(children[n].index)
      if r_node >= tips_number:
        give_descendents(tree,r_node,result)
      else:
        result.add(r_node)
  return result
# ----------------------------------------------------------------------------------------------------------------------
def my_mrca(tree,tips):
  pdm = tree.phylogenetic_distance_matrix()
  taxon = tree.taxon_namespace
  node0 = [i for i,x in enumerate(taxon) if x.label==str(tips[0])]
  node1 = [i for i,x in enumerate(taxon) if x.label==str(tips[len(tips)-1])]
  myMrca = pdm.mrca(taxon[node0[0]], taxon[node1[0]])

  return myMrca.index
# ----------------------------------------------------------------------------------------------------------------------
def my_merge_trees(recomtrees , recomnodes,clonal_tree):
  desc = []
  equ = np.zeros((len(recomtrees),4))
  for treeid in range(len(recomtrees)):
   rtree= Tree.get_from_string(recomtrees[treeid],schema='newick')
   set_index(rtree)
   equ[treeid,0] = recomnodes[treeid]
   equ[treeid,1:3] = give_equivalent_node(rtree)
   equ[treeid,3] = treeid

  s_equ = equ[equ[:,2].argsort()[::-1]]
  # print(s_equ)
  s_equ = give_older_recom(s_equ)
  # print(s_equ)

  clonaltree = Tree.get_from_string(clonal_tree,schema='newick')
  set_index(clonaltree)

  maintree = Tree.get_from_string(recomtrees[int(s_equ[0][3])],schema='newick')
  set_index(maintree)

  for i in range(1,len(s_equ)):
    temptree = Tree.get_from_string(recomtrees[int(s_equ[i][3])],schema='newick')
    set_index(temptree)
    node_maintree = int(s_equ[i][0])
    if node_maintree  >= tips_number:
      node_maintree = new_mrca(clonaltree,maintree,node_maintree,int(s_equ[i-1][0]))
    prunenode = maintree.find_node_with_label(str(node_maintree))
    node = temptree.find_node_with_label(str(int(s_equ[i][1])))
    parent = prunenode.parent_node
    if ((node.edge_length) < maintree.max_distance_from_root()) :
      # print(" *********** Stage Two ***********")
      ancestor = []
      for id,tmp_node in enumerate(prunenode.ancestor_iter()):
          ancestor.append(tmp_node)
          if node.edge_length < tmp_node.distance_from_tip() :
              attached_node = tmp_node
              attached_id = id
              break

    if (node_maintree  < tips_number):
      relocated_nodes = ancestor[attached_id-1]  # relocated node is the adjacent node of recombinant node
      maintree.prune_subtree(prunenode) # the original recombinant node was removed to reinsert in the other side
      attached_node.remove_child(relocated_nodes) # relocated node was removed to reinsert in to right side
      newborn = dendropy.datamodel.treemodel.Node()  # newborn in the new mrca of recombinant node and its sister
      newborn.add_child(prunenode)
      prunenode.edge_length = node.edge_length
      if attached_node.edge_length is None:
        attached_node.edge_length = 0
      newborn.edge_length = maintree.max_distance_from_root() - (node.edge_length + attached_node.edge_length)
      relocated_nodes.edge_length = relocated_nodes.edge_length - newborn.edge_length
      newborn.add_child(relocated_nodes)
      attached_node.add_child(newborn)
    else:
      relocated_nodes = ancestor[attached_id-1]  # relocated node is the adjacent node of recombinant node
      maintree.prune_subtree(prunenode) # the original recombinant node was removed to reinsert in the other side
      attached_node.remove_child(relocated_nodes) # relocated node was removed to reinsert in to right side
      newborn = dendropy.datamodel.treemodel.Node()  # newborn in the new mrca of recombinant node and its sister
      newborn.add_child(prunenode)
      prunenode.edge_length = node.edge_length
      newborn.edge_length = maintree.max_distance_from_root() - (node.edge_length + node.distance_from_tip())
      relocated_nodes.edge_length = relocated_nodes.edge_length - newborn.edge_length
      newborn.add_child(relocated_nodes)
      attached_node.add_child(newborn)

  return maintree.as_string(schema="newick")
# ----------------------------------------------------------------------------------------------------------------------
def merge_Recomtrees(recomtrees , recomnodes):

  clonaltree = Tree.get_from_string(clonal_tree,schema='newick')
  set_index(clonaltree)

  # if one one node is the subset of the other nodes
  conflict = common_nodes(clonaltree,recomnodes)

  # print(conflict)

  if (len(conflict) > 0)  and (conflict[0][0] == conflict[0][1]):
        # print("Same nodes")
        conflict = []

  if len(conflict) == 0 :
      # print("NOOOOOO conflict!")

      desc = []
      equ = np.zeros((len(recomtrees),4))
      for treeid in range(len(recomtrees)):
        rtree= Tree.get_from_string(recomtrees[treeid],schema='newick')
        set_index(rtree)
        equ[treeid,0] = recomnodes[treeid]
        equ[treeid,1:3] = give_equivalent_node(rtree)
        equ[treeid,3] = treeid

      s_equ = equ[equ[:,2].argsort()[::-1]]
      s_equ = give_older_recom(s_equ)
      # print(s_equ)

      maintree = Tree.get_from_string(recomtrees[int(s_equ[0][3])],schema='newick')
      set_index(maintree)

      for i in range(1,len(s_equ)):
        temptree = Tree.get_from_string(recomtrees[int(s_equ[i][3])],schema='newick')
        set_index(temptree)
        node_maintree = int(s_equ[i][0])
        if node_maintree  >= tips_number:
          node_maintree = new_mrca(clonaltree,maintree,node_maintree,int(s_equ[i-1][0]))
        prunenode = maintree.find_node_with_label(str(node_maintree))
        node = temptree.find_node_with_label(str(int(s_equ[i][1])))
        parent = prunenode.parent_node

        attached_node,attached_id,ancestor =  find_attached_node(maintree,node,prunenode)
        regraft_recomnodes(maintree,attached_node,attached_id,prunenode,node,ancestor)

      return maintree.as_string(schema="newick")

  elif len(conflict) > 0 :
    print("CONFLICT!! ")
    desc = []
    equ = np.zeros((len(recomtrees), 4))
    for treeid in range(len(recomtrees)):
      rtree = Tree.get_from_string(recomtrees[treeid], schema='newick')
      set_index(rtree)
      equ[treeid, 0] = recomnodes[treeid]
      equ[treeid, 1:3] = give_equivalent_node(rtree)
      equ[treeid, 3] = treeid

    print(equ)
    for i in range(len(conflict)):
      if conflict[i][0] < tips_number:
        # print("one node is taxa!")
        node1 = clonaltree.find_node_with_label(str(conflict[i][0]))
        node2 = clonaltree.find_node_with_label(str(conflict[i][1]))
        if node1.parent_node == node2:
          # print(node2, "is parent of", node1)
          [x,y] = np.where(equ==node2.index)
          # print(x)
          # print(y)

        equ = np.delete(equ, (x), axis=0)


      else:
        print("two nodes are internal!")

    print(equ)
    s_equ = equ[equ[:, 2].argsort()[::-1]]
    s_equ = give_older_recom(s_equ)
    # print(s_equ)

    maintree = Tree.get_from_string(recomtrees[int(s_equ[0][3])], schema='newick')
    set_index(maintree)

    for i in range(1, len(s_equ)):
      temptree = Tree.get_from_string(recomtrees[int(s_equ[i][3])], schema='newick')
      set_index(temptree)
      node_maintree = int(s_equ[i][0])
      if node_maintree >= tips_number:
        node_maintree = new_mrca(clonaltree, maintree, node_maintree, int(s_equ[i - 1][0]))
      prunenode = maintree.find_node_with_label(str(node_maintree))
      node = temptree.find_node_with_label(str(int(s_equ[i][1])))
      parent = prunenode.parent_node

      attached_node, attached_id, ancestor = find_attached_node(maintree, node, prunenode)
      regraft_recomnodes(maintree, attached_node, attached_id, prunenode, node, ancestor)

    return maintree.as_string(schema="newick")

# ----------------------------------------------------------------------------------------------------------------------
def resolve_conflict(clonaltree, recomtrees, recomnodes):
    conflict = common_nodes(clonaltree, recomnodes)
    print(recomnodes)
    print(conflict)
    if len(conflict) > 0:
      for i in range(len(conflict)):
        print(conflict[i][0], conflict[i][1])
        node1 = clonaltree.find_node_with_label(str(conflict[i][0]))
        node2 = clonaltree.find_node_with_label(str(conflict[i][1]))
        if node1.parent_node == node2:
          print(node2 , "is parent of" ,node1 )
          id_main = recomnodes.index(conflict[i][1])
          id_temp = recomnodes.index(conflict[i][0])
          maintree = Tree.get_from_string(recomtrees[id_main], schema='newick')
          set_index(maintree)
          temptree = Tree.get_from_string(recomtrees[id_temp], schema='newick')
          set_index(temptree)

          if int(node1.index) < tips_number:
            print("node is taxa")
            prunenode = node1
            recomnode = temptree.find_node_with_label(str(node1.index))
          else:
            print("node is internal")
            res = set()
            desc = give_descendents(clonaltree, node1.index, res)
            prune_index = my_mrca(maintree, list(desc))
            prunenode = maintree.find_node_with_label(str(prune_index))
            recom_index = my_mrca(temptree, list(desc))
            recomnode = temptree.find_node_with_label(str(recom_index))

          print(maintree.as_ascii_plot())
          print(maintree.as_string(schema="newick"))
          print("prunenode:")
          print(prunenode)
          print(prunenode.edge_length)
          print("recomnode:")
          print(recomnode)
          print(recomnode.edge_length)
          parent = prunenode.parent_node
          print("parent")
          print(parent)

          if ((recomnode.edge_length) > parent.distance_from_tip()) and (
                  (recomnode.edge_length) < tree.max_distance_from_root()):
            ancestor = []
            for id, tmp_node in enumerate(prunenode.ancestor_iter()):
              ancestor.append(tmp_node)
              # print(id ,"::::::",tmp_node.index)
              if recomnode.edge_length < tmp_node.distance_from_tip():
                attached_node = tmp_node
                attached_id = id
                # print(attached_node.index)
                break

            relocated_nodes = ancestor[attached_id - 1]  # relocated node is the adjacent node of recombinant node
            parent.remove_child(node)  # the original recombinant node was removed to reinsert in the other side
            attached_node.remove_child(relocated_nodes)  # relocated node was removed to reinsert in to right side
            newborn = dendropy.datamodel.treemodel.Node()  # newborn in the new mrca of recombinant node and its sister
            # newborn.edge_length = attached_node.distance_from_tip() - recom_length
            # node.edge_length = recom_length
            newborn.add_child(node)
            relocated_nodes.edge_length = relocated_nodes.edge_length - newborn.edge_length
            newborn.add_child(relocated_nodes)
            attached_node.add_child(newborn)
            print(maintree.as_string(schema="newick"))

            # if ((recomnode.edge_length) < maintree.max_distance_from_root()):
            #     print(" *********** Stage Two ***********")
            #     ancestor = []
            #     attached_id = 0
            #     for id, tmp_node in enumerate(prunenode.ancestor_iter()):
            #         # print(temp_node)
            #         ancestor.append(tmp_node)
            #         if recomnode.edge_length < tmp_node.distance_from_tip():
            #             attached_node = tmp_node
            #             attached_id = id
            #             break

            # if attached_id == 0 :
            #   parent.remove_child(prunenode)
            # else:

            # if (int(recomnode.index) < tips_number):
            # print("step child")
            #   relocated_nodes = ancestor[attached_id-1]  # relocated node is the adjacent node of recombinant node
            # maintree.prune_subtree(prunenode) # the original recombinant node was removed to reinsert in the other side
            #   attached_node.remove_child(relocated_nodes) # relocated node was removed to reinsert in to right side
            #   newborn = dendropy.datamodel.treemodel.Node()  # newborn in the new mrca of recombinant node and its sister
            #   newborn.add_child(prunenode)
            #   prunenode.edge_length = recomnode.edge_length
            #   if attached_node.edge_length is None:
            #     attached_node.edge_length = 0
            #   newborn.edge_length = maintree.max_distance_from_root() - (recomnode.edge_length + attached_node.edge_length)
            #   relocated_nodes.edge_length = relocated_nodes.edge_length - newborn.edge_length
            #   newborn.add_child(relocated_nodes)
            #   attached_node.add_child(newborn)
            # else:
            #   print("step internal")
            #   relocated_nodes = ancestor[attached_id-1]  # relocated node is the adjacent node of recombinant node
            #   if attached_node == maintree.seed_node:
            #     print("root")
            #     maintree.prune_subtree(relocated_nodes)
            #     new_tree = dendropy.Tree(taxon_namespace=taxa)
            #     other_nodes = maintree.seed_node
            #     new_tree.seed_node.add_child(other_nodes)
            #     newborn = dendropy.datamodel.treemodel.Node()
            #     newborn.add_child(recomnode)
            #     newborn.add_child(relocated_nodes)
            #     newborn.edge_length = (relocated_nodes.edge_length +relocated_nodes.distance_from_tip())  - (recomnode.edge_length + recomnode.distance_from_tip())
            #     relocated_nodes.edge_length = relocated_nodes.edge_length - newborn.edge_length
            #     new_tree.seed_node.add_child(newborn)
            #     new_tree.prune_subtree(prunenode)
            #   elif not attached_node == maintree.seed_node:
            #     print("not root")
            #     maintree.prune_subtree(prunenode) # the original recombinant node was removed to reinsert in the other side
            #     attached_node.remove_child(relocated_nodes) # relocated node was removed to reinsert in to right side
            #     newborn = dendropy.datamodel.treemodel.Node()  # newborn in the new mrca of recombinant node and its sister
            #     newborn.add_child(recomnode)
            #     newborn.edge_length = maintree.max_distance_from_root() - (recomnode.edge_length + recomnode.distance_from_tip())
            #     relocated_nodes.edge_length = relocated_nodes.edge_length - newborn.edge_length
            #     newborn.add_child(relocated_nodes)
            #     attached_node.add_child(newborn)

            print(maintree.as_ascii_plot())
            print(maintree.as_string(schema="newick"))

            # print(new_tree.as_ascii_plot())
      # print(new_tree.as_string(schema="newick"))


# ----------------------------------------------------------------------------------------------------------------------





tips_number = 10

clonal_tree = '((7:0.00564644732893243,(6:0.0018631361998130264,3:0.0018631361998130264):0.003783311129119403):0.004353552671067571,((5:0.0012983047690974484,0:0.0012983047690974484):0.0008742771902157827,((8:0.00032306378094596307,1:0.00032306378094596307):0.0009406763008812347,(4:0.0008421113229836173,(9:0.0002315542153263079,2:0.0002315542153263079):0.0006105571076573094):0.0004216287588435805):0.0009088418774860336):0.007827418040686769):0.0;'

recomnodes= [1, 8, 13]
recomtrees = ['(((7:0.00564644732893243,(6:0.0018631361998130264,3:0.0018631361998130264)10:0.003783311129119403)11:0.004353552671067571,((5:0.0012983047690974484,0:0.0012983047690974484)12:0.0008742771902157827,(8:0.0012637400818271978,(4:0.0008421113229836173,(9:0.0002315542153263079,2:0.0002315542153263079)14:0.0006105571076573094)15:0.0004216287588435805)16:0.0009088418774860336)17:0.007827418040686769)18:0.08678177784907884,1:0.09678177784907883);\n',
 '(((7:0.00564644732893243,(6:0.0018631361998130264,3:0.0018631361998130264)10:0.003783311129119403)11:0.004353552671067571,((5:0.0012983047690974484,0:0.0012983047690974484)12:0.0008742771902157827,(1:0.0012637400818271978,(4:0.0008421113229836173,(9:0.0002315542153263079,2:0.0002315542153263079)14:0.0006105571076573094)15:0.0004216287588435805)16:0.0009088418774860336)17:0.007827418040686769)18:0.0842418245159281,8:0.09424182451592811);\n',
 '(((7:0.00564644732893243,(6:0.0018631361998130264,3:0.0018631361998130264)10:0.003783311129119403)11:0.004353552671067571,((5:0.0012983047690974484,0:0.0012983047690974484)12:0.0008742771902157827,(4:0.0008421113229836173,(9:0.0002315542153263079,2:0.0002315542153263079)14:0.0006105571076573094)15:0.001330470636329614)17:0.007827418040686769)18:0.08746062965856694,(8:0.00032306378094596307,1:0.00032306378094596307)13:0.09713756587762098);\n']


# my_merge_trees(recomtrees,recomnodes,clonal_tree)

print(merge_Recomtrees(recomtrees , recomnodes))

# clonaltree = Tree.get_from_string(clonal_tree, schema='newick')
# set_index(clonaltree)
# resolve_conflict(clonaltree, recomtrees, recomnodes)



