import numpy as np
import dendropy
from dendropy.simulate import treesim
import random
from dendropy import Tree
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
import argparse




# ----------------------------------------------------------------------------------------------------------------------
def ex_recom_maker(tree ,node ,nu ,taxa):
    rand_nu = np.random.normal(nu,0.01)
    co_recom = rand_nu/2
    new_tree = dendropy.Tree(taxon_namespace=taxa)
    # changeing in topology when recombiantion is larger than tree.max_distance_from_root()
    if co_recom > tree.max_distance_from_root():
      if (node.edge_length is None):
        node.edge.length = 0
      # new_tree = dendropy.Tree(taxon_namespace=taxa)
      if node.is_leaf():
        external_len = co_recom + node.edge_length - tree.max_distance_from_root()
        tree.prune_subtree(node)
        other_nodes = tree.seed_node
        new_tree.seed_node.add_child(other_nodes)
        new_tree.seed_node.add_child(node)
        other_nodes.edge_length = external_len
        node.edge_length = co_recom + node.edge_length
      else:
        tree.prune_subtree(node)
        other_nodes = tree.seed_node
        new_tree.seed_node.add_child(other_nodes)
        new_tree.seed_node.add_child(node)
        node.edge_length = co_recom + node.edge_length
        other_nodes.edge_length = node.distance_from_tip() + node.edge_length - other_nodes.distance_from_tip()
    #*************************************************************************
    # topology does not change in this case:
    elif ((co_recom + node.edge_length) < parent.distance_from_tip()) and (node.is_leaf):
        print(" *********** Stage one --- leaf  ***********")
        node.edge.length = node.edge.length + co_recom
        sister = node.sister_nodes()
        sister[0].edge.length = sister[0].edge.length + co_recom
        parent.edge.length = parent.edge.length - co_recom
        new_tree = tree
    elif ((co_recom + node.edge_length) < parent.distance_from_tip()) and (node.is_internal):
        print(" *********** Stage one --- internal  ***********")
        node.edge.length = node.edge.length + co_recom
        new_tree = tree
    #*************************************************************************
    # changing in topology to make recombination tree:
    elif ((co_recom + node.edge_length) > parent.distance_from_tip())  and ((co_recom + node.edge_length) < tree.max_distance_from_root()):
        ancestor = []
        recom_length = co_recom + node.edge_length
        for id,tmp_node in enumerate(node.ancestor_iter()):
            ancestor.append(tmp_node)
            # print(id ,"::::::",tmp_node.index)
            if recom_length < tmp_node.distance_from_tip() :
                attached_node = tmp_node
                attached_id = id
                # print(attached_node.index)
                break

        relocated_nodes = ancestor[attached_id-1]  # relocated node is the adjacent node of recombinant node
        parent.remove_child(node)         # the original recombinant node was removed to reinsert in the other side
        attached_node.remove_child(relocated_nodes) # relocated node was removed to reinsert in to right side
        newborn = dendropy.datamodel.treemodel.Node()  # newborn in the new mrca of recombinant node and its sister
        newborn.edge_length = attached_node.distance_from_tip() - recom_length
        node.edge_length = recom_length
        newborn.add_child(node)
        relocated_nodes.edge_length = relocated_nodes.edge_length - newborn.edge_length
        newborn.add_child(relocated_nodes)
        attached_node.add_child(newborn)
        recombination_trees.append(tree.as_string(schema="newick"))
    #*************************************************************************
    return new_tree.as_string(schema="newick")
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
def give_taxon(tree,index):
    for node in tree.postorder_node_iter():
        if int(node.index) == index:
            return int(node.taxon.label)
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
        result.add(give_taxon(tree,r_node))
  return result
# ----------------------------------------------------------------------------------------------------------------------
def give_equivalent_node_byedge(recomtree):
  e = []
  for edge in recomtree.postorder_edge_iter():
    if edge.length is None:
      edge.length = 0
    e.append(edge.length)

  m = np.max(e)
  for node in recomtree.postorder_node_iter():
    if node.edge_length == m:
      return node.label,m
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
def my_mrca(tree,tips):
  pdm = tree.phylogenetic_distance_matrix()
  taxon = tree.taxon_namespace
  node0 = [i for i,x in enumerate(taxon) if x.label==str(tips[0])]
  node1 = [i for i,x in enumerate(taxon) if x.label==str(tips[len(tips)-1])]
  myMrca = pdm.mrca(taxon[node0[0]], taxon[node1[0]])

  return myMrca.index
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
def give_older_recom(s_equ):
  count = Counter(a[0] for a in s_equ)
  max = np.max([a[2] for a in s_equ])
  [a for a in s_equ if (count[a[0]] == 1) or (a[2]==max)]
  return [a for a in s_equ if (count[a[0]] == 1) or (a[2]==max)]
# ----------------------------------------------------------------------------------------------------------------------
def remove_internal_labels(strtree):
  tree = Tree.get_from_string(strtree,schema='newick')
  for node in tree.postorder_node_iter():
    if not node.label is None:
      if (int(node.label)>= tips_number):
        node.label = None

  return tree.as_string(schema="newick")

# ----------------------------------------------------------------------------------------------------------------------
def my_merge_trees(recomtrees , recomnodes):
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
def merge_Recomtrees(recomtrees , recomnodes):

  clonaltree = Tree.get_from_string(clonal_tree,schema='newick')
  set_index(clonaltree)

  # if one one node is the subset of the other nodes
  conflict = common_nodes(clonaltree,recomnodes)

  if (len(conflict) > 0)  and (conflict[0][0] == conflict[0][1]):
        print("Same nodes")
        conflict = []

  if len(conflict) == 0 :
      print("NOOOOOO conflict!")

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
    for i in range(len(conflict)):
      if conflict[i][0] < tips_number:
        print("one node is taxa!")
      else:
        print("two nodes are internal!")

# ----------------------------------------------------------------------------------------------------------------------
def resolve_conflict(clonaltree, recomtrees, recomnodes):
    conflict = common_nodes(clonaltree, recomnodes)
    print(recomnodes)
    if len(conflict) > 0:
        for i in range(len(conflict)):
            print(conflict[i][0], conflict[i][1])
            node1 = clonaltree.find_node_with_label(str(conflict[i][0]))
            node2 = clonaltree.find_node_with_label(str(conflict[i][1]))
            if node1.parent_node == node2:
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

                if ((recomnode.edge_length) > parent.distance_from_tip()) and ((recomnode.edge_length) < tree.max_distance_from_root()):
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
def edit_merge_trees(recomtrees , recomnodes):
  desc = []
  equ = np.zeros((len(recomtrees),4))

  for treeid in range(len(recomtrees)):
   rtree= Tree.get_from_string(recomtrees[treeid],schema='newick')
   set_index(rtree)
   equ[treeid,0] = recomnodes[treeid]
   equ[treeid,1:3] = give_equivalent_node(rtree)
   equ[treeid,3] = treeid

  s_equ = equ[equ[:,2].argsort()[::-1]]
  print(s_equ)
  s_equ = give_older_recom(s_equ)
  print(s_equ)

  clonaltree = Tree.get_from_string(clonal_tree,schema='newick')
  set_index(clonaltree)

  maintree = Tree.get_from_string(recomtrees[int(s_equ[0][3])],schema='newick')
  set_index(maintree)
  print(maintree.as_ascii_plot())
  print(maintree.as_string(schema="newick"))

  for i in range(1,len(s_equ)):
    temptree = Tree.get_from_string(recomtrees[int(s_equ[i][3])],schema='newick')
    set_index(temptree)
    node_maintree = int(s_equ[i][0])
    print(temptree.as_ascii_plot())
    print(temptree.as_string(schema="newick"))
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
def full_merge_trees(recomtrees , recomnodes):
  desc = []
  equ = np.zeros((len(recomtrees),4))
  for treeid in range(len(recomtrees)):
   rtree= Tree.get_from_string(recomtrees[treeid],schema='newick')
   set_index(rtree)
   equ[treeid,0] = recomnodes[treeid]
   equ[treeid,1:3] = give_equivalent_node(rtree)
   equ[treeid,3] = treeid

  s_equ = equ[equ[:,2].argsort()[::-1]]
  print(s_equ)

  s_equ = give_older_recom(s_equ)
  # print(s_equ)

  clonaltree = Tree.get_from_string(clonal_tree,schema='newick')
  set_index(clonaltree)

  maintree = Tree.get_from_string(recomtrees[int(s_equ[0][3])],schema='newick')
  set_index(maintree)
  # print(maintree.as_ascii_plot())
  # print(maintree.as_string(schema="newick"))

  for i in range(1,len(s_equ)):
    # print(maintree.as_ascii_plot())
    # print(maintree.as_string(schema="newick"))
    temptree = Tree.get_from_string(recomtrees[int(s_equ[i][3])],schema='newick')
    set_index(temptree)
    print(temptree.as_ascii_plot())
    print(temptree.as_string(schema="newick"))
    node_maintree = int(s_equ[i][0])
    print(node_maintree)
    print(s_equ[i-1][0])
    if node_maintree  >= tips_number:
      node_maintree = new_mrca(clonaltree,maintree,node_maintree,int(s_equ[i-1][0]))
      print("node_maintree")
      print(node_maintree)
    prunenode = maintree.find_node_with_label(str(node_maintree))
    print("prunenode:")
    print(prunenode)
    print(prunenode.edge_length)
    node = temptree.find_node_with_label(str(int(s_equ[i][1])))
    print("node:")
    print(node)
    print(node.edge_length)

    parent = prunenode.parent_node

    print("parent:")
    print(parent)


    if ((node.edge_length) < maintree.max_distance_from_root()) :
      # print(" *********** Stage Two ***********")
      ancestor = []
      for id,tmp_node in enumerate(prunenode.ancestor_iter()):
          ancestor.append(tmp_node)
          # print(id ,"::::::",tmp_node.index)
          if node.edge_length < tmp_node.distance_from_tip() :
              attached_node = tmp_node
              attached_id = id
              print("attached_node")
              print(attached_node)
              break

    if (node_maintree  < tips_number):
      relocated_nodes = ancestor[attached_id-1]  # relocated node is the adjacent node of recombinant node
      print("relocated_nodes:")
      print(relocated_nodes)
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
      print("relocated_nodes:")
      print(relocated_nodes)
      maintree.prune_subtree(prunenode) # the original recombinant node was removed to reinsert in the other side
      attached_node.remove_child(relocated_nodes) # relocated node was removed to reinsert in to right side
      newborn = dendropy.datamodel.treemodel.Node()  # newborn in the new mrca of recombinant node and its sister
      newborn.add_child(prunenode)
      prunenode.edge_length = node.edge_length
      newborn.edge_length = maintree.max_distance_from_root() - (node.edge_length + node.distance_from_tip())
      relocated_nodes.edge_length = relocated_nodes.edge_length - newborn.edge_length
      newborn.add_child(relocated_nodes)
      attached_node.add_child(newborn)


  print(maintree.as_ascii_plot())
  print(maintree.as_string(schema="newick"))
  return maintree.as_string(schema="newick")

# ----------------------------------------------------------------------------------------------------------------------
def resolve_parent_nodes(s_equ, tree, temptree):
    desc = []
    n = []
    length = []
    tip_len = []
    for i in range(len(s_equ)):
        my_node = int(s_equ[i][0])
        # print(my_node)
        if my_node >= tips_number:
            res = set()
            desc.append(give_descendents(tree, my_node, res))
            n.append(my_node)
            length.append(s_equ[i][2])
            tip_len.append(temptree.find_node_with_label(str(int(s_equ[i][1]))).distance_from_tip())

    for i in range(len(desc)):
        for j in range(len(s_equ)):
            if (int(s_equ[j][0]) < tips_number) and (int(s_equ[j][0]) in desc[i]):
                tmp = list(desc[i])
                # print(int(s_equ[j][0]) , list(desc[i]) , n[i])
                for k in range(len(tmp)):
                    print(tmp[k], tmp[k], length[i] + tip_len[i], n[i])

    return desc
# ----------------------------------------------------------------------------------------------------------------------
def make_clonaltree(tips_number,max_tMRCA):
    taxon_list= []
    for i in range(tips_number):
      taxon_list.append(str(i))
    taxa = dendropy.TaxonNamespace(taxon_list)
    tree = treesim.pure_kingman_tree(taxon_namespace=taxa,pop_size=3)
    t = tree.max_distance_from_root()
    normal_co = t/ max_tMRCA
    tree.scale_edges(1/normal_co)
    clonal_tree = tree.as_string(schema="newick")
    clonal_tree = clonal_tree.replace('\n',"")
    myfile = open('./clonaltree.tree', 'w')
    myfile.write(clonal_tree)
    myfile.close()

    return tree, clonal_tree, taxa
# ----------------------------------------------------------------------------------------------------------------------
def give_recom_num(tree,recom_rate,alignment_len):
    # Poisson( tree.sum() * rel_recomb_rate_per_site * alignment_length)
    recom_num =  np.random.poisson(tree.length() * recom_rate * alignment_len)
    if recom_num == 0:
        recom_num = 1
    return recom_num
# ----------------------------------------------------------------------------------------------------------------------
def make_nodes_weight(tree,status):
    set_index(tree)
    node_labels = []
    node_weight = []
    for node in tree.postorder_node_iter():
      if not node==tree.seed_node:
        if status == 1:
          node_labels.append(node.index)
          node_weight.append(node.edge_length /tree.length())
        elif status == 0:
          if node.is_leaf():
            node_labels.append(node.index)
            node_weight.append(node.edge_length /tree.length())
        elif status == 2:
          if node.is_internal():
            node_labels.append(node.index)
            node_weight.append(node.edge_length /tree.length())

    return node_labels,node_weight
# ----------------------------------------------------------------------------------------------------------------------
def recom_on_alignment(recom_num,recom_len,alignment_len,clonal_tree,node_labels,node_weight,nu_ex,taxa):
    ideal_gap = int(alignment_len/recom_num)
    my_trees = []
    nodes = []
    starts = []
    ends = []
    recomlens = []
    for recom_id in range(int(recom_num)):
      starting_falg = False
      tree = Tree.get_from_string(clonal_tree,schema='newick')
      set_index(tree)
      while not starting_falg:
        random_tip = np.random.choice(node_labels,1,node_weight)
        start_pos = random.randint(ideal_gap*recom_id, ideal_gap * (recom_id+1) )
        r_len = random.randint(threshold_len, recom_len + threshold_len)
        if (start_pos + r_len <= alignment_len) :
          nodes.append(random_tip)
          starts.append(start_pos)
          recomlens.append(r_len)
          ends.append(start_pos + r_len)
          recom_node = tree.find_node_with_label(str(int(random_tip)))
          recom_tree= ex_recom_maker(tree,recom_node,nu_ex,taxa) # make external recombination
          my_trees.append(recom_tree)
          starting_falg = True


    all_data = {'nodes':nodes , 'start':starts , 'end':ends, 'len':recomlens , 'tree':my_trees}
    df = pd.DataFrame(all_data)
    df.to_csv('./Recombination_Log.txt', sep='\t', header=True)

    return df,all_data
# ----------------------------------------------------------------------------------------------------------------------
def make_recom_fig(all_data,alignment_len,nodes_number,tips_number,clonal_tree):
    output = np.zeros((alignment_len, nodes_number))
    for i in range(nodes_number):
      for j in range(recom_num):
        if int(all_data['nodes'][j]) == i:
          s = int(all_data['start'][j])
          e = int(all_data['end'][j])
          output[s:e,i] = 1

    fig = plt.figure(figsize=(tips_number+9,tips_number/2))
    color=['red', 'green' ,'purple', 'blue','black']
    clonaltree = Tree.get_from_string(clonal_tree,schema='newick')
    set_index(clonaltree)
    # for i in range(nodes_number):
    for i in range(tips_number):
      # ax = fig.add_subplot(nodes_number,1,i+1)
      ax = fig.add_subplot(tips_number, 1, i + 1)
      # if i >= tips_number:
      #   desc = set()
      #   d = give_descendents(clonaltree,i,desc)
      #   ax.plot(output[:,i] ,label= str(i)+' is mrca:'+ str(d) ,color = color[i%5])
      # else:
      # if i < tips_number:
      ax.plot(output[:,i] ,label= i ,color = color[i%5])
      ax.legend( bbox_to_anchor=(0.045, 1.5) ,prop={'size':10} )
      ax.set_frame_on(False)
      ax.axis('off')

    ax.axis('on')
    ax.set_yticklabels([])
    plt.savefig("./BaciSim_Recombination.jpeg")
# ----------------------------------------------------------------------------------------------------------------------
def generate_final_report(df,alignment_len,clonal_tree,tips_number):
    # endpoints = df.stack().sort_values().reset_index(drop=True)
    endpoints = df[['start', 'end']].stack().sort_values().reset_index(drop=True)
    intervals = pd.DataFrame({'start': endpoints.shift().fillna(0), 'end': endpoints}).astype(int)
    # construct the list of intervals from the endpoints
    intervals['intv'] = [pd.Interval(a, b) for a, b in zip(intervals.start, intervals.end)]
    # these are the original intervals
    orig_invt = pd.arrays.IntervalArray([pd.Interval(a, b) for a, b in zip(df.start, df.end)])
    # walk through the intervals and compute the intersections
    intervals['total'] = intervals.intv.apply(lambda x: orig_invt.overlaps(x).sum())
    bounds = np.unique(df[['start', 'end']])
    if 0 not in bounds: bounds = np.insert(bounds, 0, 0)
    end = alignment_len
    bounds = np.append(bounds, end)
    total = []
    interval = []
    stat = []
    final_len = []
    nodes = []
    r_trees = []
    final_tree = []
    clonaltree = Tree.get_from_string(clonal_tree, schema='newick')
    set_index(clonaltree)
    children = []
    for i in range(len(bounds) - 1):
        # Find which intervals fit
        ix = (df['start'] <= bounds[i]) & (df['end'] >= bounds[i + 1])
        final_len.append(bounds[i + 1] - bounds[i])
        total.append(np.sum(ix))

        interval.append(df[ix].values.tolist())

        temp_node = []
        temp_tree = []
        kids = []
        temp = df[ix].values.tolist()
        for j in range(len(temp)):
            temp_node.append(int(temp[j][0]))
            if int(temp[j][0]) >= tips_number:
                desc = set()
                d = give_descendents(clonaltree, int(temp[j][0]), desc)
                kids.append(d)
            else:
                kids.append("")
            temp_tree.append(temp[j][4])
        nodes.append(temp_node)
        r_trees.append(temp_tree)
        children.append(kids)

        if (np.sum(ix) == 0):
            stat.append('clonal')
            final_tree.append(clonal_tree)
        elif (np.sum(ix) == 1):
            stat.append("recom")
            final_tree.append(r_trees[i])
        else:
            stat.append("overlap")
            # final_tree.append("")
            # final_tree.append(my_merge_trees(r_trees[i],nodes[i]))
            final_tree.append(merge_conflict(r_trees[i], nodes[i]))


    final = pd.DataFrame({'start': bounds[:-1], 'end': bounds[1:] ,'nodes':nodes, 'descendants':children,  'len':final_len , 'status':stat ,'final_tree': final_tree ,'total': total ,'tree':r_trees })
    final[['nodes', 'start', 'end', 'len', 'descendants', 'status']].to_csv('./BaciSim_Log.txt', sep='\t', header=True)

    myfile = open('./BaciSimTrees.tree', 'w')
    total = 0
    for id in range(len(final)):
        if len(final['final_tree'][id]) == 1:
            end_tree = remove_internal_labels(final['final_tree'][id][0])
        else:
            end_tree = remove_internal_labels(final['final_tree'][id])
        tmp = "[" + str(final['len'][id]) + "]" + " " + end_tree
        # print(tmp)
        myfile.write("%s" % tmp)

    myfile.close()

    return final
# ----------------------------------------------------------------------------------------------------------------------

def merge_conflict(recomtrees , recomnodes):

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
        print("one node is taxa!")
        node1 = clonaltree.find_node_with_label(str(conflict[i][0]))
        node2 = clonaltree.find_node_with_label(str(conflict[i][1]))
        if node1.parent_node == node2:
          print(node2, "is parent of", node1)
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



if __name__ == "__main__":

    parser=argparse.ArgumentParser(
        description='''You did not specify any parameters. You must at least specify the number of chromosomes sampled and the sequence length. ''',
        epilog="""All's well that ends well.""")
    parser.add_argument('-n', "--tips_number", type=int, default=10 , help='Sets the number of isolates (default is 10)')
    parser.add_argument('-g', "--alignment_len", type=int, default=1000 , help='Sets the number and lengths of fragments of genetic material (default is 5000)')
    parser.add_argument('-l', "--recom_len", type=int, default=200, help='Sets the average length of an external recombinant interval, (default is 500)')
    parser.add_argument('-r', "--recom_rate",type=float, default=0.05, help='Sets the site-specific rate of external (between species) recombination, (default is 0.05)')
    parser.add_argument('-nu',"--nu" ,  type=float, default=0.2, help='nu')
    parser.add_argument('-s',"--status" ,  type=int, default=0, help='0 is just leaves, 1 is for both internal nodes and leaves and 2 is just internal nodes')
    parser.add_argument('-t', "--tMRCA", type=float, default=0.01 ,help='tMRCA')


    # Read arguments from command line
    args = parser.parse_args()

    tips_number = args.tips_number
    alignment_len = args.alignment_len
    recom_len = args.recom_len
    recom_rate = args.recom_rate
    nu_ex = args.nu
    status = args.status
    tMRCA = args.tMRCA
    threshold_len = 200



    tree, clonal_tree, taxa = make_clonaltree(tips_number, tMRCA)
    nodes_number = len(tree.nodes())
    recom_num = give_recom_num(tree,recom_rate,alignment_len)
    node_labels,node_weight = make_nodes_weight(tree, status)
    df,all_data = recom_on_alignment(recom_num, recom_len, alignment_len, clonal_tree, node_labels, node_weight, nu_ex, taxa)
    make_recom_fig(all_data,alignment_len, nodes_number, tips_number, clonal_tree)
    final_report = generate_final_report(df, alignment_len, clonal_tree, tips_number)


    # # print(final_report)
    # plt.show()









