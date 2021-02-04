from sklearn.metrics import mean_squared_error
import numpy as np
import argparse
import math
import csv




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
def CFML_recombination(CFML_recomLog):
    CFMLData = np.zeros((alignment_len, nodes_number))
    df = pd.read_csv(CFML_recomLog,sep='\t', engine='python')
    # print(df)
    for i in range(len(df)):
        s = df['Beg'][i]
        e = df['End'][i]
        node = df['Node'][i]
        if "NODE_" in str(node):
          node = node[5:]
        CFMLData[s:e,int(give_taxon_index(tree, int(node)))] = 1

    return CFMLData
# **********************************************************************************************************************
def calc_rmse(data1,data2):
    mse = []
    rmse = []
    for t in range(tips_num):
        m = mean_squared_error(data1[:,t], data2[:,t])
        r = math.sqrt(m)
        mse.append(m)
        rmse.append(r)

    return sum(rmse)
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
def write_rmse(nu,rmse_real_predict,rmse_real_CFML):
    with open('rmse.rmse', mode='w') as rmse_file:
        rmse_writer = csv.writer(rmse_file, delimiter=';', quotechar='"', quoting=csv.QUOTE_MINIMAL)

        # rmse_writer.writerow(['nu_hmm','rmse_real_predict','rmse_clonal_predict','rmse_clonal_real'])
        rmse_writer.writerow([nu, rmse_real_predict,rmse_real_CFML])
# **********************************************************************************************************************


if __name__ == "__main__":


    parser = argparse.ArgumentParser(description='''You did not specify any parameters.''')


    parser.add_argument('-c', "--cfmlFile", type=str, required=True, help='cfmlFile')
    args = parser.parse_args()

    recomLog = args.recomlogFile
    cfml_path = args.cfmlFile


    CFMLData = CFML_recombination(cfml_path)

    rmse_real_phyloHMM= calc_rmse(realData,phyloHMMData)
    # print(rmse_real_phyloHMM)
    rmse_real_CFML = calc_rmse(realData,CFMLData)
    # print(rmse_real_CFML)

    write_rmse(nu, rmse_real_phyloHMM, rmse_real_CFML)