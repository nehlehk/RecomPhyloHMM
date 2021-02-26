import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse


x = []
y1 = []
y2 = []
# y3 = []

parser = argparse.ArgumentParser(description='''You did not specify any parameters.''')
parser.add_argument('-f', "--rmseFiles", type=str, help='tree')
args = parser.parse_args()


# rmse_path = args.rmseFiles

rmse_path = '/home/nehleh/work/results/Summary_Results/collected_rmse.csv'


# files = glob.glob(rmse_path)
# for id,f in enumerate(files):
#     x.append(id)
#     csv = pd.read_csv(f, sep=';',names=None)
#     y1.append(float(csv.iloc[:,1].name))
#     y2.append(float(csv.iloc[:,2].name))
#     # y3.append(float(csv.iloc[:,3].name))


f = open(rmse_path, "r")
df = pd.read_csv(f, sep=';')
print(df)
print(df.loc[0])
# x.append((df.iloc[:,0]))
# print(x)
# y1.append(float(csv.iloc[:,1].name))
# # print(y1)
# y2.append(float(csv.iloc[:,2].name))
# print(y2)
# y3.append(float(csv.iloc[:,3].name))


# rmse_summary = pd.DataFrame({'id': x, 'phyloHMM': y1 ,'CFML':y2 })
# rmse_summary.to_csv('./rmse_summary.csv', sep='\t', header=True)




# fig = plt.figure(figsize=(6,6))
# ax1 = fig.add_subplot(2, 1, 1)
# ax1.plot(x,y1, marker='.' ,label= "phyloHMM_realRecom_RMSE")
# ax1.plot(x,y2, marker='.' ,label= "CFML_realRecom_RMSE")
# # ax1.plot(x,y3, marker='.' ,label= "rmse_clonal_real")
# ax1.set_title('RMSE values for different dataset')
# ax1.legend(loc='best')
#
# ax2 = fig.add_subplot(2, 1, 2)
# index = ['RMSE comparison']
# X2 = np.arange(len(index))
# ax2.bar(0.00, [np.average(y1)], color = '#0F52BA', width = 0.25 , label = 'rmse_real_phyloHMM')
# ax2.bar(0.25, [np.average(y2)], color = '#6593F5', width = 0.25 , label = 'rmse_real_CFML')
# # ax2.bar(0.50, [np.average(y3)], color = '#73C2FB', width = 0.25 , label = 'rmse_clonal_real')
# ax2.set_xticks(X2 + 0.25)
# ax2.set_xticklabels(index)
# ax2.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.3)
# ax2.legend()
#
#
# plt.savefig("RMSE_comparison.jpeg")
# plt.show()

