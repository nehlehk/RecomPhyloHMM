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

rmse_path = '/home/nehleh/work/results/Summary_Results/rmse.csv'



f = open(rmse_path, "r")
df = pd.read_csv(f, sep=';', names=['nu', 'prob','Phylohmm','CFML'], header=None)
print(df)
nu_value = df['nu'].unique()
print(nu_value)
prob_value = df['prob'].unique()
print(prob_value)
case_num = df['CFML'].nunique()
print(case_num)

print(df.groupby(['nu', 'prob'], as_index=False).mean())


fig = plt.figure(figsize=(10,10))
# ax = plt.gca()
ax1 = fig.add_subplot(2, 1, 1)
for nu in nu_value:
    for prob in prob_value:
        # print("nu:",nu,"prob:",prob)
        my_df = df.loc[(df['nu'] == nu) & (df['prob'] == prob)]
        my_df = my_df.reset_index(drop=True)
        my_df.plot(kind='line', y='Phylohmm', ax=ax1, marker='*', label='RMSE values when nu =' + str(nu) + ' and prob =' + str(prob))
my_df.plot(kind='line', y='CFML',color='blue', ax=ax1, marker='*')

plt.xlabel('case_number')
plt.ylabel('RMSE-value')
ax1.legend(loc='best')

ax2 = fig.add_subplot(2, 1, 2)
# df.groupby(['nu', 'prob'], as_index=False).mean().unstack().plot(kind='bar',stacked=True)
df.groupby(['nu', 'prob'], as_index=False).mean().plot(kind='bar',stacked=True , y = ['Phylohmm','CFML'] , ax = ax2)
# df.groupby(['nu', 'prob'], as_index=False).mean().plot(kind='bar',stacked=True , y = 'CFML', ax = ax2)
ax2.legend()



plt.show()





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

