import pandas as pd
import matplotlib.pyplot as plt
import argparse
import seaborn as sns
# from pycallgraph import PyCallGraph
# from pycallgraph.output import GraphvizOutput

parser = argparse.ArgumentParser(description='''You did not specify any parameters.''')
parser.add_argument('-f', "--rmseFiles", type=str, help='rmseFiles')
args = parser.parse_args()


rmse_path = args.rmseFiles

# with PyCallGraph(output=GraphvizOutput()):
# rmse_path = '/home/nehleh/Desktop/rmse30.csv'



f = open(rmse_path, "r")
df = pd.read_csv(f, sep=';', names=['nu', 'prob','Phylohmm','CFML'], header=None)
print(df)
# nu_value = df['nu'].unique()
# prob_value = df['prob'].unique()
# case_num = df['CFML'].nunique()


# print(pd.melt(df))
# new_df = pd.melt(df, id_vars=['nu','prob'], value_vars=['Phylohmm','CFML'])
# # print(new_df)
# # new_df.boxplot( column='value' , by=['nu'] ,grid=True, rot=45, fontsize=10)
# ax = sns.boxplot(data=new_df, x='variable', y='value', hue='nu')

# sns.boxplot(x="variable", y="value", data=pd.melt(df))



fig = plt.figure(figsize=(8,8))
ax1 = fig.add_subplot(2, 1, 1)
ax1 = sns.boxplot(x = 'nu', y="Phylohmm" ,data=df , hue = 'prob' )
ax1 = sns.stripplot(x = 'nu', y="Phylohmm",data=df, hue='prob' , jitter=True, dodge=True, marker='o', color=".5")
handles, labels = ax1.get_legend_handles_labels()
l = plt.legend(handles[0:5], labels[0:5] , title= 'Probablity')


ax2 = fig.add_subplot(2, 1, 2)
ax2 = sns.boxplot(x = 'nu',y=df["CFML"] ,data=df  )
# ax2 = sns.stripplot(x = 'nu', y="CFML",data=df , jitter=True, dodge=True, marker='o', color=".5")
# handles, labels = ax2.get_legend_handles_labels()
# l = plt.legend(handles[0:5], labels[0:5] , title= 'Probablity')




# df.boxplot( column=['Phylohmm' ,'CFML'] , by=['nu', 'prob'] ,  layout=(1, 2) ,grid=True, rot=45, fontsize=10)


# fig = plt.figure(figsize=(10,10))
# ax1 = fig.add_subplot(2, 1, 1)
# for nu in nu_value:
#     for prob in prob_value:
#         # print("nu:",nu,"prob:",prob)
#         my_df = df.loc[(df['nu'] == nu) & (df['prob'] == prob)]
#         my_df = my_df.reset_index(drop=True)
#         my_df.plot(kind='line', y='Phylohmm', ax=ax1, marker='*', label='PhyloHMM-RMSE when nu =' + str(nu) + ' and prob =' + str(prob))
# my_df.plot(kind='line', y='CFML',color='blue', ax=ax1, marker='*' , label='CFML-RMSE')
#
# plt.xlabel('case_number')
# plt.ylabel('RMSE-value')
# # plt.xticks([])
# ax1.legend(bbox_to_anchor=(0.47, 0.9), prop={'size': 10})
#
# ax2 = fig.add_subplot(2, 1, 2)
# ax2.set_title('Mean RMSE values for different nu and probabilty')
# df.groupby(['nu', 'prob']).mean().plot(kind='bar',stacked=False, y = ['Phylohmm','CFML'] , ax = ax2 ,color = ['#0F52BA','#6593F5'] )
# plt.ylabel('RMSE-value')
# plt.xlabel('(nu,probabilty)')
# plt.xticks(rotation=0)
# ax2.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.3)
# ax2.legend()


# fig = plt.figure(figsize=(10,10))
# ax1 = fig.add_subplot(2, 1, 1)
# for nu in nu_value:
#     for prob in prob_value:
#         # print("nu:",nu,"prob:",prob)
#         my_df = df.loc[(df['nu'] == nu) & (df['prob'] == prob)]
#         my_df.columns = ['nu', 'prob','Phylohmm'+str(nu)+'_'+str(prob),'CFML']
#         my_df = my_df.reset_index(drop=True)
#         my_df.boxplot(column=['Phylohmm'+str(nu)+'_'+str(prob),'CFML'])
# my_df.boxplot( column=['CFML'])

# plt.xlabel('case_number')
# plt.ylabel('RMSE-value')
# # plt.xticks([])
# ax1.legend(bbox_to_anchor=(0.47, 0.9), prop={'size': 10})
#
# ax2 = fig.add_subplot(2, 1, 2)
# ax2.set_title('Mean RMSE values for different nu and probabilty')
# df.groupby(['nu', 'prob']).mean().plot(kind='bar',stacked=False, y = ['Phylohmm','CFML'] , ax = ax2 ,color = ['#0F52BA','#6593F5'] )
# plt.ylabel('RMSE-value')
# plt.xlabel('(nu,probabilty)')
# plt.xticks(rotation=0)
# ax2.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.3)
# ax2.legend()

# plt.show()
plt.savefig("RMSE_comparison.jpeg")