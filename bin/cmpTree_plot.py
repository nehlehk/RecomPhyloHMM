import pandas as pd
import matplotlib.pyplot as plt
import argparse
import seaborn as sns
import numpy as np




# parser = argparse.ArgumentParser(description='''You did not specify any parameters.''')
# parser.add_argument('-c', "--cmpTreesResult", type=str, help='cmpTreesResult')
# args = parser.parse_args()
#
# cmpTrees = args.cmpTreesResult

cmpTrees = '/home/nehleh/Desktop/all_cmpTrees.result'

f = open(cmpTrees, "r")
df = pd.read_csv(f,sep='\t', names=['No', 'clonalTree','otherTrees', 'Quartet' , 'PathDiffernce' , 'RF' , 'MatchingSplit' , 'UMAST' , 'RFWeighted' , 'GeoUnrooted' ] ,header=None)
new_df= df.loc[df['otherTrees'].between('1' , '5' , inclusive=True )].apply(pd.to_numeric)
new_df['otherTrees'].replace(to_replace=[1,2,3,4,5],value=['PhyloHMM','CFML','Gubbins','RaxML','Beast'] ,inplace=True)
# print(new_df)
g_df = new_df.loc[(new_df['otherTrees'] != 'Gubbins') ]

# print(g_df)

fig = plt.figure(figsize=(16,16))


ax1 = fig.add_subplot(2, 4, 1)
ax1 = sns.boxplot(x = 'otherTrees', y="Quartet" ,data=new_df)
ax1 = sns.stripplot(x = 'otherTrees', y="Quartet" ,data=new_df,  jitter=True, dodge=True, marker='o', color=".1")
ax1.set_title('Quartet' , fontsize=9)


ax2 = fig.add_subplot(2, 4, 2)
ax2 = sns.boxplot(x = 'otherTrees', y="PathDiffernce" ,data=new_df )
ax2 = sns.stripplot(x = 'otherTrees', y="PathDiffernce" ,data=new_df,  jitter=True, dodge=True, marker='o', color=".1")
ax2.set_title('PathDiffernce' , fontsize=9)

ax3 = fig.add_subplot(2, 4, 3)
ax3 = sns.boxplot(x = 'otherTrees', y="RF" ,data=new_df )
ax3 = sns.stripplot(x = 'otherTrees', y="RF" ,data=new_df,  jitter=True, dodge=True, marker='o', color=".1")
ax3.set_title('RF' , fontsize=9)

ax4 = fig.add_subplot(2, 4, 4)
ax4 = sns.boxplot(x = 'otherTrees', y="MatchingSplit" ,data=new_df )
ax4 = sns.stripplot(x = 'otherTrees', y="MatchingSplit" ,data=new_df,  jitter=True, dodge=True, marker='o', color=".1")
ax4.set_title('MatchingSplit' , fontsize=9)

ax5 = fig.add_subplot(2, 4, 5)
ax5 = sns.boxplot(x = 'otherTrees', y="UMAST" ,data=new_df )
ax5 = sns.stripplot(x = 'otherTrees', y="UMAST" ,data=new_df,  jitter=True, dodge=True, marker='o', color=".1")
ax5.set_title('UMAST' , fontsize=9)

ax6 = fig.add_subplot(2, 4, 6)
ax6 = sns.boxplot(x = 'otherTrees', y="RFWeighted" ,data=g_df )
ax6 = sns.stripplot(x = 'otherTrees', y="RFWeighted" ,data=g_df,  jitter=True, dodge=True, marker='o', color=".1")
ax6.set_title('RFWeighted' , fontsize=9)

ax7 = fig.add_subplot(2, 4, 7)
ax7 = sns.boxplot(x = 'otherTrees', y="GeoUnrooted" ,data=g_df )
ax7 = sns.stripplot(x = 'otherTrees', y="GeoUnrooted" ,data=g_df,  jitter=True, dodge=True, marker='o', color=".1")
ax7.set_title('GeoUnrooted' , fontsize=9)


plt.ylabel('Value')
# plt.xlabel('')
# plt.xticks([])
plt.show()


