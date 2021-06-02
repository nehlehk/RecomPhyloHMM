import pandas as pd
import matplotlib.pyplot as plt
import argparse
import seaborn as sns
import numpy as np




parser = argparse.ArgumentParser(description='''You did not specify any parameters.''')
parser.add_argument('-c', "--cmpTreesResult", type=str, help='cmpTreesResult')
args = parser.parse_args()

cmpTrees = args.cmpTreesResult

cmpTrees = '/home/nehleh/PhyloCode/Result/10nodes-20samples-len10000/short_philobacteria/Summary_Results/all_cmpTrees.result'


f = open(cmpTrees, "r")
df = pd.read_csv(f,sep='\t', names=['No', 'clonalTree','otherTrees', 'Quartet' , 'PathDiffernce' , 'RF' , 'MatchingSplit' , 'UMAST' , 'RFWeighted' , 'GeoUnrooted' ] ,header=None)
new_df= df.loc[df['otherTrees'].between('1' , '6' , inclusive=True )].apply(pd.to_numeric)
new_df['otherTrees'].replace(to_replace=[1,2,3,4,5,6],value=['PB_4' , 'PB_2','RaxML','Beast' ,'CFML','Gubbins',] ,inplace=True)

# print(new_df)

g_df = new_df.loc[(new_df['otherTrees'] != 'Gubbins') ]

c_df = new_df.loc[(new_df['otherTrees'] != 'Beast') ]


# both_df = new_df.loc[(new_df['otherTrees'].isin(['PhiloBacteria','CFML'])) ]
# print(both_df)

fig = plt.figure(figsize=(14,14))

ax1 = fig.add_subplot(4, 2, 1)
ax1 = sns.boxplot(x = 'otherTrees', y="Quartet" ,data=new_df)
ax1 = sns.stripplot(x = 'otherTrees', y="Quartet" ,data=new_df,  jitter=True, dodge=True, marker='.', color=".1")
ax1.set_title('(a)-Quartet' , fontsize=11)
plt.xlabel('')


ax2 = fig.add_subplot(4, 2, 2)
ax2 = sns.boxplot(x = 'otherTrees', y="PathDiffernce" ,data=new_df )
ax2 = sns.stripplot(x = 'otherTrees', y="PathDiffernce" ,data=new_df,  jitter=True, dodge=True, marker='.', color=".1")
ax2.set_title('(b)-PathDiffernce' , fontsize=11)
plt.xlabel('')

ax5 = fig.add_subplot(4, 2, 3)
ax5 = sns.boxplot(x = 'otherTrees', y="UMAST" ,data=new_df )
ax5 = sns.stripplot(x = 'otherTrees', y="UMAST" ,data=new_df,  jitter=True, dodge=True, marker='.', color=".1")
ax5.set_title('(c)-UMAST' , fontsize=11)
plt.xlabel('')

ax3 = fig.add_subplot(4, 2, 4)
ax3 = sns.boxplot(x = 'otherTrees', y="RF" ,data=new_df )
ax3 = sns.stripplot(x = 'otherTrees', y="RF" ,data=new_df,  jitter=True, dodge=True, marker='.', color=".1")
ax3.set_title('(d)-Robinson-Foulds(RF)' , fontsize=11)
plt.xlabel('')

ax6 = fig.add_subplot(4, 2, 5)
ax6 = sns.boxplot(x = 'otherTrees', y="RFWeighted" ,data=g_df )
ax6 = sns.stripplot(x = 'otherTrees', y="RFWeighted" ,data=g_df,  jitter=True, dodge=True, marker='.', color=".1")
ax6.set_title('(e)-RFWeighted' , fontsize=11)
plt.xlabel('')

ax4 = fig.add_subplot(4, 2, 6)
ax4 = sns.boxplot(x = 'otherTrees', y="MatchingSplit" ,data=new_df )
ax4 = sns.stripplot(x = 'otherTrees', y="MatchingSplit" ,data=new_df,  jitter=True, dodge=True, marker='.', color=".1")
ax4.set_title('(f)-MatchingSplit' , fontsize=11)
plt.xlabel('')


ax7 = fig.add_subplot(4, 2, 7)
ax7 = sns.boxplot(x = 'otherTrees', y="GeoUnrooted" ,data=g_df )
ax7 = sns.stripplot(x = 'otherTrees', y="GeoUnrooted" ,data=g_df,  jitter=True, dodge=True, marker='.', color=".1")
ax7.set_title('(g)-GeoUnrooted' , fontsize=11)
plt.xlabel('')

fig.tight_layout(pad=3.0)



# plt.ylabel('Value')
plt.xlabel('')
# plt.xticks([])
plt.show()


fig1 = plt.figure(figsize=(10,10))
ax1 = fig1.add_subplot(1, 1, 1)
ax1 = sns.boxplot(x = 'otherTrees', y="Quartet" ,data=new_df)
ax1 = sns.stripplot(x = 'otherTrees', y="Quartet" ,data=new_df,  jitter=True, dodge=True, marker='o', color=".1")
ax1.set_title('Quartet' , fontsize=12)
plt.savefig("treecmp_Quartet.jpeg")


fig2 = plt.figure(figsize=(10,10))
ax2 = fig2.add_subplot(1, 1, 1)
ax2 = sns.boxplot(x = 'otherTrees', y="PathDiffernce" ,data=new_df )
ax2 = sns.stripplot(x = 'otherTrees', y="PathDiffernce" ,data=new_df,  jitter=True, dodge=True, marker='o', color=".1")
ax2.set_title('PathDiffernce' , fontsize=9)
plt.savefig("treecmp_PathDiffernce.jpeg")


fig3 = plt.figure(figsize=(10,10))
ax3 = fig3.add_subplot(1, 1, 1)
ax3 = sns.boxplot(x = 'otherTrees', y="RF" ,data=new_df )
ax3 = sns.stripplot(x = 'otherTrees', y="RF" ,data=new_df,  jitter=True, dodge=True, marker='o', color=".1")
ax3.set_title('RF' , fontsize=9)
plt.savefig("treecmp_RF.jpeg")


fig4 = plt.figure(figsize=(10,10))
ax4 = fig4.add_subplot(1, 1, 1)
ax4 = sns.boxplot(x = 'otherTrees', y="MatchingSplit" ,data=new_df )
ax4 = sns.stripplot(x = 'otherTrees', y="MatchingSplit" ,data=new_df,  jitter=True, dodge=True, marker='o', color=".1")
ax4.set_title('MatchingSplit' , fontsize=9)
plt.savefig("treecmp_MatchingSplit.jpeg")


fig5 = plt.figure(figsize=(10,10))
ax5 = fig5.add_subplot(1, 1, 1)
ax5 = sns.boxplot(x = 'otherTrees', y="UMAST" ,data=new_df )
ax5 = sns.stripplot(x = 'otherTrees', y="UMAST" ,data=new_df,  jitter=True, dodge=True, marker='o', color=".1")
ax5.set_title('UMAST' , fontsize=9)
plt.savefig("treecmp_UMAST.jpeg")


fig6 = plt.figure(figsize=(10,10))
ax6 = fig6.add_subplot(1, 1, 1)
ax6 = sns.boxplot(x = 'otherTrees', y="RFWeighted" ,data=g_df )
ax6 = sns.stripplot(x = 'otherTrees', y="RFWeighted" ,data=g_df,  jitter=True, dodge=True, marker='o', color=".1")
ax6.set_title('RFWeighted' , fontsize=9)
plt.savefig("treecmp_RFWeighted.jpeg")


fig7 = plt.figure(figsize=(10,10))
ax7 = fig7.add_subplot(1, 1, 1)
ax7 = sns.boxplot(x = 'otherTrees', y="GeoUnrooted" ,data= g_df )
ax7 = sns.stripplot(x = 'otherTrees', y="GeoUnrooted" ,data= g_df,  jitter=True, dodge=True, marker='o', color=".1")
ax7.set_title('GeoUnrooted' , fontsize=9)
plt.savefig("treecmp_GeoUnrooted.jpeg")


plt.ylabel('Value')
plt.xlabel('')
plt.xticks([])

# plt.show()


# fig1 = plt.figure(figsize=(10,10))
# ax1 = fig1.add_subplot(1, 1, 1)
# ax1 = sns.boxplot(x = 'otherTrees', y="Quartet" ,data=c_df)
# ax1 = sns.stripplot(x = 'otherTrees', y="Quartet" ,data=c_df,  jitter=True, dodge=True, marker='o', color=".1")
# ax1.set_title('Quartet' , fontsize=12)
# plt.savefig("treecmp_Quartet.jpeg")
#
#
# fig2 = plt.figure(figsize=(10,10))
# ax2 = fig2.add_subplot(1, 1, 1)
# ax2 = sns.boxplot(x = 'otherTrees', y="PathDiffernce" ,data=c_df )
# ax2 = sns.stripplot(x = 'otherTrees', y="PathDiffernce" ,data=c_df,  jitter=True, dodge=True, marker='o', color=".1")
# ax2.set_title('PathDiffernce' , fontsize=9)
# plt.savefig("treecmp_PathDiffernce.jpeg")
#
#
# fig3 = plt.figure(figsize=(10,10))
# ax3 = fig3.add_subplot(1, 1, 1)
# ax3 = sns.boxplot(x = 'otherTrees', y="RF" ,data=c_df )
# ax3 = sns.stripplot(x = 'otherTrees', y="RF" ,data=c_df,  jitter=True, dodge=True, marker='o', color=".1")
# ax3.set_title('RF' , fontsize=9)
# plt.savefig("treecmp_RF.jpeg")
#
#
# fig4 = plt.figure(figsize=(10,10))
# ax4 = fig4.add_subplot(1, 1, 1)
# ax4 = sns.boxplot(x = 'otherTrees', y="MatchingSplit" ,data=c_df )
# ax4 = sns.stripplot(x = 'otherTrees', y="MatchingSplit" ,data=c_df,  jitter=True, dodge=True, marker='o', color=".1")
# ax4.set_title('MatchingSplit' , fontsize=9)
# plt.savefig("treecmp_MatchingSplit.jpeg")
#
#
# fig5 = plt.figure(figsize=(10,10))
# ax5 = fig5.add_subplot(1, 1, 1)
# ax5 = sns.boxplot(x = 'otherTrees', y="UMAST" ,data=c_df )
# ax5 = sns.stripplot(x = 'otherTrees', y="UMAST" ,data=c_df,  jitter=True, dodge=True, marker='o', color=".1")
# ax5.set_title('UMAST' , fontsize=9)
# plt.savefig("treecmp_UMAST.jpeg")
#
#
# fig6 = plt.figure(figsize=(10,10))
# ax6 = fig6.add_subplot(1, 1, 1)
# ax6 = sns.boxplot(x = 'otherTrees', y="RFWeighted" ,data=c_df )
# ax6 = sns.stripplot(x = 'otherTrees', y="RFWeighted" ,data=c_df,  jitter=True, dodge=True, marker='o', color=".1")
# ax6.set_title('RFWeighted' , fontsize=9)
# plt.savefig("treecmp_RFWeighted.jpeg")
#
#
# fig7 = plt.figure(figsize=(10,10))
# ax7 = fig7.add_subplot(1, 1, 1)
# ax7 = sns.boxplot(x = 'otherTrees', y="GeoUnrooted" ,data= c_df )
# ax7 = sns.stripplot(x = 'otherTrees', y="GeoUnrooted" ,data= c_df,  jitter=True, dodge=True, marker='o', color=".1")
# ax7.set_title('GeoUnrooted' , fontsize=9)
# plt.savefig("treecmp_GeoUnrooted.jpeg")
#
#
# plt.ylabel('Value')


