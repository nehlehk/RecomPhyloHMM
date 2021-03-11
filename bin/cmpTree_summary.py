import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import csv
import matplotlib.pyplot as plt
import seaborn as sns





Quartet_hmm = []
PathDiffernce_hmm = []
RF_hmm = []
MatchingSplit_hmm = []
UMAST_hmm = []
RFWeighted_hmm = []
GeoUnrooted_hmm = []

Quartet_rax = []
PathDiffernce_rax = []
RF_rax = []
MatchingSplit_rax = []
UMAST_rax = []
RFWeighted_rax = []
GeoUnrooted_rax = []

Quartet_b = []
PathDiffernce_b = []
RF_b = []
MatchingSplit_b = []
UMAST_b = []
RFWeighted_b = []
GeoUnrooted_b = []
x = []

files = glob.glob("/home/nehleh/Documents/myTemp/*.txt")
for id,f in enumerate(files):
    # print(f)
    csv = pd.read_csv(f, sep='\t', lineterminator='\n',names=None)
    x.append(id)
    # print(csv)
    Quartet_hmm.append(csv['Quartet'][0])
    PathDiffernce_hmm.append(csv['PathDiffernce'][0])
    RF_hmm.append(csv['RF(0.5)'][0])
    MatchingSplit_hmm.append(csv['MatchingSplit'][0])
    UMAST_hmm.append(csv['UMAST'][0])
    RFWeighted_hmm.append(csv['RFWeighted(0.5)'][0])
    GeoUnrooted_hmm.append(csv['GeoUnrooted'][0])

    Quartet_rax.append(csv['Quartet'][1])
    PathDiffernce_rax.append(csv['PathDiffernce'][1])
    RF_rax.append(csv['RF(0.5)'][1])
    MatchingSplit_rax.append(csv['MatchingSplit'][1])
    UMAST_rax.append(csv['UMAST'][1])
    RFWeighted_rax.append(csv['RFWeighted(0.5)'][1])
    GeoUnrooted_rax.append(csv['GeoUnrooted'][1])

    Quartet_b.append(csv['Quartet'][2])
    PathDiffernce_b.append(csv['PathDiffernce'][2])
    RF_b.append(csv['RF(0.5)'][2])
    MatchingSplit_b.append(csv['MatchingSplit'][2])
    UMAST_b.append(csv['UMAST'][2])
    RFWeighted_b.append(csv['RFWeighted(0.5)'][2])
    GeoUnrooted_b.append(csv['GeoUnrooted'][2])


final = pd.DataFrame({'Quartet_hmm': Quartet_hmm, 'PathDiffernce_hmm': PathDiffernce_hmm ,'RF_hmm':RF_hmm, 'MatchingSplit_hmm':MatchingSplit_hmm, 'UMAST_hmm':UMAST_hmm, 'RFWeighted_hmm':RFWeighted_hmm , 'GeoUnrooted_hmm':GeoUnrooted_hmm,
                      'Quartet_rax': Quartet_rax, 'PathDiffernce_rax': PathDiffernce_rax ,'RF_rax':RF_rax, 'MatchingSplit_rax':MatchingSplit_rax, 'UMAST_rax':UMAST_rax, 'RFWeighted_rax':RFWeighted_rax , 'GeoUnrooted_rax':GeoUnrooted_rax,
                      'Quartet_b': Quartet_b, 'PathDiffernce_b': PathDiffernce_b ,'RF_b':RF_b, 'MatchingSplit_b':MatchingSplit_b, 'UMAST_b':UMAST_b, 'RFWeighted_b':RFWeighted_b , 'GeoUnrooted_b':GeoUnrooted_b})


print(final)




fig = plt.figure(figsize=(16,8))

# ax1 = sns.boxplot( data=final ,hue = 'prob'  )
# ax1 = sns.stripplot(data=final,  jitter=True, dodge=True, marker='o', color=".1")


# ax = fig.add_subplot(2, 4, 1)
# ax.set_title("Quartet")
# ax.plot(x,Quartet_hmm ,label= "phyloHMM")
# ax.plot(x,Quartet_b,label= "Beast")
# ax.plot(x,Quartet_rax,label= "RaXML")
# ax.legend(loc = 'best')
#
# ax = fig.add_subplot(2, 4, 2)
# ax.set_title("PathDiffernce")
# ax.plot(x,PathDiffernce_hmm, marker='.' ,label= "phyloHMM")
# ax.plot(x,PathDiffernce_b, marker='.' ,label= "Beast")
# ax.plot(x,PathDiffernce_rax, marker='.' ,label= "RaXML")
# ax.legend(loc = 'best')
#
# ax = fig.add_subplot(2, 4, 3)
# ax.set_title("RF")
# ax.plot(x,RF_hmm, marker='.' ,label= "phyloHMM")
# ax.plot(x,RF_b, marker='.' ,label= "Beast")
# ax.plot(x,RF_rax, marker='.' ,label= "RaXML")
# ax.legend(loc = 'best')
#
# ax = fig.add_subplot(2, 4, 4)
# ax.set_title("MatchingSplit")
# ax.plot(x,MatchingSplit_hmm, marker='.' ,label= "phyloHMM")
# ax.plot(x,MatchingSplit_b, marker='.' ,label= "Beast")
# ax.plot(x,MatchingSplit_rax, marker='.' ,label= "RaXML")
# ax.legend(loc = 'best')
#
#
# ax = fig.add_subplot(2, 4, 5)
# ax.set_title("UMAST")
# ax.plot(x,UMAST_hmm, marker='.' ,label= "phyloHMM")
# ax.plot(x,UMAST_b, marker='.' ,label= "Beast")
# ax.plot(x,UMAST_rax, marker='.' ,label= "RaXML")
# ax.legend(loc = 'best')
#
# ax = fig.add_subplot(2, 4, 6)
# ax.set_title("RFWeighted")
# ax.plot(x,RFWeighted_hmm, marker='.' ,label= "phyloHMM")
# ax.plot(x,RFWeighted_b, marker='.' ,label= "Beast")
# ax.plot(x,RFWeighted_rax, marker='.' ,label= "RaXML")
# ax.legend(loc = 'best')
#
#
# ax = fig.add_subplot(2, 4, 7)
# ax.set_title("GeoUnrooted")
# ax.plot(x,GeoUnrooted_hmm, marker='.' ,label= "phyloHMM")
# ax.plot(x,GeoUnrooted_b, marker='.' ,label= "Beast")
# ax.plot(x,GeoUnrooted_rax, marker='.' ,label= "RaXML")
# ax.legend(loc = 'best')

# plt.savefig("cmpTree.jpeg")



hmm = [np.average(Quartet_hmm), np.average(PathDiffernce_hmm) ,np.average(RF_hmm) ,np.average(MatchingSplit_hmm), np.average(UMAST_hmm)]
raxml = [np.average(Quartet_rax), np.average(PathDiffernce_rax) ,np.average(RF_rax) ,np.average(MatchingSplit_rax), np.average(UMAST_rax)]
beast = [np.average(Quartet_b), np.average(PathDiffernce_b) ,np.average(RF_b) ,np.average(MatchingSplit_b), np.average(UMAST_b)]
index = ['Quartet', 'PathDiffernce', 'RF','MatchingSplit', 'UMAST']


hmm2 = [np.average(RFWeighted_hmm) ,np.average(GeoUnrooted_hmm)]
raxml2 = [ np.average(RFWeighted_rax) ,np.average(GeoUnrooted_rax)]
beast2 = [ np.average(RFWeighted_b) ,np.average(GeoUnrooted_b)]
index2 = [ 'RFWeighted', 'GeoUnrooted']



fig2 = plt.figure(figsize=(8,8))
ax = fig2.add_subplot(2, 1, 1)
X = np.arange(len(index))
ax.set_title("The distances between clonaltree and the other trees")
ax.bar(X + 0.00, hmm, color = '#0F52BA', width = 0.25 , label = 'PhyloHMM')
ax.bar(X + 0.25, raxml, color = '#6593F5', width = 0.25 , label = 'RaXml')
ax.bar(X + 0.50, beast, color = '#73C2FB', width = 0.25 , label = 'Beast')
ax.set_xticks(X + 0.25)
ax.set_xticklabels(index)
ax.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.3)
ax.legend()


ax2 = fig2.add_subplot(2, 1, 2)
X2 = np.arange(len(index2))
ax2.bar(X2 + 0.00, hmm2, color = '#0F52BA', width = 0.25 , label = 'PhyloHMM')
ax2.bar(X2 + 0.25, raxml2, color = '#6593F5', width = 0.25 , label = 'RaXml')
ax2.bar(X2 + 0.50, beast2, color = '#73C2FB', width = 0.25 , label = 'Beast')
ax2.set_xticks(X2 + 0.25)
ax2.set_xticklabels(index2)
ax2.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.3)
ax2.legend()


fig2.tight_layout()
plt.show()
