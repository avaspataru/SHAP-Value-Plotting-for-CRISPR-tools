import matplotlib.pyplot as plt
from utils import *
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

#all of the tools to be compared, including scoring methods and model types (tuscan, chop-chop)
tools = ['ssc' ,'tuscan-classification', 'wu-crispr'] #'tuscan-regression'] #,'sgrnascorer2','wu-crispr','chop-chop-xu','chop-chop-doench','chop-chop-moreno']

#the datasets these tools were ran on
datasets = ['xu']#,'doench']

labels = range(0,20)
x = np.arange(len(labels))  # the label locations
width = 0.8  # the width of the bars
totalspacer = 1 - width

spacer = totalspacer / (len(tools)*len(datasets))
barwidth = width / (len(tools)*len(datasets))

fig, ax = plt.subplots()

positions = range(0,20)

for p in range(0,21):
    plt.axvline(x=p, color="black",linewidth=0.2)

tool_positions = []
tool_labels = []
cnt = 0
for tool in tools:
    for data in datasets:
        pickleFileName = "SHAP-" + tool + "-" + data

        #get a list of tuples (avg shap val, feature name)
        featureImp = getAvgShapValues(pickleFileName)

        #creates dictionary where average[nucleotide id][position] = the average shap value of the nucleotide being at that position
        averages = getShapValsForPosFeatures(featureImp) #function def in utils

        #plot
        x_coord = x + spacer + barwidth/2 + barwidth*cnt
        pA = ax.bar(x_coord, averages[0], barwidth, color="#1f77b4")
        pC = ax.bar(x_coord, averages[1], barwidth, color="#ff7f0e")
        pG = ax.bar(x_coord, averages[2], barwidth, color="#2ca02c")
        pT = ax.bar(x_coord, averages[3], barwidth, color="#d62728")

        for xp in x_coord:
            tool_positions.append(xp)
            tool_labels.append( getShorthand(tool) )

        cnt = cnt+1

plt.xticks(tool_positions, tool_labels, rotation='vertical')

ax_t = ax.secondary_xaxis('top')
ax_t.tick_params(axis='x', direction='out', width=0, length=0, color="#D3D3D3")
ax_t.set_xticks([p + 0.5 for p in positions])
ax_t.set_xticklabels(positions)

plt.plot(range(-1,21), [0]*22, color='black', linewidth=0.2) #plot a line at 0

plt.ylabel("SHAP values")
plt.xlabel("Tools ran on XU")
plt.title("Guide positions")
plt.legend((pA[0], pC[0], pG[0], pT[0]), ('A', 'C', 'G', 'T'))


plt.show()
