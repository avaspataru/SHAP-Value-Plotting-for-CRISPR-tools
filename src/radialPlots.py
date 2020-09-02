import matplotlib.pyplot as plt
from utils import *
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import pandas as pd
from math import pi



# scales 4 separated arrays as if they are a single array
def scaleAll(l1,l2,l3,l4):
    minim = min(l1.min(),l2.min(),l3.min(),l4.min())
    maxim = max(l1.max(),l2.max(),l3.max(),l4.max())

    if minim == maxim: # no scaling to be done
        return l1,l2,l3,l4

    sl1 = (l1 - minim) / (maxim-minim)
    sl2 = (l2 - minim) / (maxim-minim)
    sl3 = (l3 - minim) / (maxim-minim)
    sl4 = (l4 - minim) / (maxim-minim)

    return sl1,sl2,sl3,sl4


# adds a single model to the radial plot as a line
def add_plot(ax,values, angles, tool, labelFlag):
    # Draw ylabels
    ax.set_rlabel_position(0)
    # Plot data
    if labelFlag:
        ax.plot(angles, values, linewidth=1, linestyle='solid', label=tool)
    else:
        ax.plot(angles, values, linewidth=1, linestyle='solid')

    # Fill area
    ax.fill(angles, values, alpha=0.1)


#all of the tools to be compared, including scoring methods and model types (tuscan, chop-chop)
tools = ['wu-crispr', 'tuscan-regression' ,'sgrnascorer2','tuscan-classification','chop-chop-xu','chop-chop-doench','ssc']

#the datasets these tools were ran on
datasets = ['xu']

#setting up for plotting
positions = range(0,20)

# number of variable
categories=list(range(0,20))
N = len(categories)

# What will be the angle of each axis in the plot? (we divide the plot / number of variable)
angles = [n / float(N) * 2 * pi for n in range(N)]
angles += angles[:1]

# Initialise the spider plot
fig = plt.figure(figsize=(10,10))
fig.suptitle("Nucleotide preference across models in each position", fontsize=14)

# Adenine subplot initialization
ax1 = plt.subplot(221, polar=True)
plt.xticks(angles[:-1], categories, color='grey', size=8)
plt.yticks(np.arange(0,1.1,step=0.25), ["0",".25",".5",".75","1"], color="grey", size=7)
plt.ylim(0,1)
plt.title("A",loc='left')

# Cythosine subplot initialization
ax2 = plt.subplot(222, polar=True)
plt.xticks(angles[:-1], categories, color='grey', size=8)
plt.yticks(np.arange(0,1.1,step=0.25), ["0",".25",".5",".75","1"], color="grey", size=7)
plt.ylim(0,1)
plt.title("C",loc='left')

# Guanine subplot initialization
ax3 = plt.subplot(223, polar=True)
plt.xticks(angles[:-1], categories, color='grey', size=8)
plt.yticks(np.arange(0,1.1,step=0.25), ["0",".25",".5",".75","1"], color="grey", size=7)
plt.ylim(0,1)
plt.title("G",loc='left')

# Thymine subplot initialization
ax4 = plt.subplot(224, polar=True)
plt.xticks(angles[:-1], categories, color='grey', size=8)
plt.yticks(np.arange(0,1.1,step=0.25), ["0",".25",".5",".75","1"], color="grey", size=7)
plt.ylim(0,1)
plt.title("T",loc='left')


for tool in tools:
    for data in datasets:
        #load data from pickle
        pickleFileName = "SHAP-" + tool + "-" + data
        #get a list of tuples (avg shap val, feature name)
        featureImp = getAvgShapValues(pickleFileName)
        #creates dictionary where average[nucleotide id][position] = the average shap value of the nucleotide being at that position
        averages = getShapValsForPosFeatures(featureImp) #function def in utils


        pA, pC, pG, pT = scaleAll(averages[0],averages[1],averages[2],averages[3])

        # We are going to plot the first line of the data frame.
        # But we need to repeat the first value to close the circular graph:

        values = [ i for i in pA]
        values.append(pA[:1])
        add_plot(ax1,values,angles, tool, True)

        values = [ i for i in pT]
        values.append(pT[:1])
        add_plot(ax4,values,angles, tool, False)

        values = [ i for i in pG]
        values.append(pG[:1])
        add_plot(ax3,values,angles, tool, False)

        values = [ i for i in pC]
        values.append(pC[:1])
        add_plot(ax2,values,angles, tool, False)



fig.legend(loc='upper right')
plt.savefig('radial.png',bbox_inches="tight")
plt.show()
