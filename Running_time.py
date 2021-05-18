#!/usr/bin/env python
# coding: utf-8

# In[1]:


import re
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import math
import scipy as sp


# In[13]:


runtime = {"ChIP-seq":{"Chromap":[3.5, 5], "BWA":[64, 100]}, 
        "10X scATAC-seq": {"Chromap":[23, 28.5], "BWA": [359, 466]}, #"BWA_cellranger1.2.0":[416, 2001]},
        "Hi-C": {"Chromap":[144, 164], "BWA":[1364, 2143]}}

k = 0
sns.set(font_scale=2, style="white")
fig, axes = plt.subplots(1, len(runtime), figsize=(16, 10))
topbar = plt.Rectangle((0,0),1,1,fc="red", edgecolor = 'none')
bottombar = plt.Rectangle((0,0),1,1,fc='blue',  edgecolor = 'none')
for title in ["ChIP-seq", "Hi-C", "10X scATAC-seq"]:
    methodv = ["Chromap", "Chromap", "BWA", "BWA"]
    timev = []
    for m in ["Chromap", "BWA"]:
        timev.extend(runtime[title][m])
    # Actually other is total running time for stacked bar plot
    categoryv = ["Alignment", "Other", "Alignment", "Other"] 
    
    df = pd.DataFrame()
    df["Method"] = methodv
    df["Time"] = timev
    df["Category"] = categoryv
    
    ax = axes[k]
    snsFig = sns.barplot(x="Method", y="Time", data=df.loc[df["Category"] == "Other"], ax=ax, color="red")
    snsFig = sns.barplot(x="Method", y="Time", data=df.loc[df["Category"] == "Alignment"], ax=ax, color="blue")
    #ax.set(yscale="log")
    if (k > 0):
        ax.set(ylabel="")
    else:
        ax.set(ylabel="Time (minute)")
    ax.set(title=title, xlabel="")
    ax.plot([0.5, 0.5], [runtime[title]["Chromap"][1], runtime[title]["BWA"][1]], c="black")
    ax.plot([0.4, 0.5], [runtime[title]["Chromap"][1], runtime[title]["Chromap"][1]], c="black")
    ax.plot([0.5, 0.6], [runtime[title]["BWA"][1], runtime[title]["BWA"][1]], c="black")
    ax.text(0.2, runtime[title]["BWA"][1] / 2, "%dX"%(round(runtime[title]["BWA"][1]/runtime[title]["Chromap"][1])),
           bbox={'facecolor':'white','alpha':1,'edgecolor':'none','pad':1})
    #if (runtime[title]["BWA"][1] > 1800):
    #    ax.axhline(1440, c="orange", ls="--")
    #    ax.text(-0.25, 1440 + 50, "1 day", c="orange")
    #s="%d"%(runtime[title]["Chromap"][1]) 
    #ax.text(x=-0.25, y=runtime[title]["Chromap"][1], s="%d"%(runtime[title]["Chromap"][1]))
    #ax.text(x=0.75, y=runtime[title]["BWA"][1], s="%d"%(runtime[title]["BWA"][1]))
    if (k == 2):
        ax.legend([bottombar, topbar], ['Mapping', 'Other\n(Dedup,...)'], loc="upper left", bbox_to_anchor=(1.01, 0.5))
    k += 1 
plt.tight_layout() 
plt.savefig("chromap_bwa_runtime.pdf", bbox_inches="tight", format="pdf" )
plt.show()

