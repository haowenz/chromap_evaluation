#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import math
import numpy as np
import scipy as sp
import seaborn as sns
import pandas as pd
from itertools import combinations 
from scipy.cluster import hierarchy
#import matplotlib.style as mpl
import string
import re

import csv
import matplotlib.patheffects as PathEffects


# In[144]:

# cell,sub1, sub2 is cellranger_v2. Others are from cellranger v1.2.0
data = {"chromap_cell":{}, "bwa_cell":{}, "chromap":{}, "bwa":{}, "bowtie":{}, "bwa_sub1":{}, "bwa_sub2":{}} # cellranger is for atac-seq
barcodes = {}
for platform in ["chromap_cell", "bwa_cell", "chromap", "bwa", "bowtie", "bwa_sub1", "bwa_sub2"]:
    fp = open("./" + platform + "_metadata.tsv")
    header = fp.readline()
    for line in fp:
        cols = line.rstrip().split("\t")
        barcode = cols[0]
        cluster = 0
        cellType = 0
        cellType = cols[14]
        cluster = cols[4]
        if (cluster == "NA"):
            continue
        #if (cellType == "Treg"):
        #    cellType = "CD4Tconv"
        cluster = int(cluster)
        data[platform][barcode] = {"celltype":cellType, "cluster":cluster,
                                   "UMAP_1":float(cols[-2]), "UMAP_2":float(cols[-1]),
                                  "ncount":float(cols[2]), "nfeature":float(cols[3]),
                                  "assignscore":float(cols[-4])}
        if (platform in ["chromap_cell", "bwa_cell", "bwa"]):
            if (barcode not in barcodes ):
                barcodes[barcode] = 0
            barcodes[barcode] += 1
    fp.close()


# In[175]:


# Find the consensus cell type for each cell
barcodeCellType = {}
for b in barcodes:
    if (b not in barcodes or barcodes[b] < 3):
        continue
    types = {}
    for platform in ["chromap_cell", "bwa_cell", "bwa"]: #, "bwa_sub1", "bwa_sub2"]:
        if (data[platform][b]["celltype"] not in types):
            types[data[platform][b]["celltype"]] = 0
        types[data[platform][b]["celltype"]] += 1
        #if (platform == "bwa"): #and "B" in data[platform][b]["celltype"]):
        #    types[data[platform][b]["celltype"]] += 10
    maxType = -1
    maxVote = 0
    for t in types:
        if (types[t] > maxVote):
            maxVote = types[t]
            maxType = t
    barcodeCellType[b] = maxType    


# In[176]:


# Reannotate the clusters for each platform
for platform in ["chromap_cell", "bwa_cell", "bwa"]:
    clusterCellTypeDist = {}
    origClusterCellType = {}
    for b in data[platform]:
        cluster = data[platform][b]["cluster"]
        if (cluster not in clusterCellTypeDist):
            clusterCellTypeDist[cluster] = {}
            origClusterCellType[cluster] = data[platform][b]["celltype"]
        if (b not in barcodeCellType):
            continue
        consensusCellType = barcodeCellType[b]
        if (consensusCellType not in clusterCellTypeDist[cluster]):
            clusterCellTypeDist[cluster][consensusCellType] = 0
        clusterCellTypeDist[cluster][consensusCellType] += 1
        
    clusterCellType = {}
    for cluster in clusterCellTypeDist:
        if (len(clusterCellTypeDist[cluster]) == 0):
            clusterCellType[cluster] = origClusterCellType[cluster]
        else:
            maxType = -1
            maxVote = 0
            for t in clusterCellTypeDist[cluster]:
                if (clusterCellTypeDist[cluster][t] > maxVote):
                    maxVote = clusterCellTypeDist[cluster][t]
                    maxType = t
            #f ("CD8T" in clusterCellTypeDist[cluster] 
            #    and clusterCellTypeDist[cluster]["CD8T"] * 2 > maxVote):
            #    maxType = "CD8T"
            #if ("Treg" in clusterCellTypeDist[cluster] 
            #    and clusterCellTypeDist[cluster]["CD8T"] * 2 > maxVote):
            #    maxType = "CD8T"
            if ("CD8T" in maxType):
                maxType = "CD8T"
            if ("CD4T" in maxType or "Treg" in maxType):
                maxType = "CD4T"
            if (maxType == "Neutrophils"):
                maxType = "Mono/Macro"
            clusterCellType[cluster] = maxType
    for b in data[platform]:
        data[platform][b]["newCellType"] = clusterCellType[ data[platform][b]["cluster"] ] 
    #display(clusterCellTypeDist)



# Visualization
colors = ["r", "g", "b", "orange", "cyan", "fuchsia", "pink", "slategrey"]
cellTypeToColor = {}
k = 0
tag = "newCellType"
platforms = {"chromap_cell":"Chromap", "bwa_cell":"CellRanger_v2.0.0", "bwa":"CellRanger_v1.2.0"}
for platform in platforms:
    for b in data[platform]:
        if (data[platform][b][tag] not in cellTypeToColor):
            cellTypeToColor[ data[platform][b][tag] ] = colors[k]
            k += 1
            
for platform in platforms:
    df = pd.DataFrame.from_dict(data[platform], orient="index")
    sns.set(font_scale=1.5, style="white")
    fig = plt.figure(figsize=(7, 7))
    plt.grid(False)
    #plt.set_facecolor([1, 1, 1])
    
    for cellType in cellTypeToColor:
        subdf = df.loc[ df["newCellType"] == cellType ]
        plt.plot(subdf["UMAP_1"], subdf["UMAP_2"], "o", marker="o", ms=1,
                color=cellTypeToColor[cellType])

    # Add the annotation
    for cellType in cellTypeToColor:    
        subdf = df.loc[ df[tag] == cellType ]
        if (len(subdf) == 0):
            continue
        qx = 0.5
        qy = 0.5
        if (cellType == "NK" and platform != "bwa_cell"):
            qx = 0.8
            qy = 0.25
        elif ((cellType == "B" or cellType == "CD8T") and platform == "chromap"):
            if (cellType == "CD8T"):
                qy = 0.3
            qx = 0.7
        x = np.quantile(subdf["UMAP_1"], qx) - 1
        y = np.quantile(subdf["UMAP_2"], qy) 

        txt=plt.text(x, y, cellType, fontdict={"size":15, "weight":"bold"})
        txt.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='w')])
    plt.xlabel("UMAP_1")
    plt.ylabel("UMAP_2")
    title = platforms[platform]
    #title = "Chromap" #(0.5 hour)"
    #if (platform == "bwa"):
    #    title = "CellRanger_BWA" #(26 hour)"
    #elif (platform == "bowtie"):
    #    title = "Cellranger_Bowtie2" #(29 hour)"
    plt.title(title, fontsize=25)
    
    # Add legend
    if (platform == "chromap_cell"):
        patches = []
        for cellType in cellTypeToColor:
            patches.append(plt.plot([], [], marker="o", color=cellTypeToColor[cellType], 
                                    ms=5, ls="", label=cellType)[0])
        plt.legend(handles=patches, loc="upper right")   
    plt.savefig( "maestro_" + platform + "_umap" + ".pdf", bbox_inches="tight", format="pdf" )
    plt.show()


# In[79]:


import sklearn.metrics.cluster as clusterMetrics


# In[178]:


# Compute NMIs
#alignerIdx = ["chromap", "bwa", "bowtie", "bwa_sub1", "bwa_sub2"]
alignerIdx = ["chromap_cell", "bwa_cell", "bwa"]
for i in range(3):
    for j in range(i + 1, 3):
        alignA = alignerIdx[i]
        alignB = alignerIdx[j]
        labelA = []
        labelB = []
        for b in data[alignA]:
            if (b not in data[alignB]):
                continue
            labelA.append(data[alignA][b]["newCellType"])
            labelB.append(data[alignB][b]["newCellType"])
        nmi = clusterMetrics.normalized_mutual_info_score(labelA, labelB)
        #nmi = clusterMetrics.adjusted_rand_score(labelA, labelB)
        print("NMI", alignA, alignB, nmi)
        nmi = clusterMetrics.adjusted_rand_score(labelA, labelB)
        print("ARI", alignA, alignB, nmi)


# In[179]:


# Compute NMIs
alignerIdx = ["chromap", "bwa", "bowtie", "bwa_sub1", "bwa_sub2"]
alignerIdx = ["chromap_cell", "bwa_cell", "bwa"]
for i in range(3):
    for j in range(i + 1, 3):
        alignA = alignerIdx[i]
        alignB = alignerIdx[j]
        labelA = []
        labelB = []
        for b in data[alignA]:
            if (b not in data[alignB]):
                continue
            labelA.append(data[alignA][b]["cluster"])
            labelB.append(data[alignB][b]["cluster"])
        nmi = clusterMetrics.normalized_mutual_info_score(labelA, labelB)
        #nmi = clusterMetrics.adjusted_rand_score(labelA, labelB)
        print("NMI", alignA, alignB, nmi)
        nmi = clusterMetrics.adjusted_rand_score(labelA, labelB)
        print("ARI", alignA, alignB, nmi)

# In[22]:


import csv


# In[38]:


# Compare archr clustering
archrData = {"chromap":{}, "bwa":{}, "bowtie":{}, "bwa_sub1":{}, "bwa_sub2":{}} # cellranger is for atac-seq
for platform in ["chromap", "bwa", "bowtie2", "bwa_sub1", "bwa_sub2"]:
    with open("archr/archr_" + platform + ".csv", newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        #header = reader.readline()
        for cols in reader:
            #cols = line[0].split("")
            if (cols[0] ==""):
                continue
            if (platform == "bowtie2"):
                platform = "bowtie"
            archrData[platform][cols[-1]] = {"cluster":cols[-2]}


# In[40]:


# Compute NMIs
#alignerIdx = ["chromap", "bwa", "bowtie", "bwa_sub1", "bwa_sub2"]
alignerIdx = ["chromap_cell", "bwa_cell", "bwa"]
for i in range(3):
    for j in range(i + 1, 3):
        alignA = alignerIdx[i]
        alignB = alignerIdx[j]
        labelA = []
        labelB = []
        for b in archrData[alignA]:
            if (b not in archrData[alignB]):
                continue
            labelA.append(archrData[alignA][b]["cluster"])
            labelB.append(archrData[alignB][b]["cluster"])
        nmi = clusterMetrics.normalized_mutual_info_score(labelA, labelB)
        #nmi = clusterMetrics.adjusted_rand_score(labelA, labelB)
        print("NMI", alignA, alignB, nmi)
        nmi = clusterMetrics.adjusted_rand_score(labelA, labelB)
        print("ARI", alignA, alignB, nmi)


# In[103]:


# Compute NMIs
#alignerIdx = ["chromap", "bwa", "bowtie", "bwa_sub1", "bwa_sub2"]
alignerIdx = ["chromap_cell", "bwa_cell", "bwa"]
for i in range(3):
    alignA = alignerIdx[i]
    alignB = alignerIdx[i]
    labelA = []
    labelB = []
    for b in data[alignA]:
        if (b not in archrData[alignB]):
            continue
        labelA.append(data[alignA][b]["cluster"])
        labelB.append(archrData[alignB][b]["cluster"])
    nmi = clusterMetrics.normalized_mutual_info_score(labelA, labelB)
    #nmi = clusterMetrics.adjusted_rand_score(labelA, labelB)
    print("NMI", alignA, alignB, nmi)
    nmi = clusterMetrics.adjusted_rand_score(labelA, labelB)
    print("ARI", alignA, alignB, nmi)


# In[77]:


