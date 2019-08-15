from __future__ import unicode_literals
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.lines as line

import matplotlib.patches as patches
from tqdm import tqdm


plt.rcParams["font.size"] = 13

import matplotlib as mpl
import matplotlib.font_manager as font_manager
path = '/home/adyprat/.fonts/Proxima Nova Reg.ttf'
prop = font_manager.FontProperties(fname=path)
prop = font_manager.FontProperties(fname=path)

    
def plot(inputDF, height = 7, 
         levels = [], rotation = []):
    """
    Script to produce Figure 

    :param inputDF: Dataframe containing ratios to be visualized either as slanted arrows (> 1.25 or < 0.75)or a double squiggly arrow (>0.75 and < 1.25)
    :type inputDF: :obj:`pandas DataFrame`
    :param height: Height of final image
    :type height: float
    :param levels: Which columns in Level 1 to use for plotting, specified using column names
    :type levels: list
    :param rotation: specify which column titles have to be rotated
    :type rotation: list
    """
    rowNames = inputDF.index
    maxRows = len(inputDF.index)
    maxCols = len(inputDF.columns)
    pad = 2
    aspRatio = (maxCols+pad)/(maxRows+pad)
    fSize = (height*aspRatio,height)

    f = plt.figure(figsize=fSize)

    
    ax = plt.gca()
    ax.set_yticks(np.arange(0,maxRows+pad+1))
    ax.set_xticks(np.arange(0,maxCols+pad))
    

    Indices = [""] + list(rowNames)  +  [""] +  [""]
    ax.set_yticklabels(Indices, fontproperties = prop)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    aspRatio = fSize[0]/fSize[1]
    altRow = True
    for rowIdx in range(len(rowNames)):
        # For every alternative row,
        # add a rectangular grey colored patch       
        
        if altRow:
            ptch = patches.Rectangle((0,rowIdx+0.5),
                               width = maxCols + 1,
                               height = 1,  
                               edgecolor = (1,1,1),
                               facecolor = (0.9,0.9,0.9),)

        else:
            ptch = patches.Rectangle((0,rowIdx+0.5),
                   width = maxCols,
                   height = 1,  
                   edgecolor = (1,1,1),
                   facecolor = (1,1,1),)
        altRow = not(altRow)
        ax.add_artist(ptch)
        
    colCnt = 0
    colStart = colCnt

    for levl1Idx in range(len(levels)):
        ax.add_artist(ptch)
        colStart = colCnt
        plt.text(colStart+2, maxRows + 2, 
                        levels[levl1Idx], size=12, rotation=0,
                        ha="center", va="center", fontproperties = prop,
                        bbox=dict(boxstyle="round",
                        ec=(1,1,1), fc=(1,1,1)))

        colNames = inputDF[levels[levl1Idx]].columns

        for colIdx in range(len(colNames)):
                colCnt += 1
                plt.text(colStart + colIdx + 1, maxRows + 1, colNames[colIdx],
                         size=10, rotation=rotation[levl1Idx],
                         ha="center", va="center",  fontproperties = prop,
                         bbox=dict(boxstyle="round",
                         ec=(1,1,1,0),
                         fc=(1,1,1,0)))

                for rowIdx in range(len(rowNames)):

                    value  =  inputDF.loc[rowNames[rowIdx],levels[levl1Idx]][colNames[colIdx]]
                    
                    # Set font size
                    fSize = 23
                    
                    # Set desired unicode values for plotting as follows:
                    # u2197 is an angled upward poining arrow, to indicate increase
                    # u2198 is a angled downward poining arrow, to indicate decrease
                    # u21AD is a left-right wave arrow, to indicate no change
                    
                    valList = ['\u2197','\u2198','\u21AD']

                    if value >= 1.25:
                        txt = valList[0]
                        col = '#4CAF50'
                        fSize = 24

                    elif value <= 0.75:
                        txt = valList[1]
                        col = '#E53935'
                        fSize = 24

                    elif value < 1.25 and value > 0.75:
                        txt = valList[2]
                        col = '#2E4053'
                        fSize = 24

                    else:
                        # Missing values?
                        txt = '-'
                        col = '#2E4053'
                        fSize = 24

                    plt.text(colStart+colIdx+1, rowIdx+1, 
                        txt, fontsize= fSize, rotation=0,
                        ha="center", va="center", color = col,
                        bbox=dict(boxstyle="round", 
                        ec=(1,1,1,0), fc=(1,1,1,0)))

    # Clean-up plot
    ax.yaxis.set_ticks_position('none') 
    ax.xaxis.set_ticks_position('none') 
    ax.set_xticklabels([])
    
    # Show the plot
    plt.savefig('outputs/Simulated/boolNetMotifs.pdf')
