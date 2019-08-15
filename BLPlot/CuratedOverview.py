from __future__ import unicode_literals
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.lines as line

import matplotlib.patches as patches
#%matplotlib inline

from tqdm import tqdm

plt.rcParams["font.size"] = 12

import matplotlib as mpl
import matplotlib.font_manager as font_manager
## We recommend using the font Proxima Nova Reg
## After downloading the .ttf file, it can be configured to
## be used in this script by uncommenting the following lines

# path = '/path/to/Proxima Nova Reg.ttf'
# prop = font_manager.FontProperties(fname=path)

def plot(inputDF, height = 7, randValues = [], shape = [], 
            palettes = [], text = [], levels = [], rotation = []):
    """

    :param inputDF: Multilevel dataframe to be plotted
    :type inputDF: :obj:`pandas DataFrame`

    :param height: Height of image to be generated
    :type height: float

    :param randValues: Cutoffs below which a pre-defined shape is set, useful for indicating less-than-significant values
    :type randValues: list

    :param shape: Shape to be drawn, choose from 's' (square), 'rs' (rounded square), 'w' (wedge), 'b' (fixed width rectangle), 'f' (full width, fixed color square), 'text'
    :type shape: list

    :param palettes: Color palette. Define this object using sns.color_palette()
    :type palettes: list

    :param text: Specify if column is to be treated as text
    :type text: list

    :param levels: Which columns to plot, specified using column names
    :type levels: list

    :param rotation: Whether to rotate the text of a column name.
    :type rotation: list

    """
    
    levls = levels
    print(levls)
    rowNames = inputDF.index
    maxRows = len(inputDF.index)
    maxCols = len(inputDF.columns)
    pad = 2
    aspRatio = (maxCols+pad)/(maxRows+pad)
    fSize = (height*aspRatio+0.5,height)

    f = plt.figure(figsize=fSize)

    
    ax = plt.gca()
    ax.set_yticks(np.arange(0,maxRows+pad))
    ax.set_xticks(np.arange(0,maxCols+pad))
    

    Indices = [""] + list(rowNames)  +  [""]
    ax.set_yticklabels(Indices, fontproperties = prop)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    aspRatio = fSize[0]/fSize[1]
    print(aspRatio)
    alt = True
    for rowIdx in range(len(rowNames)):
        if alt:
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
        alt = not(alt)
        ax.add_artist(ptch)

    colCnt = 0
    for levlIdx in range(len(levls)):
        randValue = randValues[levlIdx]
        colStart = colCnt
        plt.text(colStart+2, maxRows + 2, 
                        levls[levlIdx], size=12, rotation=0,
                        ha="center", va="center", fontproperties = prop,
                        bbox=dict(boxstyle="round",
                        ec=(1,1,1), fc=(1,1,1)))
        
        colNames = inputDF[levels[levlIdx]].columns

        for colIdx in range(len(colNames)):
                colCnt += 1

                plt.text(colStart + colIdx + 1, maxRows + 1, colNames[colIdx],
                         size=10, rotation=rotation[levlIdx],
                         ha="center", va="center",  fontproperties = prop,
                         bbox=dict(boxstyle="round",
                         ec=(1,1,1,0),
                         fc=(1,1,1,0)))

                for rowIdx in range(len(rowNames)):

                    value  =  inputDF.loc[rowNames[rowIdx],levls[levlIdx]][colNames[colIdx]]
                    if text[levlIdx]:
                        fSize = 23

                        if value == 'Y':
                            txt = '\u2713'
                            col = '#4CAF50'
                        elif value == 'N':
                            txt = '\u2717'
                            col = '#E53935'
                        else:
                            valList = ['\u2B24','\u25EF','\u25B0','\u25B1','\u25AE','\u25AF','\u2593','\u2591']
                            
                            if value < 1:
                                txt = valList[1]+valList[1]+valList[1]
                            elif value < 2:
                                txt = valList[0]+valList[1]+valList[1]
                            elif value < 3:
                                txt = valList[0]+valList[0]+valList[1]
                            else:
                                txt = valList[0]+valList[0]+valList[0]

                            col = '#2E4053'
                            fSize = 9
                        plt.text(colStart+colIdx+1, rowIdx+1, 
                            txt, fontsize= fSize, rotation=0,
                            ha="center", va="center", color = col,
                            bbox=dict(boxstyle="round", 
                            ec=(1,1,1,0), fc=(1,1,1,0)))
                    
                    if shape[levlIdx] != 'text':
                        #randValue = inputDF.loc['dummy',levls[levlIdx]][colNames[colIdx]]
                        Oldvalue = 1
                        if value < randValue:
                            col = '#2E4053'
                        elif value >= 1:
                            Oldvalue = value

                            value = min(value,5)
                            col = palettes[levlIdx][int(np.floor((value/5)*10))]

                            value = 1


                                
                        elif np.isnan(value):
                            value = 0
                        else:
                            col = palettes[levlIdx][int(np.floor((value)*10))]


                        if shape[levlIdx] == 'c':
                            circle1=patches.Circle((colStart+colIdx+1,rowIdx+1),
                                               radius = np.sqrt(value)/2.5,
                                               #height = value,    
                                               facecolor=col,
                                               edgecolor = 'k',)

                        elif shape[levlIdx] == 's':
                            value = value*0.8
                            circle1=patches.Rectangle((colStart+colIdx+1-(value/2),rowIdx+1-(value)/2),
                                       width = value,
                                       height = value,    
                                       facecolor=col,
                                       edgecolor = 'k',)

                        elif shape[levlIdx] == 'rs':
                            value = value*0.8
                            if value <= 0.15:
                                value = 0.15
                            boxPad = 0.075
                            newVal = value - boxPad*2
                            circle1=patches.FancyBboxPatch((colStart+colIdx+1-(newVal/2),rowIdx+1-(newVal)/2),
                                       width = newVal,
                                       height = newVal,    
                                       facecolor=col,
                                       edgecolor = 'k',
                                       boxstyle=patches.BoxStyle("Round", pad=boxPad))
                                        #patches.BoxStyle("Round4", pad=0.05))
                                
                        elif shape[levlIdx] == 'w':
                            circle1=patches.Wedge((colStart+colIdx+1,rowIdx+1),
                                       r = 0.4,
                                       theta1 = 0,
                                       theta2 = round(value*360,2),    
                                       facecolor = col,
                                       edgecolor = 'k',)
                        elif shape[levlIdx] == 'b':
                            circle1=patches.Rectangle((colStart+colIdx+0.6,rowIdx+0.65),
                                       width = 0.75,
                                       height = value,    
                                       facecolor=col,
                                       edgecolor = 'k',)
                            
                        elif shape[levlIdx] == 'f':
                            # flat color
                            circle1=patches.Rectangle((colStart+colIdx+1,rowIdx+0.6),
                                       width = 1,
                                       height = 1,    
                                       facecolor = col,
                                       edgecolor = col,)

                            

                        ax.add_artist(circle1)
                        if Oldvalue > 1:
                            plt.text(colStart+colIdx+1, rowIdx+1, 
                            round(Oldvalue,1), fontsize= 10, rotation=0,
                            ha="center", va="center", color = 'white',
                            bbox=dict(boxstyle="round", 
                            ec=(1,1,1,0), fc=(1,1,1,0)))
                            

    ax.yaxis.set_ticks_position('none') 
    ax.xaxis.set_ticks_position('none') 
    ax.set_xticklabels([])

    #plt.grid()
#PTData = pd.read_csv("outputs/dyn-BF/dyn-BF-AUROCscores.csv", header = 0, index_col = 0)

