# -*- coding: utf-8 -*-
"""
Created on Thu May 14 09:37:11 2020

@author: René Wilbers
"""

#user settings
morphdir = r'D:\FSNeuron\SWC' # directory of morphologies. Put only SWC & ASC files in this directory!

#morphology properties

import os
import neurom as nm
from neurom import viewer
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from morph_analysis_tools import *

files=os.listdir(morphdir)


for index,file in enumerate(files):
    fnsplit=file[:-4].split('_')
    nrnid =int(fnsplit[2])
    species=fnsplit[0][:-2]
    layer=fnsplit[1]
    
    nrn = nm.load_neuron(os.path.join(morphdir, file))
    #get neuron stats
    axon=get_sections(nrn, 'axon')
    dendBasal=get_sections(nrn, 'basal_dendrite')
    dendApical=get_sections(nrn, 'apical_dendrite')
    dend = dendBasal + dendApical
    
    leafs = get_leafs(dend)
    for leaf in leafs:
        DLL, DBL_total, DBL_mean, order, pathlength, branchlengths = get_leaf_features(leaf)
        rowdict={'id':nrnid,'species':species,'layer':layer,'DLL':DLL, 'DBL_total':DBL_total, 'DBL_mean':DBL_mean, 'order':order, 'pathlength':pathlength, 'branchlengths':branchlengths}
        if index==0:
            rowdict["branchlengths"]=[rowdict["branchlengths"]]
            df = pd.DataFrame(rowdict, index=[0])
        else:
            df=df.append(rowdict, ignore_index=True)

df.to_pickle('D:\Human inhibition\FS paper\Apr2020\Figure 2\df_per_terminal.pkl')

plt.figure()
for specie, color in zip(['Human', 'Mouse'][::-1], ['C1', 'C0']):
    data=df.loc[df['species']==specie]
    plt.scatter(data['DBL_mean'], data['DLL'], color=color)
plt.xlabel('Mean length of branching segments (\u03BCm)')
plt.ylabel('Terminal segment length (\u03BCm)')
plt.legend(['Mouse', 'Human'])
plt.savefig(os.path.join(r'D:\Human inhibition\FS paper\Apr2020\Figure 2\DBLvsDLL_per_terminal.pdf'), dpi=300, transparent=True)
