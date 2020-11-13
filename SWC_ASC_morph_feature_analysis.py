# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 16:03:39 2019

@author: René Wilbers
"""

#user settings
output_name = 'SWC_test' # will be used for excel file names
morphdir = r'C:\YourInputDirectory' # directory of morphologies. optional .marker files can be added in this folder to indicate truncation points
outputdir= r'C:\YourOutputDirectory'  #output directory

#morphology properties

import os
import neurom as nm
from neurom import viewer
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from morph_analysis_tools import *

allfiles=os.listdir(morphdir)
files=[]
markerfiles=[]
for file in allfiles:
    if file.endswith(".swc"):
        files.append(file)
    if file.endswith(".marker"):
        markerfiles.append(file)


df=pd.DataFrame(columns=['id','species','layer',
                         'TAL', 'axdiam','TAV', 'axstart','ABP','axreachvol','TAA',
                         'TDL', 'denddiam', 'TDV', 'DBP', 'dendreachvol', 'TDA', 'DSC',
                         'dendleaforder', 'dendpathlength', 'dendleaflength', 'dendbranchlength', 'dendmaxpathlength',
                         'dendmeanleafdiam','dendmeanbranchdiam',
                         'TDL_A', 'denddiam_A', 'TDV_A', 'DBP_A', 'dendreachvol_A', 'TDA_A', 'DSC_A',
                         'dendleaforder_A', 'dendpathlength_A', 'dendleaflength_A', 'dendbranchlength_A', 'dendmaxpathlength_A',
                         'dendmeanleafdiam_A','dendmeanbranchdiam_A',
                         'TDL_B', 'denddiam_B', 'TDV_B', 'DBP_B', 'dendreachvol_B', 'TDA_B', 'DSC_B',
                         'dendleaforder_B', 'dendpathlength_B', 'dendleaflength_B', 'dendbranchlength_B', 'dendmaxpathlength_B',
                         'dendmeanleafdiam_B','dendmeanbranchdiam_B',
                         'axleaforder', 'axpathlength', 'axleaflength', 'axbranchlength', 'axmaxpathlength',
                         'full_dendrites', 'BPs_per_full_dendrite'])
orderdf=pd.DataFrame(columns=['id','species','layer','order',
                         'TDL', 'denddiam', 'TDV', 'DBL', 'DLL','TDA'])
thicknesspath=pd.DataFrame(columns=['id','species','layer',
                                    'pathlengths', 'thicknesses'])

for index,file in enumerate(files):
    print('Analyzing file #'+str(index+1)+' out of '+str(len(files)))
#    if index>0:
#        fig.add_subplot(7,10,index+1, sharex=ax, sharey=ax)
#        plt.axis('off')
    fnsplit=file[:-4].split('_')
    nrnid =int(fnsplit[2])
    species=fnsplit[0][:-2]
    layer=fnsplit[1]
    try:
        nrn = nm.load_neuron(os.path.join(morphdir, file))
        markerfn =[a for a in markerfiles if str(nrnid) in a] 
        if markerfn:
            markerdata=pd.read_csv(os.path.join(morphdir, markerfn[0]))
            if 'name' in markerdata.columns:
               markerdata=markerdata.drop(markerdata[markerdata['name']==30].index)
        else:
            print('no markerfile found for neuron ' + str(nrnid))
                
        #get neuron stats
        axon=get_sections(nrn, 'axon')
        dendBasal=get_sections(nrn, 'basal_dendrite')
        dendApical=get_sections(nrn, 'apical_dendrite')
        dend = dendBasal + dendApical
    
    
        if axon:
            TAL, axdiam, TAV, ABP, ABL, TAA  = get_seclist_measures(axon)
            ASC = get_stem_count(nrn, 'axon', ['soma', 'basal_dendrite'])
            axreachvol=get_sections_hullvolume(axon)
            axleaforder, axpathlength, axleaflength, axbranchlength, axmaxpathlength, meanleafdiam, meanbranchdiam = get_endleaf_data(axon, markerdata)
        else:
            TAL= axdiam= TAV= ABP= ABL= TAA= ASC =axreachvol =axleaforder= axpathlength= axleaflength= axbranchlength= axmaxpathlength= meanleafdiam= meanbranchdiamaxreachvol= float('nan')
        
        DSC_A = get_stem_count(nrn, 'apical_dendrite', ['soma'])
        DSC_B = get_stem_count(nrn, 'basal_dendrite', ['soma'])
        
        if dend:
            TDL, denddiam, TDV, DBP, DBL, TDA = get_seclist_measures(dend) 
            DSC = DSC_A+DSC_B
            dendreachvol=get_sections_hullvolume(dend)
            dendleaforder, dendpathlength, dendleaflength, dendbranchlength, dendmaxpathlength, dendmeanleafdiam, dendmeanbranchdiam = get_endleaf_data(dend,markerdata)
            full_dendrites, BPs_per_full_dendrite = get_full_dendrite_measures(dend, markerdata)
        else:
            TDL= denddiam= TDV= DBP= DBL= TDA= DSC= dendreachvol= dendleaforder= dendpathlength= dendleaflength= dendbranchlength= dendmaxpathlength= dendmeanleafdiam= dendmeanbranchdiam= full_dendrites= BPs_per_full_dendrite = float('nan')
            
        if dendApical:
            TDL_A, denddiam_A, TDV_A, DBP_A, DBL_A, TDA_A = get_seclist_measures(dendApical)
            dendreachvol_A=get_sections_hullvolume(dendApical)
            dendleaforder_A, dendpathlength_A, dendleaflength_A, dendbranchlength_A, dendmaxpathlength_A, dendmeanleafdiam_A, dendmeanbranchdiam_A = get_endleaf_data(dendApical, markerdata)
        else:
            TDL_A= denddiam_A= TDV_A= DBP_A= DBL_A= TDA_A=DSC_A= dendreachvol_A= dendleaforder_A= dendpathlength_A= dendleaflength_A= dendbranchlength_A= dendmaxpathlength_A= dendmeanleafdiam_A= dendmeanbranchdiam_A = float('nan')
        
        if dendBasal:
            TDL_B, denddiam_B, TDV_B, DBP_B, DBL_B, TDA_B = get_seclist_measures(dendBasal) 
            dendreachvol_B=get_sections_hullvolume(dendBasal)
            dendleaforder_B, dendpathlength_B, dendleaflength_B, dendbranchlength_B, dendmaxpathlength_B, dendmeanleafdiam_B, dendmeanbranchdiam_B = get_endleaf_data(dendBasal, markerdata)
        else:
            TDL_B= denddiam_B= TDV_B= DBP_B= DBL_B= TDA_B=DSC_B= dendreachvol_B= dendleaforder_B= dendpathlength_B= dendleaflength_B= dendbranchlength_B= dendmaxpathlength_B= dendmeanleafdiam_B= dendmeanbranchdiam_B = float('nan')
            
        
        #axon start distance
        sec=find_remote_axons(nrn)
        if len(sec)==1:
            axstart=nm.fst.sectionfunc.section_path_length(sec[0])-sec[0].length
        if not sec:
            axstart=0
        if len(sec)>1:
            print('Multiple axons found in', file)
            axstart=nm.fst.sectionfunc.section_path_length(sec[0])-sec[0].length
            
            
        df=df.append({'id': nrnid, 'species': species,'layer':layer,
                      'TAL':TAL,'axdiam':axdiam,'TAV':TAV,'axstart':axstart,'ABP':ABP,'ABL':ABL,'axreachvol':axreachvol,'TAA':TAA,
                      'TDL':TDL, 'denddiam':denddiam, 'TDV':TDV, 'DBP':DBP, 'DBL':DBL, 'dendreachvol':dendreachvol, 'TDA':TDA,
                      'DSC':DSC, 
                      'dendleaforder':dendleaforder, 'dendpathlength':dendpathlength, 'dendleaflength':dendleaflength,
                      'dendbranchlength':dendbranchlength, 'dendmaxpathlength':dendmaxpathlength,
                      'dendmeanleafdiam':dendmeanleafdiam,'dendmeanbranchdiam':dendmeanbranchdiam,
                      'TDL_A':TDL_A, 'denddiam_A':denddiam_A, 'TDV_A':TDV_A, 'DBP_A':DBP_A, 'DBL_A':DBL_A, 'dendreachvol_A':dendreachvol_A, 'TDA_A':TDA_A,
                      'DSC_A':DSC_A, 
                      'dendleaforder_A':dendleaforder_A, 'dendpathlength_A':dendpathlength_A, 'dendleaflength_A':dendleaflength_A,
                      'dendbranchlength_A':dendbranchlength_A, 'dendmaxpathlength_A':dendmaxpathlength_A,
                      'dendmeanleafdiam_A':dendmeanleafdiam_A,'dendmeanbranchdiam_A':dendmeanbranchdiam_A,
                      'TDL_B':TDL_B, 'denddiam_B':denddiam_B, 'TDV_B':TDV_B, 'DBP_B':DBP_B, 'DBL_B':DBL_B, 'dendreachvol_B':dendreachvol_B, 'TDA_B':TDA_B,
                      'DSC_B':DSC_B, 
                      'dendleaforder_B':dendleaforder_B, 'dendpathlength_B':dendpathlength_B, 'dendleaflength_B':dendleaflength_B,
                      'dendbranchlength_B':dendbranchlength_B, 'dendmaxpathlength_B':dendmaxpathlength_B,
                      'dendmeanleafdiam_B':dendmeanleafdiam_B,'dendmeanbranchdiam_B':dendmeanbranchdiam_B,
                      'axleaforder':axleaforder, 'axpathlength':axpathlength, 'axleaflength':axleaflength,
                      'axbranchlength':axbranchlength, 'axmaxpathlength':axmaxpathlength,
                      'full_dendrites':full_dendrites, 'BPs_per_full_dendrite':BPs_per_full_dendrite},
                        ignore_index=True)
        # mean dendritic parameters per dendritic branch order
        TL, meanDIAM, TV, BL, LL, TA = get_seclist_measures_by_order(dend)
        
        for i in range(len(TL)):
            orderdf=orderdf.append({'id': nrnid, 'species': species,'layer':layer,'order':i+1,
                          'TDL':TL[i], 'denddiam':meanDIAM[i], 'TDV':TV[i], 'DBL':BL[i],'DLL':LL[i], 'TDA':TA[i],},
                            ignore_index=True)
        #dendritic pathlength vs thickness:
        pathlengths, thicknesses = get_dend_thicknesspath(dend)
        thicknesspath=thicknesspath.append({'id': nrnid, 'species': species,'layer':layer,
                          'pathlengths':pathlengths, 'thicknesses':thicknesses},
                            ignore_index=True)
    except:
        break
        print('Could not load', file)

df['axtodend']=df.TDA/df.TAA
#removed because zero BP can result in division by 0:
#df['dend_um_per_BP']=df.TDL/df.DBP
#df['axon_um_per_BP']=df.TAL/df.ABP
df['axonal_efficiency']=np.power(df.axreachvol, 1/3)/df.TAL
df['dendritic_efficiency']=np.power(df.dendreachvol, 1/3)/df.TDL
df['mean_total_stem_length']=df.TDL/df.DSC
df['mean_total_stem_branchpoints']=df.DBP/df.DSC


df.to_excel(os.path.join(outputdir, output_name+'.xlsx'))
orderdf.to_excel(os.path.join(outputdir, output_name+'_ByBranchOrder.xlsx'))
thicknesspath.to_excel(os.path.join(outputdir, output_name+'_DendriticThicknessPath.xlsx'))

# uncomment below code for saving histograms of segment thicknesses
#for cellid, diams in zip(thicknesspath['id'], thicknesspath['thicknesses']):
#    plt.figure()
#    plt.hist(diams)
#    plt.xlabel('Diameter(um)')
#    plt.ylabel('Segment Count')
#    plt.title('cell ' + str(cellid))
#    plt.savefig(os.path.join(outputdir, 'Diameterhist_cellid_'+str(cellid)), dpi=300 )
