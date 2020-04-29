# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 10:58:32 2019

@author: René Wilbers
"""
import numpy as np
import neurom as nm
from neurom import iter_sections, iter_segments, NeuriteType
import scipy
import matplotlib.pyplot as plt

def somaplot(ax, nrn):
#    from neurom import iter_sections
#    for sec in iter_sections(nrn):
#        if sec.type == NeuriteType.axon:
##            delete section from nrn
    nm.view.view.plot_neuron(ax, nrn, neurite_type=NeuriteType.basal_dendrite)
    ax.set_xlim(nrn.soma.center[0]-30, nrn.soma.center[0]+30)
    ax.set_ylim(nrn.soma.center[1]-30, nrn.soma.center[1]+30)
def nrnplot(ax, nrn):
    nm.view.view.plot_neuron(ax, nrn)
    ax.autoscale()

def find_remote_axons(neuron):
    remote_axons = list()
    for section in iter_sections(neuron):
        if section.is_root():
            continue
        if section.type == NeuriteType.axon and section.parent.type != NeuriteType.axon:
            remote_axons.append(section)
    return remote_axons

def get_sections(neuron, secname):
    secs = list()
    for section in iter_sections(neuron):
        if section.type == getattr(NeuriteType, secname):
            secs.append(section)
    return secs

def get_endleaf_data(secs):
    endleafs=list()
    leaflengths=list()
    leafdiams=list()
    leaforders=list()
    branchlengths=list()
    branchdiams=list()
    pathlengths=list()
    for i,sec in enumerate(secs):
        if sec.is_leaf():
            endleafs.append(sec)
            leaflengths.append(sec.length)
            pathlengths.append(0)
            leafdiams.append( np.sqrt((sec.volume/sec.length)/np.pi)*2 )
            leaforders.append(0)
            parsec=sec
            while not parsec.is_root():
                leaforders[-1]+=1
                pathlengths[-1]+=parsec.length
                parsec=parsec.parent
            pathlengths[-1]+=parsec.length
        else:
            branchlengths.append(sec.length)
            branchdiams.append( np.sqrt((sec.volume/sec.length)/np.pi)*2 )
    meanleaforder=np.array(leaforders).mean()
    meanpathlength=np.array(pathlengths).mean()
    meanleaflength=np.array(leaflengths).mean()
    meanbranchlength=np.array(branchlengths).mean()
    maxpathlength=np.array(pathlengths).max()
    meanleafdiam=np.array(leafdiams).mean()
    meanbranchdiam=np.array(branchdiams).mean()
    
    return meanleaforder, meanpathlength, meanleaflength, meanbranchlength, maxpathlength, meanleafdiam, meanbranchdiam
    
#    plt.figure()
#    plt.hist(np.array(pathlengths))
                
def get_dend_thicknesspath(dend):
    pathlengths=list()
    thicknesses=list()
    for sec in dend:
        thicknesses.append( np.sqrt((sec.volume/sec.length)/np.pi)*2 )
        pathlengths.append(sec.length/2)
        parsec=sec
        while not parsec.is_root():
            if not parsec is sec:
                pathlengths[-1]+=parsec.length
            parsec=parsec.parent
        pathlengths[-1]+=parsec.length
    return pathlengths, thicknesses
        
    
def get_seclist_measures_by_order(secs):
    orders=np.ones(len(secs)).astype('int')
    for i,sec in enumerate(secs):
        tmp=sec
        while tmp.parent:
            orders[i]+=1
            tmp=tmp.parent
    n=orders.max()
    TA=np.zeros(n)
    TL=np.zeros(n)
    BL=np.zeros(n)
    LL=np.zeros(n)
    TV=np.zeros(n)
    meanDIAM=np.zeros(n)
    secs=np.array(secs)
    for order in np.unique(orders):
        i=order-1
        ordersecs=secs[orders==order]
        lengths=list()
        leaflengths=list()
        branchlengths=list()
        diams=list()
        vols=list()
        for sec in ordersecs:
            if sec.is_leaf():
                leaflengths.append(sec.length)
            else:
                branchlengths.append(sec.length)
            lengths.append(sec.length)
            diams.append( np.sqrt((sec.volume/sec.length)/np.pi)*2 )
            vols.append(sec.volume)
            TA[i]+=sec.area
        TL[i]=np.array(lengths).sum()
        BL[i]=np.array(branchlengths).mean()
        LL[i]=np.array(leaflengths).mean()
        TV[i]=np.array(vols).sum()
        meanDIAM[i]=(np.array(diams)*np.array(lengths)).sum()/TL[i]
        
    return TL, meanDIAM, TV, BL, LL, TA
        

def get_stem_count(neuron, secname,parnames):
    partypes=[]
    for par in parnames:
        partypes.append(getattr(NeuriteType, par))
    SC=0
    for section in iter_sections(neuron):
        if section.type ==  getattr(NeuriteType, secname) and section.is_root():
            SC+=1
            continue
        if section.type ==  getattr(NeuriteType, secname) and section.parent.type in partypes:
            SC+=1
    return SC

def get_seclist_measures(secs):
    lengths=list()
    diams=list()
    vols=list()
    TA=0 #Total Area
    BP=0 # Branchpoints
    for sec in secs:
        lengths.append(sec.length)
        diams.append( np.sqrt((sec.volume/sec.length)/np.pi)*2 )
        vols.append(sec.volume)
        if sec.children:
            BP+=1
        TA+=sec.area
    lengths=np.array(lengths)
    diams=np.array(diams)
    vols=np.array(vols)
    TL=lengths.sum()
    BL=lengths.mean()
    TV=vols.sum()
    meanDIAM=(diams*lengths).sum()/TL
        
    return TL, meanDIAM, TV, BP, BL, TA # Total length, mean diam, total volume, branch points, mean branch length, total surface area

def get_sections_hullvolume(secs):
    cloud=np.empty((0,3))
    for sec in secs:
        cloud=np.vstack([cloud, sec.points[:,0:3]])
    hull=scipy.spatial.ConvexHull(cloud)
    volume=hull.volume
    return volume
#    plt.scatter(cloud[:,0], cloud[:,1])
#    plt.plot(cloud[hull.vertices,0], cloud[hull.vertices,1], 'r--', lw=2)