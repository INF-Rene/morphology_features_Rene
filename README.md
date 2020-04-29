# morphology_features_Rene
 Script for extracting many dendritic and axonal features from SWC or ASC files using the NeuroM package

Author: René Wilbers, renewilbers@gmail.com

# how to use:
 * Only tested with python 3.x . Needs the NeuroM package installed (https://neurom.readthedocs.io)
 * Open SWC_ASC_morph_feature_analysis.py and indicate the location of input & output folders
 * Make sure filenames are in the following format:
  [Species][xx]_[Layer]_[neuron_id], for example: HumanIN_L3_123456.SWC or MousePN_L23_123.ASC
 * Run Script
 * Features are exported to excel files.

# feature explanations:
* TAL/TDL/TDL_A/TDL_B --> total axonal length, total dendritic length, total dendritic length apicals, total dendritic length basals
* TAA/TDA etc --> surface area
* TAV/TDV etc --> volume
* ABP/DBP etc --> branchpoints
* leaf order --> leaf = final segment of a dendritic or axonal branch. Order depends on how often a segment has branched. So, a segment coming directly from the soma is 1st order, the segments branching from there are 2nd order, etc... Leaf order is calculated as the mean order of the leafs.
* reachvol--> volume of a polygonal 3D shape defined by the leaf endings  
* DSC --> dendritic stem count. number of primary dendrites stemming from the soma
* dendpathlength --> mean path length of soma to ends of the dendrites
* dendleaflength --> segment length of the leafs (terminal segments)
* dendbranchlength --> segment length of the branches (nonterminal segments)
* dendmaxpathlength --> length of the longest path from soma to the end of a dendrite
* dendmeanbranchdiam --> mean diameter of the branches (nonterminal segments)
* dendmeanleafdiam --> mean diameter of the leafs (terminal segments)
* axstart --> indicates how far away from the soma the axon starts (if it emerges from dendrite), 0 if it emerges from soma.
* axtodend --> TDA/TAA
* dend_um_per_BP --> TDL/DBP
* axon_um_per_BP --> TAL/ABP
* axonal efficiency --> Calculated as axreachvol^(1/3) / TAL. Says something about how efficiently the axon reaches around (so large reach with small amount of cable is considered efficient).
* dendric efficiency --> Same story, but  dendreachvol^(1/3) / TDL
* mean_total_stem_length --> TDL/DSC
* mean_total_stem_ branchpoints --> DBP/DSC 
