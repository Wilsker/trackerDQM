##########################################
#       clusterDigiOccupancies.py
##########################################
# Author: Joshuha Thomas-Wilsker
##########################################
# Short script that reads DQM ROOT files
# namely 3 files taken during a high
# luminosity run that have 3 different
# bunch patterns, loads cluster/digi
# occupancy TkHMaps and calculates the
# maximum and median values for the
# distributions in the different layers
# of the SiStrip tracker.
##########################################
from ROOT import TFile, TTree, gDirectory, gROOT, TH1, TF1, TProfile, TProfile2D
from array import array
from numpy import median, amax
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt


def get_median_number_clusters(TkHMap_clus):

    nbinsx =  TkHMap_clus.GetXaxis().GetNbins()
    nbinsy = TkHMap_clus.GetYaxis().GetNbins()
    numXXX_entries = []

    for lentry in range(nbinsx):
        for mentry in range(nbinsy):
            if TkHMap_clus.GetBinContent(lentry,mentry) != 0:
                numXXX_entries.append(TkHMap_clus.GetBinContent(lentry,mentry))
    median_num_clusters = median(numXXX_entries)
    return median_num_clusters

def main():
    #To run on lxplus
    #infile_8b4e = TFile('/afs/cern.ch/work/j/jthomasw/private/IHEP/trackerDQM/highlumiruns_rootfiles/DQM_V0001_R000302674__ZeroBias8b4e1__Run2017D-PromptReco-v1__DQMIO.root')
    #infile_IsolatedBunches = TFile('/afs/cern.ch/work/j/jthomasw/private/IHEP/trackerDQM/highlumiruns_rootfiles/DQM_V0001_R000302674__ZeroBiasIsolatedBunches1__Run2017D-PromptReco-v1__DQMIO.root')
    #infile_NominalTrains = TFile('/afs/cern.ch/work/j/jthomasw/private/IHEP/trackerDQM/highlumiruns_rootfiles/DQM_V0001_R000302674__ZeroBiasNominalTrains1__Run2017D-PromptReco-v1__DQMIO.root')

    #To run on MAC
    infile_8b4e = TFile('/Users/joshuhathomas-wilsker/Documents/work/lxplus_remote/work/private/IHEP/trackerDQM/cluster_digi_occupancy_study/DQM_V0001_R000302674__ZeroBias8b4e1__Run2017D-PromptReco-v1__DQMIO.root')
    infile_IsolatedBunches = TFile('/Users/joshuhathomas-wilsker/Documents/work/lxplus_remote/work/private/IHEP/trackerDQM/cluster_digi_occupancy_study/DQM_V0001_R000302674__ZeroBiasIsolatedBunches1__Run2017D-PromptReco-v1__DQMIO.root')
    infile_NominalTrains = TFile('/Users/joshuhathomas-wilsker/Documents/work/lxplus_remote/work/private/IHEP/trackerDQM/cluster_digi_occupancy_study/DQM_V0001_R000302674__ZeroBiasNominalTrains1__Run2017D-PromptReco-v1__DQMIO.root')

    files_list = [infile_8b4e,infile_NominalTrains,infile_IsolatedBunches]
    outputFile = open('clusterOccupancyStats.xml','w+')
    subdetectors = ['TIB','TOB','TID','TEC']
    fwd_bwd = ['MINUS','PLUS']
    TIB_layers = ['layer_1','layer_2','layer_3','layer_4']
    TIB_suffix = ['TIB_L1','TIB_L2','TIB_L3','TIB_L4']
    TOB_layers = ['layer_1','layer_2','layer_3','layer_4','layer_5','layer_6']
    TOB_suffix = ['TOB_L1','TOB_L2','TOB_L3','TOB_L4','TOB_L5','TOB_L6']
    TID_layers = ['wheel_1','wheel_2','wheel_3']
    TIDM_suffix = ['TIDM_D1','TIDM_D2','TIDM_D3']
    TIDP_suffix = ['TIDP_D1','TIDP_D2','TIDP_D3']
    TEC_layers = ['wheel_1','wheel_2','wheel_3','wheel_4','wheel_5','wheel_6','wheel_7','wheel_8','wheel_9']
    TECM_suffix = ['TECM_W1','TECM_W2','TECM_W3','TECM_W4','TECM_W5','TECM_W6','TECM_W7','TECM_W8','TECM_W9']
    TECP_suffix = ['TECP_W1','TECP_W2','TECP_W3','TECP_W4','TECP_W5','TECP_W6','TECP_W7','TECP_W8','TECP_W9']

    for ientry in range(len(subdetectors)):
        if subdetectors[ientry] == "TID":
            outputFile.write(subdetectors[ientry] + '\n')
            for jentry in range(len(fwd_bwd)):
                outputFile.write(fwd_bwd[jentry] + '\n')
                clusters_median_8b4e_all_layers = []
                clusters_median_NomTrains_all_layers = []
                clusters_median_IsoBunches_all_layers = []
                clusters_max_8b4e_all_layers = []
                clusters_max_NomTrains_all_layers = []
                clusters_max_IsoBunches_all_layers = []
                digis_median_8b4e_all_layers = []
                digis_median_NomTrains_all_layers = []
                digis_median_IsoBunches_all_layers = []
                digis_max_8b4e_all_layers = []
                digis_max_NomTrains_all_layers = []
                digis_max_IsoBunches_all_layers = []

                for kentry in range(len(TID_layers)):
                    outputFile.write(TID_layers[kentry] + '\n')
                    directoryPath = 'DQMData/Run 302674/SiStrip/Run summary/MechanicalView/'
                    directoryPath+=(subdetectors[ientry])+'/'
                    directoryPath+=(fwd_bwd[jentry])+'/'
                    directoryPath+=(TID_layers[kentry])+'/'
                    median_values = []
                    maximum_values = []

                    for pentry in range(len(files_list)):
                        tempFile = files_list[pentry]
                        tempFile.cd(directoryPath)

                        TIDsuff = ''
                        if fwd_bwd[jentry] == 'MINUS':
                            TIDsuff = TIDM_suffix[kentry]
                        elif fwd_bwd[jentry] == 'PLUS':
                            TIDsuff = TIDP_suffix[kentry]
                        cluster_histo_name = 'TkHMap_NumberOfCluster_'+TIDsuff
                        digi_histo_name = 'TkHMap_NumberOfDigi_'+TIDsuff
                        TkHMap_clusters = gDirectory.Get(cluster_histo_name)


                        nbinsx =  TkHMap_clusters.GetXaxis().GetNbins()
                        nbinsy = TkHMap_clusters.GetYaxis().GetNbins()
                        numXXX_entries = []

                        for lentry in range(nbinsx):
                            for mentry in range(nbinsy):
                                if TkHMap_clusters.GetBinContent(lentry,mentry) != 0:
                                    numXXX_entries.append(TkHMap_clusters.GetBinContent(lentry,mentry))
                        if pentry == 0:
                            fileslist = '8b4e'
                            #new_median = get_median_number_clusters(TkHMap_clusters)
                            clusters_median_8b4e_all_layers.append(median(numXXX_entries))
                            clusters_max_8b4e_all_layers.append(amax(numXXX_entries))
                        elif pentry == 1:
                            fileslist = 'NominalTrains'
                            clusters_median_NomTrains_all_layers.append(median(numXXX_entries))
                            clusters_max_NomTrains_all_layers.append(amax(numXXX_entries))
                        elif pentry == 2:
                            fileslist = 'IsolatedBunches'
                            clusters_median_IsoBunches_all_layers.append(median(numXXX_entries))
                            clusters_max_IsoBunches_all_layers.append(amax(numXXX_entries))
                        outputFile.write(fileslist + '\n')
                        outputFile.write('Clusters ' + fileslist + ': median = ' + str(median(numXXX_entries)) + ' maximum = ' + str(amax(numXXX_entries)) + '\n')
                        median_values.append(median(numXXX_entries))
                        maximum_values.append(amax(numXXX_entries))
                        del numXXX_entries[:]
                        TkHMap_digis = gDirectory.Get(digi_histo_name)
                        nbinsx = TkHMap_digis.GetXaxis().GetNbins()
                        nbinsy = TkHMap_digis.GetYaxis().GetNbins()

                        for lentry in range(nbinsx):
                            for mentry in range(nbinsy):
                                if TkHMap_digis.GetBinContent(lentry,mentry) != 0:
                                    numXXX_entries.append(TkHMap_digis.GetBinContent(lentry,mentry))
                        if pentry == 0:
                            digis_median_8b4e_all_layers.append(median(numXXX_entries))
                            digis_max_8b4e_all_layers.append(amax(numXXX_entries))
                        elif pentry == 1:
                            fileslist = 'NominalTrains'
                            digis_median_NomTrains_all_layers.append(median(numXXX_entries))
                            digis_max_NomTrains_all_layers.append(amax(numXXX_entries))
                        elif pentry == 2:
                            fileslist = 'IsolatedBunches'
                            digis_median_IsoBunches_all_layers.append(median(numXXX_entries))
                            digis_max_IsoBunches_all_layers.append(amax(numXXX_entries))
                        outputFile.write('Digis ' + fileslist + ': median = ' + str(median(numXXX_entries)) + ' maximum = ' + str(amax(numXXX_entries)) + '\n')

                index = np.arange(len(TID_layers))
                xTickMarks = [TID_layers[i] for i in range(len(TID_layers))]
                median_clusters_fig = plt.figure()
                median_clusters_ax = median_clusters_fig.add_subplot(111)
                median_clusters_bar_8b4b = median_clusters_ax.bar(index, clusters_median_8b4e_all_layers, 0.2, color='b')
                median_clusters_bar_nomtrains = median_clusters_ax.bar(index+0.2, clusters_median_NomTrains_all_layers, 0.2, color='g')
                median_clusters_bar_isobunches = median_clusters_ax.bar(index+0.4, clusters_median_IsoBunches_all_layers, 0.2, color='r')
                median_clusters_ax.set_xlabel('Det. Layer')
                median_clusters_ax.set_ylabel('Cluster occupancy')
                median_clusters_ax.set_xticks(index+0.3)
                median_clusters_ax.set_xticklabels(xTickMarks)
                median_clusters_ax.legend( (median_clusters_bar_8b4b[0] ,median_clusters_bar_nomtrains[0], median_clusters_bar_isobunches[0]), ('8b4e', 'Nom. Trains', 'Iso. Bunches'), loc='upper center', bbox_to_anchor=(1,1), prop={'size': 10})
                median_clusters_fig.subplots_adjust(right=0.9)
                median_clusters_title_string = 'Median cluster occupancy per detector layer: TID ' + fwd_bwd[jentry]
                median_clusters_outfile_name = 'median_clus_occupancy_TID_'+fwd_bwd[jentry]
                plt.title(median_clusters_title_string)
                plt.savefig(median_clusters_outfile_name)

                max_clusters_fig = plt.figure()
                max_clusters_ax = max_clusters_fig.add_subplot(111)
                max_clusters_bar_8b4b = max_clusters_ax.bar(index, clusters_max_8b4e_all_layers, 0.2, color='b')
                max_clusters_bar_nomtrains = max_clusters_ax.bar(index+0.2, clusters_max_NomTrains_all_layers, 0.2, color='g')
                max_clusters_bar_isobunches = max_clusters_ax.bar(index+0.4, clusters_max_IsoBunches_all_layers, 0.2, color='r')
                max_clusters_ax.set_xlabel('Det. Layer')
                max_clusters_ax.set_ylabel('Cluster occupancy')
                max_clusters_ax.set_xticks(index+0.3)
                max_clusters_ax.set_xticklabels(xTickMarks)
                max_clusters_ax.legend( (max_clusters_bar_8b4b[0] ,max_clusters_bar_nomtrains[0], max_clusters_bar_isobunches[0]), ('8b4e', 'Nom. Trains', 'Iso. Bunches'), loc='upper center', bbox_to_anchor=(1,1), prop={'size': 10})
                max_clusters_fig.subplots_adjust(right=0.9)
                max_clusters_title_string = 'Maximum cluster occupancy per detector layer: TID ' + fwd_bwd[jentry]
                max_clusters_outfile_name = 'max_clus_occupancy_TID_'+fwd_bwd[jentry]
                plt.title(max_clusters_title_string)
                plt.savefig(max_clusters_outfile_name)

                median_digis_fig = plt.figure()
                median_digis_ax = median_digis_fig.add_subplot(111)
                median_digis_bar_8b4b = median_digis_ax.bar(index, digis_median_8b4e_all_layers, 0.2, color='b')
                median_digis_bar_nomtrains = median_digis_ax.bar(index+0.2, digis_median_NomTrains_all_layers, 0.2, color='g')
                median_digis_bar_isobunches = median_digis_ax.bar(index+0.4, digis_median_IsoBunches_all_layers, 0.2, color='r')
                median_digis_ax.set_xlabel('Det. Layer')
                median_digis_ax.set_ylabel('Digi occupancy')
                median_digis_ax.set_xticks(index+0.3)
                median_digis_ax.set_xticklabels(xTickMarks)
                median_digis_ax.legend( (median_digis_bar_8b4b[0] ,median_digis_bar_nomtrains[0], median_digis_bar_isobunches[0]), ('8b4e', 'Nom. Trains', 'Iso. Bunches'), loc='upper center', bbox_to_anchor=(1,1), prop={'size': 10})
                median_digis_fig.subplots_adjust(right=0.9)
                median_digis_title_string = 'Median digi occupancy per detector layer: TID ' + fwd_bwd[jentry]
                median_digis_outfile_name = 'median_digi_occupancy_TID_'+fwd_bwd[jentry]
                plt.title(median_digis_title_string)
                plt.savefig(median_digis_outfile_name)

                max_digis_fig = plt.figure()
                max_digis_ax = max_digis_fig.add_subplot(111)
                max_digis_bar_8b4b = max_digis_ax.bar(index, digis_max_8b4e_all_layers, 0.2, color='b')
                max_digis_bar_nomtrains = max_digis_ax.bar(index+0.2, digis_max_NomTrains_all_layers, 0.2, color='g')
                max_digis_bar_isobunches = max_digis_ax.bar(index+0.4, digis_max_IsoBunches_all_layers, 0.2, color='r')
                max_digis_ax.set_xlabel('Det. Layer')
                max_digis_ax.set_ylabel('Digi occupancy')
                max_digis_ax.set_xticks(index+0.3)
                max_digis_ax.set_xticklabels(xTickMarks)
                max_digis_ax.legend( (max_digis_bar_8b4b[0] , max_digis_bar_nomtrains[0], max_digis_bar_isobunches[0]), ('8b4e', 'Nom. Trains', 'Iso. Bunches'), loc='upper center', bbox_to_anchor=(1,1), prop={'size': 10})
                max_digis_fig.subplots_adjust(right=0.9)
                max_digis_title_string = 'Maximum digi occupancy per detector layer: TID ' + fwd_bwd[jentry]
                max_digis_outfile_name = 'max_digi_occupancy_TID_'+fwd_bwd[jentry]
                plt.title(max_digis_title_string)
                plt.savefig(max_digis_outfile_name)

        elif subdetectors[ientry] == 'TEC':
            outputFile.write(subdetectors[ientry] + '\n')
            for jentry in range(len(fwd_bwd)):
                outputFile.write(fwd_bwd[jentry] + '\n')
                clusters_median_8b4e_all_layers = []
                clusters_median_NomTrains_all_layers = []
                clusters_median_IsoBunches_all_layers = []
                clusters_max_8b4e_all_layers = []
                clusters_max_NomTrains_all_layers = []
                clusters_max_IsoBunches_all_layers = []
                digis_median_8b4e_all_layers = []
                digis_median_NomTrains_all_layers = []
                digis_median_IsoBunches_all_layers = []
                digis_max_8b4e_all_layers = []
                digis_max_NomTrains_all_layers = []
                digis_max_IsoBunches_all_layers = []
                for kentry in range(len(TEC_layers)):
                    outputFile.write(TEC_layers[kentry] + '\n')
                    directoryPath = 'DQMData/Run 302674/SiStrip/Run summary/MechanicalView/'
                    directoryPath+=(subdetectors[ientry])+'/'
                    directoryPath+=(fwd_bwd[jentry])+'/'
                    directoryPath+=(TEC_layers[kentry])+'/'
                    for pentry in range(len(files_list)):
                        tempFile = files_list[pentry]
                        tempFile.cd(directoryPath)
                        TECsuff = ''
                        if fwd_bwd[jentry] == 'MINUS':
                            TECsuff = TECM_suffix[kentry]
                        elif fwd_bwd[jentry] == 'PLUS':
                            TECsuff = TECP_suffix[kentry]
                        cluster_histo_name = 'TkHMap_NumberOfCluster_'+TECsuff
                        digi_histo_name = 'TkHMap_NumberOfDigi_'+TECsuff
                        TkHMap_clusters = gDirectory.Get(cluster_histo_name)
                        nbinsx =  TkHMap_clusters.GetXaxis().GetNbins()
                        nbinsy = TkHMap_clusters.GetYaxis().GetNbins()
                        numXXX_entries = []
                        for lentry in range(nbinsx):
                            for mentry in range(nbinsy):
                                if TkHMap_clusters.GetBinContent(lentry,mentry) != 0:
                                    numXXX_entries.append(TkHMap_clusters.GetBinContent(lentry,mentry))
                        if pentry == 0:
                            fileslist = '8b4e'
                            clusters_median_8b4e_all_layers.append(median(numXXX_entries))
                            clusters_max_8b4e_all_layers.append(amax(numXXX_entries))
                        elif pentry == 1:
                            fileslist = 'NominalTrains'
                            clusters_median_NomTrains_all_layers.append(median(numXXX_entries))
                            clusters_max_NomTrains_all_layers.append(amax(numXXX_entries))
                        elif pentry == 2:
                            fileslist = 'IsolatedBunches'
                            clusters_median_IsoBunches_all_layers.append(median(numXXX_entries))
                            clusters_max_IsoBunches_all_layers.append(amax(numXXX_entries))

                        outputFile.write(fileslist + '\n')
                        outputFile.write('Clusters median = ' + str(median(numXXX_entries)) + ' maximum = ' + str(amax(numXXX_entries)) + '\n')
                        del numXXX_entries[:]
                        TkHMap_digis = gDirectory.Get(digi_histo_name)
                        nbinsx = TkHMap_digis.GetXaxis().GetNbins()
                        nbinsy = TkHMap_digis.GetYaxis().GetNbins()
                        for lentry in range(nbinsx):
                            for mentry in range(nbinsy):
                                if TkHMap_digis.GetBinContent(lentry,mentry) != 0:
                                    numXXX_entries.append(TkHMap_digis.GetBinContent(lentry,mentry))
                        if pentry == 0:
                            digis_median_8b4e_all_layers.append(median(numXXX_entries))
                            digis_max_8b4e_all_layers.append(amax(numXXX_entries))
                        elif pentry == 1:
                            fileslist = 'NominalTrains'
                            digis_median_NomTrains_all_layers.append(median(numXXX_entries))
                            digis_max_NomTrains_all_layers.append(amax(numXXX_entries))
                        elif pentry == 2:
                            fileslist = 'IsolatedBunches'
                            digis_median_IsoBunches_all_layers.append(median(numXXX_entries))
                            digis_max_IsoBunches_all_layers.append(amax(numXXX_entries))
                        outputFile.write('Digis median = ' + str(median(numXXX_entries)) + ' maximum = ' + str(amax(numXXX_entries)) + '\n')

                index = np.arange(len(TEC_layers))
                xTickMarks = [TEC_layers[i] for i in range(len(TEC_layers))]
                median_clusters_fig = plt.figure()
                median_clusters_ax = median_clusters_fig.add_subplot(111)
                median_clusters_bar_8b4b = median_clusters_ax.bar(index, clusters_median_8b4e_all_layers, 0.2, color='b')
                median_clusters_bar_nomtrains = median_clusters_ax.bar(index+0.2, clusters_median_NomTrains_all_layers, 0.2, color='g')
                median_clusters_bar_isobunches = median_clusters_ax.bar(index+0.4, clusters_median_IsoBunches_all_layers, 0.2, color='r')
                median_clusters_ax.set_xlabel('Det. Layer')
                median_clusters_ax.set_ylabel('Cluster occupancy')
                median_clusters_ax.set_xticks(index+0.3)
                median_clusters_ax.set_xticklabels(xTickMarks)
                median_clusters_ax.legend( (median_clusters_bar_8b4b[0] ,median_clusters_bar_nomtrains[0], median_clusters_bar_isobunches[0]), ('8b4e', 'Nom. Trains', 'Iso. Bunches'), loc='upper center', bbox_to_anchor=(1,1), prop={'size': 10})
                median_clusters_fig.subplots_adjust(right=0.9)
                median_clusters_title_string = 'Median cluster occupancy per detector layer: TEC ' + fwd_bwd[jentry]
                median_clusters_outfile_name = 'median_clus_occupancy_TEC_'+fwd_bwd[jentry]
                plt.title(median_clusters_title_string)
                plt.savefig(median_clusters_outfile_name)

                max_clusters_fig = plt.figure()
                max_clusters_ax = max_clusters_fig.add_subplot(111)
                max_clusters_bar_8b4b = max_clusters_ax.bar(index, clusters_max_8b4e_all_layers, 0.2, color='b')
                max_clusters_bar_nomtrains = max_clusters_ax.bar(index+0.2, clusters_max_NomTrains_all_layers, 0.2, color='g')
                max_clusters_bar_isobunches = max_clusters_ax.bar(index+0.4, clusters_max_IsoBunches_all_layers, 0.2, color='r')
                max_clusters_ax.set_xlabel('Det. Layer')
                max_clusters_ax.set_ylabel('Cluster occupancy')
                max_clusters_ax.set_xticks(index+0.3)
                max_clusters_ax.set_xticklabels(xTickMarks)
                max_clusters_ax.legend( (max_clusters_bar_8b4b[0] ,max_clusters_bar_nomtrains[0], max_clusters_bar_isobunches[0]), ('8b4e', 'Nom. Trains', 'Iso. Bunches'), loc='upper center', bbox_to_anchor=(1,1), prop={'size': 10})
                max_clusters_fig.subplots_adjust(right=0.9)
                max_clusters_title_string = 'Maximum cluster occupancy per detector layer: TEC ' + fwd_bwd[jentry]
                max_clusters_outfile_name = 'max_clus_occupancy_TEC_'+fwd_bwd[jentry]
                plt.title(max_clusters_title_string)
                plt.savefig(max_clusters_outfile_name)

                median_digis_fig = plt.figure()
                median_digis_ax = median_digis_fig.add_subplot(111)
                median_digis_bar_8b4b = median_digis_ax.bar(index, digis_median_8b4e_all_layers, 0.2, color='b')
                median_digis_bar_nomtrains = median_digis_ax.bar(index+0.2, digis_median_NomTrains_all_layers, 0.2, color='g')
                median_digis_bar_isobunches = median_digis_ax.bar(index+0.4, digis_median_IsoBunches_all_layers, 0.2, color='r')
                median_digis_ax.set_xlabel('Det. Layer')
                median_digis_ax.set_ylabel('Digi occupancy')
                median_digis_ax.set_xticks(index+0.3)
                median_digis_ax.set_xticklabels(xTickMarks)
                median_digis_ax.legend( (median_digis_bar_8b4b[0] ,median_digis_bar_nomtrains[0], median_digis_bar_isobunches[0]), ('8b4e', 'Nom. Trains', 'Iso. Bunches'), loc='upper center', bbox_to_anchor=(1,1), prop={'size': 10})
                median_digis_fig.subplots_adjust(right=0.9)
                median_digis_title_string = 'Median digi occupancy per detector layer: TEC ' + fwd_bwd[jentry]
                median_digis_outfile_name = 'median_digi_occupancy_TEC_'+fwd_bwd[jentry]
                plt.title(median_digis_title_string)
                plt.savefig(median_digis_outfile_name)

                max_digis_fig = plt.figure()
                max_digis_ax = max_digis_fig.add_subplot(111)
                max_digis_bar_8b4b = max_digis_ax.bar(index, digis_max_8b4e_all_layers, 0.2, color='b')
                max_digis_bar_nomtrains = max_digis_ax.bar(index+0.2, digis_max_NomTrains_all_layers, 0.2, color='g')
                max_digis_bar_isobunches = max_digis_ax.bar(index+0.4, digis_max_IsoBunches_all_layers, 0.2, color='r')
                max_digis_ax.set_xlabel('Det. Layer')
                max_digis_ax.set_ylabel('Digi occupancy')
                max_digis_ax.set_xticks(index+0.3)
                max_digis_ax.set_xticklabels(xTickMarks)
                max_digis_ax.legend( (max_digis_bar_8b4b[0] , max_digis_bar_nomtrains[0], max_digis_bar_isobunches[0]), ('8b4e', 'Nom. Trains', 'Iso. Bunches'), loc='upper center', bbox_to_anchor=(1,1), prop={'size': 10})
                max_digis_fig.subplots_adjust(right=0.9)
                max_digis_title_string = 'Maximum digi occupancy per detector layer: TEC ' + fwd_bwd[jentry]
                max_digis_outfile_name = 'max_digi_occupancy_TEC_'+fwd_bwd[jentry]
                plt.title(max_digis_title_string)
                plt.savefig(max_digis_outfile_name)

        elif subdetectors[ientry] == 'TIB':
            outputFile.write(subdetectors[ientry] + '\n')
            clusters_median_8b4e_all_layers = []
            clusters_median_NomTrains_all_layers = []
            clusters_median_IsoBunches_all_layers = []
            clusters_max_8b4e_all_layers = []
            clusters_max_NomTrains_all_layers = []
            clusters_max_IsoBunches_all_layers = []
            digis_median_8b4e_all_layers = []
            digis_median_NomTrains_all_layers = []
            digis_median_IsoBunches_all_layers = []
            digis_max_8b4e_all_layers = []
            digis_max_NomTrains_all_layers = []
            digis_max_IsoBunches_all_layers = []
            for jentry in range(len(TIB_layers)):
                outputFile.write(TIB_layers[jentry] + '\n')
                directoryPath = 'DQMData/Run 302674/SiStrip/Run summary/MechanicalView/'
                directoryPath+=(subdetectors[ientry])+'/'
                directoryPath+=(TIB_layers[jentry])+'/'
                for pentry in range(len(files_list)):
                    tempFile = files_list[pentry]
                    tempFile.cd(directoryPath)
                    cluster_histo_name = 'TkHMap_NumberOfCluster_'+TIB_suffix[jentry]
                    digi_histo_name = 'TkHMap_NumberOfDigi_'+TIB_suffix[jentry]
                    TkHMap_clusters = gDirectory.Get(cluster_histo_name)
                    nbinsx =  TkHMap_clusters.GetXaxis().GetNbins()
                    nbinsy = TkHMap_clusters.GetYaxis().GetNbins()
                    numXXX_entries = []
                    for lentry in range(nbinsx):
                        for mentry in range(nbinsy):
                            if TkHMap_clusters.GetBinContent(lentry,mentry) != 0:
                                numXXX_entries.append(TkHMap_clusters.GetBinContent(lentry,mentry))
                    if pentry == 0:
                        fileslist = '8b4e'
                        clusters_median_8b4e_all_layers.append(median(numXXX_entries))
                        clusters_max_8b4e_all_layers.append(amax(numXXX_entries))
                    elif pentry == 1:
                        fileslist = 'NominalTrains'
                        clusters_median_NomTrains_all_layers.append(median(numXXX_entries))
                        clusters_max_NomTrains_all_layers.append(amax(numXXX_entries))
                    elif pentry == 2:
                        fileslist = 'IsolatedBunches'
                        clusters_median_IsoBunches_all_layers.append(median(numXXX_entries))
                        clusters_max_IsoBunches_all_layers.append(amax(numXXX_entries))
                    outputFile.write(fileslist + '\n')
                    outputFile.write('Clusters median = ' + str(median(numXXX_entries)) + ' maximum = ' + str(amax(numXXX_entries)) + '\n')
                    del numXXX_entries[:]
                    TkHMap_digis = gDirectory.Get(digi_histo_name)
                    nbinsx = TkHMap_digis.GetXaxis().GetNbins()
                    nbinsy = TkHMap_digis.GetYaxis().GetNbins()
                    for lentry in range(nbinsx):
                        for mentry in range(nbinsy):
                            if TkHMap_digis.GetBinContent(lentry,mentry) != 0:
                                numXXX_entries.append(TkHMap_digis.GetBinContent(lentry,mentry))
                    if pentry == 0:
                        digis_median_8b4e_all_layers.append(median(numXXX_entries))
                        digis_max_8b4e_all_layers.append(amax(numXXX_entries))
                    elif pentry == 1:
                        fileslist = 'NominalTrains'
                        digis_median_NomTrains_all_layers.append(median(numXXX_entries))
                        digis_max_NomTrains_all_layers.append(amax(numXXX_entries))
                    elif pentry == 2:
                        fileslist = 'IsolatedBunches'
                        digis_median_IsoBunches_all_layers.append(median(numXXX_entries))
                        digis_max_IsoBunches_all_layers.append(amax(numXXX_entries))
                    outputFile.write('Digis median = ' + str(median(numXXX_entries)) + ' maximum = ' + str(amax(numXXX_entries)) + '\n')

            index = np.arange(len(TIB_layers))
            xTickMarks = [TIB_layers[i] for i in range(len(TIB_layers))]
            median_clusters_fig = plt.figure()
            median_clusters_ax = median_clusters_fig.add_subplot(111)
            median_clusters_bar_8b4b = median_clusters_ax.bar(index, clusters_median_8b4e_all_layers, 0.2, color='b')
            median_clusters_bar_nomtrains = median_clusters_ax.bar(index+0.2, clusters_median_NomTrains_all_layers, 0.2, color='g')
            median_clusters_bar_isobunches = median_clusters_ax.bar(index+0.4, clusters_median_IsoBunches_all_layers, 0.2, color='r')
            median_clusters_ax.set_xlabel('Det. Layer')
            median_clusters_ax.set_ylabel('Cluster occupancy')
            median_clusters_ax.set_xticks(index+0.3)
            median_clusters_ax.set_xticklabels(xTickMarks)
            median_clusters_ax.legend( (median_clusters_bar_8b4b[0] ,median_clusters_bar_nomtrains[0], median_clusters_bar_isobunches[0]), ('8b4e', 'Nom. Trains', 'Iso. Bunches'), loc='upper center', bbox_to_anchor=(1,1), prop={'size': 10})
            median_clusters_fig.subplots_adjust(right=0.9)
            median_clusters_title_string = 'Median cluster occupancy per detector layer: TIB '
            median_clusters_outfile_name = 'median_clus_occupancy_TIB'
            plt.title(median_clusters_title_string)
            plt.savefig(median_clusters_outfile_name)

            max_clusters_fig = plt.figure()
            max_clusters_ax = max_clusters_fig.add_subplot(111)
            max_clusters_bar_8b4b = max_clusters_ax.bar(index, clusters_max_8b4e_all_layers, 0.2, color='b')
            max_clusters_bar_nomtrains = max_clusters_ax.bar(index+0.2, clusters_max_NomTrains_all_layers, 0.2, color='g')
            max_clusters_bar_isobunches = max_clusters_ax.bar(index+0.4, clusters_max_IsoBunches_all_layers, 0.2, color='r')
            max_clusters_ax.set_xlabel('Det. Layer')
            max_clusters_ax.set_ylabel('Cluster occupancy')
            max_clusters_ax.set_xticks(index+0.3)
            max_clusters_ax.set_xticklabels(xTickMarks)
            max_clusters_ax.legend( (max_clusters_bar_8b4b[0] ,max_clusters_bar_nomtrains[0], max_clusters_bar_isobunches[0]), ('8b4e', 'Nom. Trains', 'Iso. Bunches'), loc='upper center', bbox_to_anchor=(1,1), prop={'size': 10})
            max_clusters_fig.subplots_adjust(right=0.9)
            max_clusters_title_string = 'Maximum cluster occupancy per detector layer: TIB '
            max_clusters_outfile_name = 'max_clus_occupancy_TIB'
            plt.title(max_clusters_title_string)
            plt.savefig(max_clusters_outfile_name)

            median_digis_fig = plt.figure()
            median_digis_ax = median_digis_fig.add_subplot(111)
            median_digis_bar_8b4b = median_digis_ax.bar(index, digis_median_8b4e_all_layers, 0.2, color='b')
            median_digis_bar_nomtrains = median_digis_ax.bar(index+0.2, digis_median_NomTrains_all_layers, 0.2, color='g')
            median_digis_bar_isobunches = median_digis_ax.bar(index+0.4, digis_median_IsoBunches_all_layers, 0.2, color='r')
            median_digis_ax.set_xlabel('Det. Layer')
            median_digis_ax.set_ylabel('Digi occupancy')
            median_digis_ax.set_xticks(index+0.3)
            median_digis_ax.set_xticklabels(xTickMarks)
            median_digis_ax.legend( (median_digis_bar_8b4b[0] ,median_digis_bar_nomtrains[0], median_digis_bar_isobunches[0]), ('8b4e', 'Nom. Trains', 'Iso. Bunches'), loc='upper center', bbox_to_anchor=(1,1), prop={'size': 10})
            median_digis_fig.subplots_adjust(right=0.9)
            median_digis_title_string = 'Median digi occupancy per detector layer: TIB '
            median_digis_outfile_name = 'median_digi_occupancy_TIB'
            plt.title(median_digis_title_string)
            plt.savefig(median_digis_outfile_name)

            max_digis_fig = plt.figure()
            max_digis_ax = max_digis_fig.add_subplot(111)
            max_digis_bar_8b4b = max_digis_ax.bar(index, digis_max_8b4e_all_layers, 0.2, color='b')
            max_digis_bar_nomtrains = max_digis_ax.bar(index+0.2, digis_max_NomTrains_all_layers, 0.2, color='g')
            max_digis_bar_isobunches = max_digis_ax.bar(index+0.4, digis_max_IsoBunches_all_layers, 0.2, color='r')
            max_digis_ax.set_xlabel('Det. Layer')
            max_digis_ax.set_ylabel('Digi occupancy')
            max_digis_ax.set_xticks(index+0.3)
            max_digis_ax.set_xticklabels(xTickMarks)
            max_digis_ax.legend( (max_digis_bar_8b4b[0] , max_digis_bar_nomtrains[0], max_digis_bar_isobunches[0]), ('8b4e', 'Nom. Trains', 'Iso. Bunches'), loc='upper center', bbox_to_anchor=(1,1), prop={'size': 10})
            max_digis_fig.subplots_adjust(right=0.9)
            max_digis_title_string = 'Maximum digi occupancy per detector layer: TIB '
            max_digis_outfile_name = 'max_digi_occupancy_TIB'
            plt.title(max_digis_title_string)
            plt.savefig(max_digis_outfile_name)

        elif subdetectors[ientry] == 'TOB':
            outputFile.write(subdetectors[ientry] + '\n')
            clusters_median_8b4e_all_layers = []
            clusters_median_NomTrains_all_layers = []
            clusters_median_IsoBunches_all_layers = []
            clusters_max_8b4e_all_layers = []
            clusters_max_NomTrains_all_layers = []
            clusters_max_IsoBunches_all_layers = []
            digis_median_8b4e_all_layers = []
            digis_median_NomTrains_all_layers = []
            digis_median_IsoBunches_all_layers = []
            digis_max_8b4e_all_layers = []
            digis_max_NomTrains_all_layers = []
            digis_max_IsoBunches_all_layers = []
            for jentry in range(len(TOB_layers)):
                outputFile.write(TOB_layers[jentry] + '\n')
                directoryPath = 'DQMData/Run 302674/SiStrip/Run summary/MechanicalView/'
                directoryPath+=(subdetectors[ientry])+'/'
                directoryPath+=(TOB_layers[jentry])+'/'
                for pentry in range(len(files_list)):
                    tempFile = files_list[pentry]
                    tempFile.cd(directoryPath)
                    cluster_histo_name = 'TkHMap_NumberOfCluster_'+TOB_suffix[jentry]
                    digi_histo_name = 'TkHMap_NumberOfDigi_'+TOB_suffix[jentry]
                    TkHMap_clusters = gDirectory.Get(cluster_histo_name)
                    nbinsx =  TkHMap_clusters.GetXaxis().GetNbins()
                    nbinsy = TkHMap_clusters.GetYaxis().GetNbins()
                    numXXX_entries = []
                    for lentry in range(nbinsx):
                        for mentry in range(nbinsy):
                            if TkHMap_clusters.GetBinContent(lentry,mentry) != 0:
                                numXXX_entries.append(TkHMap_clusters.GetBinContent(lentry,mentry))
                    if pentry == 0:
                        fileslist = '8b4e'
                        clusters_median_8b4e_all_layers.append(median(numXXX_entries))
                        clusters_max_8b4e_all_layers.append(amax(numXXX_entries))
                    elif pentry == 1:
                        fileslist = 'NominalTrains'
                        clusters_median_NomTrains_all_layers.append(median(numXXX_entries))
                        clusters_max_NomTrains_all_layers.append(amax(numXXX_entries))
                    elif pentry == 2:
                        fileslist = 'IsolatedBunches'
                        clusters_median_IsoBunches_all_layers.append(median(numXXX_entries))
                        clusters_max_IsoBunches_all_layers.append(amax(numXXX_entries))
                    outputFile.write(fileslist + '\n')
                    outputFile.write('Clusters median = ' + str(median(numXXX_entries)) + ' maximum = ' + str(amax(numXXX_entries)) + '\n')
                    del numXXX_entries[:]
                    TkHMap_digis = gDirectory.Get(digi_histo_name)
                    nbinsx = TkHMap_digis.GetXaxis().GetNbins()
                    nbinsy = TkHMap_digis.GetYaxis().GetNbins()
                    for lentry in range(nbinsx):
                        for mentry in range(nbinsy):
                            if TkHMap_digis.GetBinContent(lentry,mentry) != 0:
                                numXXX_entries.append(TkHMap_digis.GetBinContent(lentry,mentry))
                    if pentry == 0:
                        digis_median_8b4e_all_layers.append(median(numXXX_entries))
                        digis_max_8b4e_all_layers.append(amax(numXXX_entries))
                    elif pentry == 1:
                        fileslist = 'NominalTrains'
                        digis_median_NomTrains_all_layers.append(median(numXXX_entries))
                        digis_max_NomTrains_all_layers.append(amax(numXXX_entries))
                    elif pentry == 2:
                        fileslist = 'IsolatedBunches'
                        digis_median_IsoBunches_all_layers.append(median(numXXX_entries))
                        digis_max_IsoBunches_all_layers.append(amax(numXXX_entries))
                    outputFile.write('Digis median = ' + str(median(numXXX_entries)) + ' maximum = ' + str(amax(numXXX_entries)) + '\n')

            index = np.arange(len(TOB_layers))
            xTickMarks = [TOB_layers[i] for i in range(len(TOB_layers))]
            median_clusters_fig = plt.figure()
            median_clusters_ax = median_clusters_fig.add_subplot(111)
            median_clusters_bar_8b4b = median_clusters_ax.bar(index, clusters_median_8b4e_all_layers, 0.2, color='b')
            median_clusters_bar_nomtrains = median_clusters_ax.bar(index+0.2, clusters_median_NomTrains_all_layers, 0.2, color='g')
            median_clusters_bar_isobunches = median_clusters_ax.bar(index+0.4, clusters_median_IsoBunches_all_layers, 0.2, color='r')
            median_clusters_ax.set_xlabel('Det. Layer')
            median_clusters_ax.set_ylabel('Cluster occupancy')
            median_clusters_ax.set_xticks(index+0.3)
            median_clusters_ax.set_xticklabels(xTickMarks)
            median_clusters_ax.legend( (median_clusters_bar_8b4b[0] ,median_clusters_bar_nomtrains[0], median_clusters_bar_isobunches[0]), ('8b4e', 'Nom. Trains', 'Iso. Bunches'), loc='upper center', bbox_to_anchor=(1,1), prop={'size': 10})
            median_clusters_fig.subplots_adjust(right=0.9)
            median_clusters_title_string = 'Median cluster occupancy per detector layer: TOB '
            median_clusters_outfile_name = 'median_clus_occupancy_TOB'
            plt.title(median_clusters_title_string)
            plt.savefig(median_clusters_outfile_name)

            max_clusters_fig = plt.figure()
            max_clusters_ax = max_clusters_fig.add_subplot(111)
            max_clusters_bar_8b4b = max_clusters_ax.bar(index, clusters_max_8b4e_all_layers, 0.2, color='b')
            max_clusters_bar_nomtrains = max_clusters_ax.bar(index+0.2, clusters_max_NomTrains_all_layers, 0.2, color='g')
            max_clusters_bar_isobunches = max_clusters_ax.bar(index+0.4, clusters_max_IsoBunches_all_layers, 0.2, color='r')
            max_clusters_ax.set_xlabel('Det. Layer')
            max_clusters_ax.set_ylabel('Cluster occupancy')
            max_clusters_ax.set_xticks(index+0.3)
            max_clusters_ax.set_xticklabels(xTickMarks)
            max_clusters_ax.legend( (max_clusters_bar_8b4b[0] ,max_clusters_bar_nomtrains[0], max_clusters_bar_isobunches[0]), ('8b4e', 'Nom. Trains', 'Iso. Bunches'), loc='upper center', bbox_to_anchor=(1,1), prop={'size': 10})
            max_clusters_fig.subplots_adjust(right=0.9)
            max_clusters_title_string = 'Maximum cluster occupancy per detector layer: TOB '
            max_clusters_outfile_name = 'max_clus_occupancy_TOB'
            plt.title(max_clusters_title_string)
            plt.savefig(max_clusters_outfile_name)

            median_digis_fig = plt.figure()
            median_digis_ax = median_digis_fig.add_subplot(111)
            median_digis_bar_8b4b = median_digis_ax.bar(index, digis_median_8b4e_all_layers, 0.2, color='b')
            median_digis_bar_nomtrains = median_digis_ax.bar(index+0.2, digis_median_NomTrains_all_layers, 0.2, color='g')
            median_digis_bar_isobunches = median_digis_ax.bar(index+0.4, digis_median_IsoBunches_all_layers, 0.2, color='r')
            median_digis_ax.set_xlabel('Det. Layer')
            median_digis_ax.set_ylabel('Digi occupancy')
            median_digis_ax.set_xticks(index+0.3)
            median_digis_ax.set_xticklabels(xTickMarks)
            median_digis_ax.legend( (median_digis_bar_8b4b[0] ,median_digis_bar_nomtrains[0], median_digis_bar_isobunches[0]), ('8b4e', 'Nom. Trains', 'Iso. Bunches'), loc='upper center', bbox_to_anchor=(1,1), prop={'size': 10})
            median_digis_fig.subplots_adjust(right=0.9)
            median_digis_title_string = 'Median digi occupancy per detector layer: TOB '
            median_digis_outfile_name = 'median_digi_occupancy_TOB'
            plt.title(median_digis_title_string)
            plt.savefig(median_digis_outfile_name)

            max_digis_fig = plt.figure()
            max_digis_ax = max_digis_fig.add_subplot(111)
            max_digis_bar_8b4b = max_digis_ax.bar(index, digis_max_8b4e_all_layers, 0.2, color='b')
            max_digis_bar_nomtrains = max_digis_ax.bar(index+0.2, digis_max_NomTrains_all_layers, 0.2, color='g')
            max_digis_bar_isobunches = max_digis_ax.bar(index+0.4, digis_max_IsoBunches_all_layers, 0.2, color='r')
            max_digis_ax.set_xlabel('Det. Layer')
            max_digis_ax.set_ylabel('Digi occupancy')
            max_digis_ax.set_xticks(index+0.3)
            max_digis_ax.set_xticklabels(xTickMarks)
            max_digis_ax.legend( (max_digis_bar_8b4b[0] , max_digis_bar_nomtrains[0], max_digis_bar_isobunches[0]), ('8b4e', 'Nom. Trains', 'Iso. Bunches'), loc='upper center', bbox_to_anchor=(1,1), prop={'size': 10})
            max_digis_fig.subplots_adjust(right=0.9)
            max_digis_title_string = 'Maximum digi occupancy per detector layer: TOB '
            max_digis_outfile_name = 'max_digi_occupancy_TOB'
            plt.title(max_digis_title_string)
            plt.savefig(max_digis_outfile_name)

main()
