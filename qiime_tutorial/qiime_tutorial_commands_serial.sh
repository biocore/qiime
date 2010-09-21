#!/bin/bash

# Pre-processing
echo "Check Mapping"
rm -rf mapping_output ; check_id_map.py -m Fasting_Map.txt -o mapping_output

echo "De-multiplexing"
rm -rf split_library_output ; split_libraries.py -m Fasting_Map.txt -f Fasting_Example.fna -q Fasting_Example.qual -o split_library_output

# Data analysis
echo "Pick OTUs through OTU table"
rm -rf wf_da ; pick_otus_through_otu_table.py -i split_library_output/seqs.fna -p custom_parameters.txt -o wf_da

#OTU table visualization
echo "Even-sampling rarefaction"
single_rarefaction.py -i wf_da/uclust_picked_otus/rep_set/rdp_assigned_taxonomy/otu_table/seqs_otu_table.txt -d 146 -o wf_da/uclust_picked_otus/rep_set/rdp_assigned_taxonomy/otu_table/rarified_otu_table.txt

#OTU Heatmap
echo "OTU Heatmap"
make_otu_heatmap_html.py -i wf_da/uclust_picked_otus/rep_set/rdp_assigned_taxonomy/otu_table/seqs_otu_table.txt -o wf_da/uclust_picked_otus/rep_set/rdp_assigned_taxonomy/otu_table/OTU_Heatmap

#OTU Network
echo "OTU Network"
make_otu_network.py -m Fasting_Map.txt -i wf_da/uclust_picked_otus/rep_set/rdp_assigned_taxonomy/otu_table/seqs_otu_table.txt -o wf_da/uclust_picked_otus/rep_set/rdp_assigned_taxonomy/otu_table/OTU_Network

#Make Pie Charts
echo "Summarize taxa"
summarize_taxa.py -i wf_da/uclust_picked_otus/rep_set/rdp_assigned_taxonomy/otu_table/seqs_otu_table.txt -o wf_da/uclust_picked_otus/rep_set/rdp_assigned_taxonomy/otu_table/otu_table_Level3.txt -L 3 -r 0

echo "Make Pie Charts"
make_pie_charts.py -i wf_da/uclust_picked_otus/rep_set/rdp_assigned_taxonomy/otu_table/otu_table_Level3.txt -l Phylum -o wf_da/uclust_picked_otus/rep_set/rdp_assigned_taxonomy/otu_table/Pie_Charts -k white -s

echo "Alpha rarefaction"
rm -rf wf_arare ; alpha_rarefaction.py -i wf_da/uclust_picked_otus/rep_set/rdp_assigned_taxonomy/otu_table/seqs_otu_table.txt -m Fasting_Map.txt -o wf_arare/ -p custom_parameters.txt -t wf_da/uclust_picked_otus/rep_set/pynast_aligned_seqs/fasttree_phylogeny/seqs_rep_set.tre

echo "Beta diversity/3d plots"
rm -rf wf_bdiv ; beta_diversity_through_3d_plots.py -i wf_da/uclust_picked_otus/rep_set/rdp_assigned_taxonomy/otu_table/seqs_otu_table.txt -m Fasting_Map.txt -o wf_bdiv/ -p custom_parameters.txt -t wf_da/uclust_picked_otus/rep_set/pynast_aligned_seqs/fasttree_phylogeny/seqs_rep_set.tre

echo "Make 2D Plots - Unweighted Unifrac"
make_2d_plots.py -i wf_bdiv/unweighted_unifrac_pc.txt -m Fasting_Map.txt -o wf_bdiv/unweighted_unifrac_2d -k white -p wf_bdiv/prefs.txt

echo "Make Distance Histograms - Unweighted Unifrac"
make_distance_histograms.py -d wf_bdiv/unweighted_unifrac_seqs_otu_table.txt -m Fasting_Map.txt -o wf_bdiv/Distance_Histograms -p wf_bdiv/prefs.txt

echo "Jackknifed beta diversity"
rm -rf wf_jack/; jackknifed_beta_diversity.py -i wf_da/uclust_picked_otus/rep_set/rdp_assigned_taxonomy/otu_table/seqs_otu_table.txt -o wf_jack -p custom_parameters.txt -e 110 -t wf_da/uclust_picked_otus/rep_set/pynast_aligned_seqs/fasttree_phylogeny/seqs_rep_set.tre -m Fasting_Map.txt

echo "Make Bootstrapped Tree"
make_bootstrapped_tree.py -m wf_jack/unweighted_unifrac/upgma_cmp/master_tree.tre -s wf_jack/unweighted_unifrac/upgma_cmp/jackknife_support.txt -o wf_jack/unweighted_unifrac/upgma_cmp/jackknife_named_nodes.pdf



#rm -rf mapping_output
#rm -rf split_library_output
#rm -rf wf_da
#rm -rf wf_arare
#rm -rf wf_bdiv
#rm -rf wf_jack
#rm -rf jobs