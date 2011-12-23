#!/bin/bash
# interactive commands are commented out
#print_qiime_config.py

# Pre-processing
echo "Check mapping file"
rm -rf mapping_output ; check_id_map.py -m Fasting_Map.txt -o mapping_output -v

echo "Demultiplexing"
rm -rf split_library_output ; split_libraries.py -m Fasting_Map.txt -f Fasting_Example.fna -q Fasting_Example.qual -o split_library_output

# otus
echo "Pick OTUs through OTU table"
rm -rf otus ; pick_otus_through_otu_table.py -i split_library_output/seqs.fna -o otus -a

#per_library_stats.py -i otus/otu_table.biom

#OTU Heatmap
echo "OTU Heatmap"
make_otu_heatmap_html.py -i otus/otu_table.biom -o otus/OTU_Heatmap/

#OTU Network
echo "OTU Network"
make_otu_network.py -m Fasting_Map.txt -i otus/otu_table.biom -o otus/OTU_Network

#Make Taxa Summary Charts
echo "Summarize taxa"
rm -rf wf_taxa_summary ; summarize_taxa_through_plots.py -i otus/otu_table.biom -o wf_taxa_summary -m Fasting_Map.txt

echo "Alpha rarefaction"
#alpha_diversity.py -h
echo "alpha_diversity:metrics shannon,PD_whole_tree,chao1,observed_species" > alpha_params.txt

rm -rf wf_arare ; alpha_rarefaction.py -i otus/otu_table.biom -m Fasting_Map.txt -o wf_arare/ -p alpha_params.txt -t otus/rep_set.tre -a

echo "Beta diversity and plots"
rm -rf wf_bdiv_even146 ; beta_diversity_through_plots.py -i otus/otu_table.biom -m Fasting_Map.txt -o wf_bdiv_even146/ -t otus/rep_set.tre -e 146 -a

echo "Jackknifed beta diversity"
rm -rf wf_jack ; jackknifed_beta_diversity.py -i otus/otu_table.biom -t otus/rep_set.tre -m Fasting_Map.txt -o wf_jack -e 110 -a

echo "Make Bootstrapped Tree"
make_bootstrapped_tree.py -m wf_jack/unweighted_unifrac/upgma_cmp/master_tree.tre -s wf_jack/unweighted_unifrac/upgma_cmp/jackknife_support.txt -o wf_jack/unweighted_unifrac/upgma_cmp/jackknife_named_nodes.pdf

echo "Make Bi-Plots"
rm -rf 3d_biplot ; make_3d_plots.py -i wf_bdiv_even146/unweighted_unifrac_pc.txt -m Fasting_Map.txt -t wf_taxa_summary/otu_table_L3.txt --n_taxa_keep 5 -o 3d_biplot
