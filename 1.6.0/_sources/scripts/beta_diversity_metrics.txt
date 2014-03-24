.. index:: beta_diversity_metrics

*beta_diversity_metrics* -- List of available metrics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*Non-phylogenetic beta diversity metrics.* These are count based metrics which is based on the OTU table:

	* abund_jaccard: abundance weighted Jaccard distance
	* binary_dist_chisq: Binary chi-square distance
	* binary_dist_chord: Binary chord distance
	* binary_dist_euclidean: Binary euclidean distance
	* binary_dist_hamming: Binary Hamming distance (binary Manhattan distance)
	* binary_dist_jaccard: Binary Jaccard distance (binary Soergel distance)
	* binary_dist_lennon: Binary Lennon distance
	* binary_dist_ochiai: Binary Ochiai distance
	* binary_dist_pearson: Binary Pearson distance
	* binary_dist_sorensen_dice: Binary Sörensen-Dice distance (binary Bray-Curtis distance or binary Whittaker distance)
	* dist_bray_curtis: Bray-Curtis distance (normalized Manhattan distance)
	* dist_canberra: Canberra distance
	* dist_chisq: Chi-square distance
	* dist_chord: Chord distance
	* dist_euclidean: Euclidean distance
	* dist_gower: Gower distance
	* dist_hellinger: Hellinger distance
	* dist_kulczynski: Kulczynski distance
	* dist_manhattan: Manhattan distance
	* dist_morisita_horn: Morisita-Horn distance
	* dist_pearson: Pearson distance
	* dist_soergel: Soergel distance
	* dist_spearman_approx: Spearman rank distance
	* dist_specprof: Species profile distance

*Phylogenetic beta diversity metrics.* These metrics are based on UniFrac, which takes into account the evolutionary relationship between sequences:

	* dist_unifrac_G: The G metric calculates the fraction branch length in the sample i + sample j tree that is exclusive to sample i and it is asymmetric.
	* dist_unifrac_G_full_tree: The full_tree version calculates the fraction of branch length in the full tree that is exclusive to sample i and it is asymmetric.
	* dist_unweighted_unifrac: This is the standard unweighted UniFrac, which is used to assess 'who's there' without taking in account the relative abundance of identical sequences.
	* dist_unweighted_unifrac_full_tree: Typically, when computing the dissimilarity between two samples, unifrac considers only the parts of the phylogenetic tree contained by otus in either sample. The full_tree modification considers the entire supplied tree.
	* dist_weighted_normalized_unifrac: Weighted UniFrac with normalized values and is used to include abundance information. The normalization adjusts for varying root-to-tip distances.
	* dist_weighted_unifrac: Weighted UniFrac.
