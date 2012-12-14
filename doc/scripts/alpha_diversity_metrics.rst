.. index:: alpha_diversity_metrics

*alpha_diversity_metrics* -- List of available metrics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*Non-phylogeny based metrics:*

	* berger_parker_d
	* brillouin_d
	* chao1
	* chao1_confidence
	* dominance
	* doubles (# otus with exactly two individuals in sample)
	* equitability
	* fisher_alpha
	* gini index
	* goods coverage
	* heip_e (note, using heip_e at low (<5) individuals may cause errors
	* kempton_taylor_q
	* margalef
	* mcintosh_d
	* mcintosh_e
	* menhinick
	* michaelis_menten_fit
	* observed_species
	* osd (observed # otus, singleton OTUs, doubleton OTUs)
	* robbins
	* shannon (base 2 is used in the logarithms)
	* simpson (1 - Dominance)
	* simpson_reciprocal (1 / Dominance)
	* simpson_e
	* singles (# OTUs with exactly one individual present in sample)
	* strong
	
*Phylogeny based metrics:*
	
	* PD_whole_tree
