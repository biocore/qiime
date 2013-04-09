# uclust --maxrejects 500 --input /Users/caporaso/temp/UclustExactMatchFilterY1UiiWNdGoDSc1QAeVgH.fasta --id 0.97 --w 12 --stepwords 20 --usersort --maxaccepts 20 --libonly --stable_sort --lib /Users/caporaso/code/QiimeUtils/qiime_test_data/pick_iterative_otus_through_otu_table/ucrss_iter//0//new_refseqs.fna --uc /Users/caporaso/code/QiimeUtils/qiime_test_data/pick_iterative_otus_through_otu_table/ucrss_iter//1//step1_otus/prefiltered_seqs2_clusters.uc
# version=1.2.22
# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
N	*	252	*	*	*	*	*	QiimeExactMatch.PC.607_1000	*
L	0	248	*	*	*	*	*	r0 PC.635_779	*
H	0	248	100.0	+	0	0	248M	QiimeExactMatch.sample1_1	r0 PC.635_779
L	1	248	*	*	*	*	*	r1 PC.636_263	*
H	1	248	100.0	+	0	0	248M	QiimeExactMatch.PC.607_1	r1 PC.636_263
N	*	224	*	*	*	*	*	QiimeExactMatch.PC.607_2	*
D	0	2	*	*	*	*	100.0	r0 PC.635_779	*
D	1	2	*	*	*	*	100.0	r1 PC.636_263	*
