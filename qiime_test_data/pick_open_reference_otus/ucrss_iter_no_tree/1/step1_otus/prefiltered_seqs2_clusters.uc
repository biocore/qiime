# uclust --maxrejects 500 --input /Users/caporaso/temp/UclustExactMatchFilterHBSfp0erKamAea0afVrH.fasta --id 0.97 --w 12 --stepwords 20 --rev --usersort --maxaccepts 20 --libonly --stable_sort --lib /Users/caporaso/code/Qiime/qiime_test_data/pick_open_reference_otus/ucrss_iter_no_tree/0/new_refseqs.fna --uc /Users/caporaso/code/Qiime/qiime_test_data/pick_open_reference_otus/ucrss_iter_no_tree//1//step1_otus/prefiltered_seqs2_clusters.uc
# version=1.2.22
# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
L	3	1389	*	*	*	*	*	295053	*
H	3	133	100.0	+	0	0	519I133M737I	QiimeExactMatch.f5_1370	295053
N	*	133	*	*	*	*	*	QiimeExactMatch.f5_1300	*
L	8	152	*	*	*	*	*	New.0.CleanUp.ReferenceOTU4 f1_4866	*
H	8	152	100.0	+	0	0	152M	QiimeExactMatch.f5_1360	New.0.CleanUp.ReferenceOTU4 f1_4866
H	3	143	100.0	+	0	0	519I143M727I	QiimeExactMatch.f5_1470	295053
H	3	141	100.0	+	0	0	519I141M729I	QiimeExactMatch.f5_1510	295053
H	8	152	98.7	+	0	0	152M	QiimeExactMatch.f5_1400	New.0.CleanUp.ReferenceOTU4 f1_4866
D	3	4	*	*	*	*	100.0	295053	*
D	8	3	*	*	*	*	99.3	New.0.CleanUp.ReferenceOTU4 f1_4866	*
