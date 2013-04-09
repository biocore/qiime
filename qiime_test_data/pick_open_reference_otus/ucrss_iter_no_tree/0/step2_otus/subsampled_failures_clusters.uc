# uclust --maxrejects 500 --input /Users/caporaso/temp/UclustExactMatchFilter49DBnlYkTUM6epLPvU9q.fasta --id 0.97 --w 12 --stepwords 20 --rev --usersort --maxaccepts 20 --stable_sort --uc /Users/caporaso/code/Qiime/qiime_test_data/pick_open_reference_otus/ucrss_iter_no_tree//0//step2_otus/subsampled_failures_clusters.uc
# version=1.2.22
# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
S	0	128	*	*	*	*	*	QiimeExactMatch.t2_53194	*
H	0	116	100.0	+	0	0	116M12I	QiimeExactMatch.t1_54144	QiimeExactMatch.t2_53194
C	0	2	100.0	*	*	*	*	QiimeExactMatch.t2_53194	*
