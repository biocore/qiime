# uclust --maxrejects 500 --input /Users/caporaso/temp/UclustExactMatchFilterybGdwGFopNg7eqqOczEE.fasta --id 0.97 --w 12 --stepwords 20 --rev --usersort --maxaccepts 20 --stable_sort --uc /Users/caporaso/code/qiime/qiime_test_data/pick_open_reference_otus/ucrss_iter_no_tree//0//step2_otus/subsampled_failures_clusters.uc
# version=1.2.22
# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
S	0	126	*	*	*	*	*	QiimeExactMatch.p2_33176	*
S	1	137	*	*	*	*	*	QiimeExactMatch.p2_35810	*
S	2	152	*	*	*	*	*	QiimeExactMatch.f3_141	*
H	0	143	100.0	+	0	0	126M17D	QiimeExactMatch.t1_55836	QiimeExactMatch.p2_33176
H	0	128	100.0	+	0	0	126M2D	QiimeExactMatch.t2_53194	QiimeExactMatch.p2_33176
H	0	121	100.0	+	0	0	121M5I	QiimeExactMatch.t1_56191	QiimeExactMatch.p2_33176
H	0	120	100.0	+	0	0	120M6I	QiimeExactMatch.t1_57363	QiimeExactMatch.p2_33176
C	0	5	100.0	*	*	*	*	QiimeExactMatch.p2_33176	*
C	1	1	*	*	*	*	*	QiimeExactMatch.p2_35810	*
C	2	1	*	*	*	*	*	QiimeExactMatch.f3_141	*
