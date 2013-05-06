# uclust --maxrejects 500 --input /Users/caporaso/temp/UclustExactMatchFilterTIohweRWyzP5hhXb2O72.fasta --id 0.97 --w 12 --stepwords 20 --rev --usersort --maxaccepts 20 --stable_sort --uc /Users/caporaso/code/Qiime/qiime_test_data/pick_open_reference_otus/ucrss_iter_no_tax//0//step2_otus/subsampled_failures_clusters.uc
# version=1.2.22
# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
S	0	152	*	*	*	*	*	QiimeExactMatch.f1_19304	*
S	1	137	*	*	*	*	*	QiimeExactMatch.p2_35810	*
S	2	143	*	*	*	*	*	QiimeExactMatch.t1_53967	*
H	2	126	100.0	+	0	0	126M17I	QiimeExactMatch.p2_33176	QiimeExactMatch.t1_53967
H	2	120	100.0	+	0	0	120M23I	QiimeExactMatch.t1_56790	QiimeExactMatch.t1_53967
C	0	1	*	*	*	*	*	QiimeExactMatch.f1_19304	*
C	1	1	*	*	*	*	*	QiimeExactMatch.p2_35810	*
C	2	3	100.0	*	*	*	*	QiimeExactMatch.t1_53967	*
