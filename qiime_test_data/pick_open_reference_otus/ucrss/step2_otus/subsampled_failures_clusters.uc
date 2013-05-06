# uclust --maxrejects 500 --input /Users/caporaso/temp/UclustExactMatchFilterbRMuLlFSNMZPN9cWMFBk.fasta --id 0.97 --w 12 --stepwords 20 --rev --usersort --maxaccepts 20 --stable_sort --uc /Users/caporaso/code/Qiime/qiime_test_data/pick_open_reference_otus/ucrss//step2_otus/subsampled_failures_clusters.uc
# version=1.2.22
# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
S	0	120	*	*	*	*	*	QiimeExactMatch.t1_56713	*
S	1	140	*	*	*	*	*	QiimeExactMatch.f2_10228	*
H	0	143	100.0	+	0	0	120M23D	QiimeExactMatch.t1_55836	QiimeExactMatch.t1_56713
H	0	141	100.0	+	0	0	120M21D	QiimeExactMatch.t2_56059	QiimeExactMatch.t1_56713
H	0	132	100.0	+	0	0	120M12D	QiimeExactMatch.t2_52739	QiimeExactMatch.t1_56713
C	0	4	100.0	*	*	*	*	QiimeExactMatch.t1_56713	*
C	1	1	*	*	*	*	*	QiimeExactMatch.f2_10228	*
