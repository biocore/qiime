# uclust --maxrejects 500 --input /Users/caporaso/temp/UclustExactMatchFiltertvLxSgdNQgs2P3mMcwY3.fasta --id 0.97 --w 12 --stepwords 20 --rev --usersort --maxaccepts 20 --stable_sort --uc /Users/caporaso/code/qiime/qiime_test_data/pick_open_reference_otus/ucrss//step2_otus/subsampled_failures_clusters.uc
# version=1.2.22
# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
S	0	143	*	*	*	*	*	QiimeExactMatch.t1_53967	*
S	1	137	*	*	*	*	*	QiimeExactMatch.t2_52630	*
S	2	152	*	*	*	*	*	QiimeExactMatch.f2_4913	*
H	2	131	99.2	+	0	0	131M21I	QiimeExactMatch.p2_32628	QiimeExactMatch.f2_4913
S	3	152	*	*	*	*	*	QiimeExactMatch.f3_141	*
H	0	118	100.0	+	0	0	118M25I	QiimeExactMatch.t1_52964	QiimeExactMatch.t1_53967
H	0	143	99.3	+	0	0	143M	QiimeExactMatch.t2_55465	QiimeExactMatch.t1_53967
C	0	3	99.7	*	*	*	*	QiimeExactMatch.t1_53967	*
C	1	1	*	*	*	*	*	QiimeExactMatch.t2_52630	*
C	2	2	99.2	*	*	*	*	QiimeExactMatch.f2_4913	*
C	3	1	*	*	*	*	*	QiimeExactMatch.f3_141	*
