# uclust --maxrejects 500 --input /Users/caporaso/temp/UclustExactMatchFilterPsDnGaTTlonxXszz944q.fasta --id 0.97 --w 12 --stepwords 20 --rev --usersort --maxaccepts 20 --stable_sort --uc /Users/caporaso/code/qiime/qiime_test_data/pick_open_reference_otus/ucrss//step4_otus/failures_failures_clusters.uc
# version=1.2.22
# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
S	0	152	*	*	*	*	*	QiimeExactMatch.p1_38017	*
H	0	128	98.4	+	0	0	128M24I	QiimeExactMatch.p2_36441	QiimeExactMatch.p1_38017
S	1	131	*	*	*	*	*	QiimeExactMatch.p2_33268	*
C	0	2	98.4	*	*	*	*	QiimeExactMatch.p1_38017	*
C	1	1	*	*	*	*	*	QiimeExactMatch.p2_33268	*
