# uclust --maxrejects 500 --input /Users/caporaso/temp/UclustExactMatchFilterUMYk3pHYdVXFez1hEsWS.fasta --id 0.97 --w 12 --stepwords 20 --rev --usersort --maxaccepts 20 --stable_sort --uc /Users/caporaso/code/qiime/qiime_test_data/pick_open_reference_otus/ucrss_iter_no_tax//0//step2_otus/subsampled_failures_clusters.uc
# version=1.2.22
# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
S	0	137	*	*	*	*	*	QiimeExactMatch.p2_35810	*
S	1	119	*	*	*	*	*	QiimeExactMatch.p1_36300	*
H	1	152	98.3	+	0	0	119M33D	QiimeExactMatch.f2_4913	QiimeExactMatch.p1_36300
S	2	131	*	*	*	*	*	QiimeExactMatch.p2_33268	*
S	3	140	*	*	*	*	*	QiimeExactMatch.t1_56590	*
H	3	118	100.0	+	0	0	118M22I	QiimeExactMatch.t1_54579	QiimeExactMatch.t1_56590
C	0	1	*	*	*	*	*	QiimeExactMatch.p2_35810	*
C	1	2	98.3	*	*	*	*	QiimeExactMatch.p1_36300	*
C	2	1	*	*	*	*	*	QiimeExactMatch.p2_33268	*
C	3	2	100.0	*	*	*	*	QiimeExactMatch.t1_56590	*
