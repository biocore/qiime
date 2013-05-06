# uclust --maxrejects 500 --input /Users/caporaso/temp/UclustExactMatchFilterGIDc0lt00d0SQD2fQWiX.fasta --id 0.97 --w 12 --stepwords 20 --rev --usersort --maxaccepts 20 --stable_sort --uc /Users/caporaso/code/Qiime/qiime_test_data/pick_open_reference_otus/ucrss//step4_otus/failures_failures_clusters.uc
# version=1.2.22
# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
S	0	137	*	*	*	*	*	QiimeExactMatch.p2_35810	*
S	1	119	*	*	*	*	*	QiimeExactMatch.p1_36300	*
H	0	133	100.0	+	0	0	133M4I	QiimeExactMatch.t1_53038	QiimeExactMatch.p2_35810
H	1	152	98.3	+	0	0	119M33D	QiimeExactMatch.f2_4913	QiimeExactMatch.p1_36300
H	1	131	99.2	+	0	0	119M12D	QiimeExactMatch.p2_32628	QiimeExactMatch.p1_36300
S	2	152	*	*	*	*	*	QiimeExactMatch.p1_38017	*
H	2	128	98.4	+	0	0	128M24I	QiimeExactMatch.p2_36441	QiimeExactMatch.p1_38017
S	3	131	*	*	*	*	*	QiimeExactMatch.p2_33268	*
S	4	136	*	*	*	*	*	QiimeExactMatch.p2_31582	*
C	0	2	100.0	*	*	*	*	QiimeExactMatch.p2_35810	*
C	1	3	98.7	*	*	*	*	QiimeExactMatch.p1_36300	*
C	2	2	98.4	*	*	*	*	QiimeExactMatch.p1_38017	*
C	3	1	*	*	*	*	*	QiimeExactMatch.p2_33268	*
C	4	1	*	*	*	*	*	QiimeExactMatch.p2_31582	*
