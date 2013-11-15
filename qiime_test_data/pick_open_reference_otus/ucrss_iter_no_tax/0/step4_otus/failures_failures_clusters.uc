# uclust --maxrejects 500 --input /Users/caporaso/temp/UclustExactMatchFilterSNGE0qlhGzvv9AutBq0T.fasta --id 0.97 --w 12 --stepwords 20 --rev --usersort --maxaccepts 20 --stable_sort --uc /Users/caporaso/code/qiime/qiime_test_data/pick_open_reference_otus/ucrss_iter_no_tax//0//step4_otus/failures_failures_clusters.uc
# version=1.2.22
# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
S	0	152	*	*	*	*	*	QiimeExactMatch.f1_4866	*
H	0	152	98.7	+	0	0	152M	QiimeExactMatch.f1_8307	QiimeExactMatch.f1_4866
S	1	152	*	*	*	*	*	QiimeExactMatch.p1_38017	*
H	1	128	98.4	+	0	0	128M24I	QiimeExactMatch.p2_36441	QiimeExactMatch.p1_38017
H	0	140	98.6	+	0	0	140M12I	QiimeExactMatch.f2_10228	QiimeExactMatch.f1_4866
C	0	3	98.6	*	*	*	*	QiimeExactMatch.f1_4866	*
C	1	2	98.4	*	*	*	*	QiimeExactMatch.p1_38017	*
