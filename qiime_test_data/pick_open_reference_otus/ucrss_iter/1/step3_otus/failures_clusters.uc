# uclust --maxrejects 500 --input /Users/caporaso/temp/UclustExactMatchFilterqAk7PHnguOEkcDdAivG7.fasta --id 0.97 --w 12 --stepwords 20 --rev --usersort --maxaccepts 20 --libonly --stable_sort --lib /Users/caporaso/code/qiime/qiime_test_data/pick_open_reference_otus/ucrss_iter/1/step2_otus/step2_rep_set.fna --uc /Users/caporaso/code/qiime/qiime_test_data/pick_open_reference_otus/ucrss_iter//1//step3_otus/failures_clusters.uc
# version=1.2.22
# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
L	0	133	*	*	*	*	*	New.1.ReferenceOTU0 f5_1330	*
H	0	133	100.0	+	0	0	133M	QiimeExactMatch.f5_1300	New.1.ReferenceOTU0 f5_1330
D	0	2	*	*	*	*	100.0	New.1.ReferenceOTU0 f5_1330	*
