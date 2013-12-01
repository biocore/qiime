# uclust --maxrejects 500 --input /Users/caporaso/temp/UclustExactMatchFilterdiKUoTREeKB2CgB2Ao0l.fasta --id 0.6 --w 12 --stepwords 20 --rev --usersort --maxaccepts 20 --libonly --stable_sort --lib /Users/caporaso/code/qiime/qiime_test_data/pick_open_reference_otus/refseqs.fna --uc /Users/caporaso/code/qiime/qiime_test_data/pick_open_reference_otus/ucrss_iter_no_tax//1//prefilter_otus/seqs2_clusters.uc
# version=1.2.22
# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
L	3	1389	*	*	*	*	*	295053	*
H	3	133	100.0	+	0	0	519I133M737I	QiimeExactMatch.f5_1370	295053
H	3	133	96.2	+	0	0	519I133M737I	QiimeExactMatch.f5_1300	295053
L	0	1486	*	*	*	*	*	426848	*
H	0	152	78.9	+	0	0	519I93MI59M814I	QiimeExactMatch.f5_1360	426848
H	3	143	100.0	+	0	0	519I143M727I	QiimeExactMatch.f5_1470	295053
H	3	141	100.0	+	0	0	519I141M729I	QiimeExactMatch.f5_1510	295053
H	0	152	77.6	+	0	0	519I93MI59M814I	QiimeExactMatch.f5_1400	426848
D	0	3	*	*	*	*	78.3	426848	*
D	3	5	*	*	*	*	99.1	295053	*
