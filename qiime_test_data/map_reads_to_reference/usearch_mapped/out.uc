# usearch --maxaccepts 1 --id 0.75 --queryalnfract 0.35 --blast6out /Users/caporaso/Dropbox/code/qiime_test_data/assign_reads_to_database/usearch_mapped/out.bl6 --uc /Users/caporaso/Dropbox/code/qiime_test_data/assign_reads_to_database/usearch_mapped/out.uc --query /Users/caporaso/Dropbox/code/qiime_test_data/assign_reads_to_database/query_nt.fasta --maxrejects 32 --evalue 1e-10 --db /Users/caporaso/Dropbox/code/qiime_test_data/assign_reads_to_database/refseqs_pr.fasta --targetalnfract 0.0
# version=5.2.32
# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
N	*	21	*	*	*	*	*	s1_1	*
N	*	22	*	*	*	*	*	s1_1	*
N	*	21	*	*	*	*	*	s1_1	*
N	*	66	*	*	*	*	*	s2_2 some comments...	*
N	*	23	*	*	*	*	*	s2_2 some comments...	*
N	*	120	*	*	*	*	*	s2_2 some comments...	*
N	*	61	*	*	*	*	*	s2_2 some comments...	*
H	0	376	100.0	.	0	0	376M	s2_2 some comments...	eco:b0015-pr dnaJ
N	*	63	*	*	*	*	*	s2_2 some comments...	*
N	*	27	*	*	*	*	*	s2_2 some comments...	*
N	*	90	*	*	*	*	*	s2_2 some comments...	*
H	1	115	100.0	.	0	0	115M	s1_3	eco:b0122-pr
N	*	25	*	*	*	*	*	s1_3	*
H	1	116	98.3	.	0	0	115M	s1_4	eco:b0122-pr
N	*	25	*	*	*	*	*	s1_4	*
N	*	66	*	*	*	*	*	s1_5	*
N	*	23	*	*	*	*	*	s1_5	*
N	*	120	*	*	*	*	*	s1_5	*
N	*	61	*	*	*	*	*	s1_5	*
H	0	376	100.0	.	0	0	376M	s1_5	eco:b0015-pr dnaJ
N	*	63	*	*	*	*	*	s1_5	*
N	*	27	*	*	*	*	*	s1_5	*
N	*	90	*	*	*	*	*	s1_5	*
N	*	66	*	*	*	*	*	s1_6 some comments...	*
N	*	23	*	*	*	*	*	s1_6 some comments...	*
N	*	120	*	*	*	*	*	s1_6 some comments...	*
N	*	61	*	*	*	*	*	s1_6 some comments...	*
H	0	376	99.7	.	0	0	376M	s1_6 some comments...	eco:b0015-pr dnaJ
N	*	63	*	*	*	*	*	s1_6 some comments...	*
N	*	27	*	*	*	*	*	s1_6 some comments...	*
N	*	90	*	*	*	*	*	s1_6 some comments...	*
