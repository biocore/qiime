#!/usr/bin/env python
#file test_filter_by_metadata.py
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project" 
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Release"

from cogent.util.unit_test import TestCase, main

from StringIO import StringIO
from qiime.pool_by_metadata import (pool_otu_table, pool_map)

class TopLevelTests(TestCase):
    """Tests of top-level functions"""

    def setUp(self):
        self.fasting_otu_table_str = \
"""#Full OTU Counts
#OTU ID	PC.354	PC.355	PC.356	PC.481	PC.593	PC.607	PC.634	PC.635	PC.636	Consensus Lineage
0	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
1	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
2	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Porphyromonadaceae;Parabacteroides
3	2	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
4	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
5	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
6	0	0	0	0	0	0	0	1	0	Root;Bacteria;Actinobacteria;Actinobacteria
7	0	0	2	0	0	0	0	0	2	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
8	1	1	0	2	4	0	0	0	0	Root;Bacteria;Firmicutes;"Bacilli";"Lactobacillales";Lactobacillaceae;Lactobacillus
9	0	0	2	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
10	0	1	0	0	0	0	0	0	0	Root;Bacteria
11	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
12	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes
13	1	0	0	1	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
14	0	0	1	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
15	0	0	0	0	1	0	0	0	0	Root;Bacteria
16	1	0	2	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
17	0	0	0	1	0	0	4	10	37	Root;Bacteria;Bacteroidetes
18	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
19	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes
20	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
21	0	0	0	0	0	0	2	3	2	Root;Bacteria;Bacteroidetes
22	0	0	0	0	2	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
23	14	1	14	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Bacilli";"Lactobacillales";Lactobacillaceae;Lactobacillus
24	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
25	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
26	0	0	0	0	0	0	0	1	1	Root;Bacteria;Bacteroidetes
27	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes
28	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
29	6	0	4	0	2	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
30	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes
31	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
32	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
33	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
34	0	0	0	0	0	0	8	10	2	Root;Bacteria
35	1	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
36	1	0	1	0	0	0	0	1	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
37	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
38	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
39	0	0	0	0	0	0	0	1	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
40	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
41	0	0	1	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
42	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes
43	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
44	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
45	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Erysipelotrichi";"Erysipelotrichales";Erysipelotrichaceae;Coprobacillus
46	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes
47	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
48	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
49	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
50	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
51	0	1	0	0	0	0	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
52	0	2	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
53	0	0	0	0	0	0	2	0	1	Root;Bacteria;Proteobacteria;Deltaproteobacteria
54	0	0	0	0	0	0	5	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Porphyromonadaceae;Parabacteroides
55	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
56	0	0	0	0	0	1	0	0	0	Root;Bacteria;Bacteroidetes
57	0	0	0	0	0	0	0	1	0	Root;Bacteria
58	1	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
59	0	0	0	0	0	0	0	0	1	Root;Bacteria;Deferribacteres;Deferribacteres;Deferribacterales;Deferribacteraceae;Mucispirillum
60	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
61	0	0	1	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
62	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
63	1	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
64	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
65	0	0	0	6	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
66	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
67	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
68	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
69	0	0	1	0	0	0	0	0	0	Root;Bacteria
70	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
71	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
72	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
73	0	0	0	0	0	5	0	0	0	Root;Bacteria;Bacteroidetes
74	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
75	1	0	1	0	0	0	0	0	0	Root;Bacteria;Bacteroidetes
76	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
77	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
78	1	0	1	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
79	2	3	8	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
80	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Porphyromonadaceae;Parabacteroides
81	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
82	0	0	0	0	0	2	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
83	0	0	0	1	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
84	1	0	0	0	0	0	0	2	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae";Ruminococcus
85	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
86	0	0	0	0	0	0	0	1	0	Root;Bacteria
87	0	0	1	0	0	2	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
88	0	0	0	0	0	0	0	1	0	Root;Bacteria
89	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
90	0	0	0	9	0	0	3	0	0	Root;Bacteria;Firmicutes;"Erysipelotrichi";"Erysipelotrichales";Erysipelotrichaceae;Turicibacter
91	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";Butyrivibrio
92	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
93	0	0	0	0	0	0	2	1	0	Root;Bacteria;Bacteroidetes
94	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
95	0	0	0	2	0	0	0	0	0	Root;Bacteria;Bacteroidetes
96	0	0	0	1	0	1	0	1	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
97	0	0	0	0	0	1	0	0	0	Root;Bacteria
98	0	0	0	0	0	0	0	1	0	Root;Bacteria
99	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
100	0	0	0	1	0	0	0	0	0	Root;Bacteria
101	0	0	0	3	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
102	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
103	0	1	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
104	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
105	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
106	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
107	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
108	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;Incertae Sedis XIII;Anaerovorax
109	0	0	0	1	0	0	1	5	2	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
110	0	0	0	0	0	2	0	0	0	Root;Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae;Olsenella
111	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
112	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
113	0	0	0	0	0	1	0	0	0	Root;Bacteria
114	0	0	0	0	0	1	0	0	0	Root;Bacteria
115	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes
116	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
117	1	0	2	0	0	6	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
118	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
119	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
120	1	3	1	2	1	9	2	4	5	Root;Bacteria;Bacteroidetes
121	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
122	0	0	0	1	0	2	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
123	0	0	0	0	0	0	1	0	0	Root;Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae
124	0	0	0	0	0	0	1	0	0	Root;Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae
125	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes
126	0	0	2	0	0	0	0	1	0	Root;Bacteria
127	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
128	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
129	0	0	0	1	0	0	0	0	0	Root;Bacteria
130	0	0	0	0	5	2	0	0	0	Root;Bacteria;Proteobacteria;Epsilonproteobacteria;Campylobacterales;Helicobacteraceae;Helicobacter
131	0	0	1	3	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
132	0	0	0	0	1	0	0	0	0	Root;Bacteria
133	0	0	1	0	0	0	0	0	0	Root;Bacteria
134	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
135	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
136	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
137	0	0	0	0	0	0	0	1	0	Root;Bacteria
138	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
139	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
140	0	0	0	0	0	0	1	3	0	Root;Bacteria
141	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
142	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
143	0	0	1	0	0	0	0	0	0	Root;Bacteria
144	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
145	0	0	2	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
146	1	0	0	0	2	0	2	0	3	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
147	0	1	0	1	1	0	0	0	3	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
148	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes
149	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
150	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
151	0	0	0	1	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
152	0	0	0	1	0	0	1	2	19	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
153	0	2	1	2	0	0	1	1	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
154	2	18	0	1	0	0	21	4	4	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
155	0	0	0	0	0	5	9	5	3	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
156	0	0	1	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
157	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
158	1	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
159	0	0	0	0	0	0	0	1	1	Root;Bacteria;Bacteroidetes
160	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
161	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
162	0	0	0	0	0	3	5	2	6	Root;Bacteria;Deferribacteres;Deferribacteres;Deferribacterales;Deferribacteraceae;Mucispirillum
163	0	0	0	0	0	0	0	0	1	Root;Bacteria
164	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
165	2	1	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
166	0	0	0	0	0	0	0	1	0	Root;Bacteria
167	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
168	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
169	0	2	0	7	0	0	0	2	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
170	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
171	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
172	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
173	0	0	0	0	0	1	0	0	0	Root;Bacteria
174	1	0	0	0	10	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Peptostreptococcaceae";"Peptostreptococcaceae Incertae Sedis"
175	0	0	0	0	1	0	0	0	0	Root;Bacteria;Bacteroidetes
176	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
177	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia"
178	0	0	0	2	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
179	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
180	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
181	1	4	2	6	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
182	0	0	0	0	0	1	0	0	0	Root;Bacteria
183	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia"
184	0	0	0	1	0	0	3	1	0	Root;Bacteria;Bacteroidetes
185	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
186	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
187	0	1	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
188	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
189	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
190	0	0	0	0	0	0	0	1	0	Root;Bacteria
191	2	1	10	2	24	0	0	1	1	Root;Bacteria
192	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Bacilli";"Lactobacillales";Streptococcaceae;Streptococcus
193	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";Butyrivibrio
194	0	0	2	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae";Acetanaerobacterium
195	0	0	0	0	0	1	0	0	0	Root;Bacteria
196	0	0	0	0	0	1	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
197	0	1	0	0	0	0	0	0	0	Root;Bacteria
198	0	2	0	0	0	1	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales
199	0	0	0	0	0	1	1	0	0	Root;Bacteria
200	0	0	0	2	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
201	0	0	0	1	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
202	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
203	0	2	2	4	0	5	1	5	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
204	1	4	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
205	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
206	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
207	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
208	0	2	0	2	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
209	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
210	0	0	0	0	0	0	0	0	1	Root;Bacteria
211	1	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
212	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
213	0	0	0	0	0	0	0	2	0	Root;Bacteria;Firmicutes
214	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
215	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
216	0	0	0	0	0	0	0	1	0	Root;Bacteria;Bacteroidetes
217	0	0	0	0	0	2	0	1	0	Root;Bacteria
218	0	0	0	0	9	1	0	0	0	Root;Bacteria;Bacteroidetes
219	0	0	0	0	1	0	0	0	0	Root;Bacteria
220	1	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
221	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes
222	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
223	0	0	0	0	0	0	0	2	2	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
224	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
225	0	2	1	0	0	0	0	0	0	Root;Bacteria;Bacteroidetes
226	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
227	0	1	2	0	9	1	1	1	3	Root;Bacteria;Bacteroidetes
228	16	0	0	0	12	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
229	0	0	0	0	0	1	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;Incertae Sedis XIII
230	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
231	0	19	2	0	2	0	3	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
232	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
233	0	0	0	0	1	0	0	0	0	Root;Bacteria;Bacteroidetes
234	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Bacilli";"Lactobacillales";Lactobacillaceae;Lactobacillus
235	0	1	1	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
236	0	0	0	0	0	2	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales
237	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
238	0	0	0	0	0	0	0	1	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
239	0	0	0	0	0	1	0	0	0	Root;Bacteria
240	0	0	0	0	0	1	0	0	0	Root;Bacteria
241	0	0	0	0	0	0	2	0	0	Root;Bacteria;TM7;TM7_genera_incertae_sedis
242	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
243	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
244	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes
245	0	0	0	1	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
246	0	0	0	0	0	0	0	1	0	Root;Bacteria
247	0	0	1	0	0	0	0	0	0	Root;Bacteria;Bacteroidetes
248	1	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Bacilli";"Lactobacillales";Lactobacillaceae;Lactobacillus
249	1	0	0	0	0	0	0	0	0	Root;Bacteria
250	1	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
251	0	0	0	1	4	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
252	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
253	0	0	0	0	2	0	0	5	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
254	11	13	6	13	2	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
255	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
256	0	0	0	0	0	0	1	0	0	Root;Bacteria
257	0	0	0	0	0	0	5	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
258	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
259	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
260	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
261	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
262	0	1	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";Bryantella
263	0	0	0	0	1	0	0	0	0	Root;Bacteria
264	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
265	0	0	0	0	0	2	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
266	0	0	0	2	0	0	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
267	1	0	0	5	17	20	0	0	0	Root;Bacteria
268	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
269	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
270	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
271	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
272	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
273	0	0	0	0	0	0	1	0	0	Root;Bacteria
274	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
275	0	0	0	0	0	0	1	0	0	Root;Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Verrucomicrobiaceae;Akkermansia
276	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
277	1	0	0	0	0	0	0	0	0	Root;Bacteria
278	0	0	0	0	0	1	0	0	0	Root;Bacteria
279	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
280	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
281	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
282	0	0	0	0	0	0	2	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Porphyromonadaceae;Parabacteroides
283	0	0	0	0	0	0	2	1	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
284	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
285	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
286	0	2	3	1	4	0	5	0	4	Root;Bacteria;Bacteroidetes
287	0	0	0	0	0	0	1	1	1	Root;Bacteria;Bacteroidetes
288	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
289	0	0	0	0	3	0	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
290	0	0	0	0	0	0	0	0	2	Root;Bacteria;Firmicutes;"Bacilli";Bacillales;"Staphylococcaceae";Staphylococcus
291	0	0	0	0	1	0	0	0	0	Root;Bacteria
292	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
293	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
294	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
295	29	1	10	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
296	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
297	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
298	0	0	0	0	0	0	1	0	0	Root;Bacteria;Actinobacteria;Actinobacteria
299	0	0	0	0	0	0	1	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
300	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia"
301	0	0	0	0	0	0	2	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
302	0	0	0	0	0	1	0	0	0	Root;Bacteria
303	0	0	0	0	0	0	0	0	1	Root;Bacteria
304	0	0	0	0	0	0	0	1	0	Root;Bacteria;Bacteroidetes
305	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
306	0	0	0	0	0	0	0	0	1	Root;Bacteria
307	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
308	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae";"Ruminococcaceae Incertae Sedis"
309	0	0	0	1	0	0	0	0	0	Root;Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae;Denitrobacterium
310	0	0	1	0	0	0	0	0	0	Root;Bacteria
311	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
312	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
313	0	1	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Porphyromonadaceae;Parabacteroides
314	0	0	1	0	0	0	0	0	0	Root;Bacteria;Bacteroidetes
315	1	3	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
316	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
317	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
318	0	0	0	0	0	1	0	0	0	Root;Bacteria;Proteobacteria
319	0	2	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
320	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
321	0	0	0	0	0	0	0	0	1	Root;Bacteria
322	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
323	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
324	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
325	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
326	0	0	0	0	4	0	0	0	2	Root;Bacteria;Firmicutes;"Erysipelotrichi";"Erysipelotrichales";Erysipelotrichaceae;Erysipelotrichaceae Incertae Sedis
327	0	0	0	0	0	0	0	1	0	Root;Bacteria;Bacteroidetes
328	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
329	2	2	0	1	0	0	0	0	0	Root;Bacteria;Bacteroidetes
330	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes
331	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes
332	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
333	0	0	0	0	0	6	0	3	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
334	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
335	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
336	0	0	1	0	0	0	0	0	0	Root;Bacteria
337	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
338	0	0	0	0	0	0	0	1	0	Root;Bacteria
339	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
340	0	0	2	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
341	0	0	1	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
342	0	0	0	0	0	1	0	0	0	Root;Bacteria
343	0	0	0	0	0	0	0	0	1	Root;Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae
344	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
345	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
346	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
347	0	0	0	1	0	0	0	0	0	Root;Bacteria
348	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
349	0	0	0	0	0	0	1	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
350	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
351	0	0	0	0	2	2	1	4	1	Root;Bacteria;Bacteroidetes
352	3	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
353	0	4	4	0	1	2	0	2	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
354	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
355	0	0	0	0	0	0	0	1	0	Root;Bacteria
356	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
357	0	0	0	4	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
358	0	0	1	0	0	0	0	0	0	Root;Bacteria
359	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
360	0	0	1	0	0	0	0	1	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
361	2	0	2	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
362	1	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
363	0	0	0	0	0	1	0	1	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae
364	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
365	0	0	0	0	0	2	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
366	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";Roseburia
367	0	0	0	0	1	0	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
368	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
369	0	0	0	0	0	1	0	0	0	Root;Bacteria
370	2	1	0	5	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
371	1	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
372	0	1	0	0	0	0	0	0	0	Root;Bacteria
373	0	1	0	0	0	0	3	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;Clostridiaceae;"Clostridiaceae 1";Clostridium
374	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
375	0	0	0	0	0	0	4	0	0	Root;Bacteria;Firmicutes;"Erysipelotrichi";"Erysipelotrichales";Erysipelotrichaceae;Erysipelotrichaceae Incertae Sedis
376	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
377	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
378	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes
379	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
380	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Bacilli";Bacillales;"Staphylococcaceae";Staphylococcus
381	0	0	2	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
382	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
383	4	9	0	2	0	0	0	2	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
384	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
385	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Bacilli";"Lactobacillales";"Carnobacteriaceae";"Carnobacteriaceae 1"
386	0	0	1	0	0	0	0	0	0	Root;Bacteria
387	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
388	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
389	0	1	0	0	0	0	0	0	0	Root;Bacteria
390	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
391	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes
392	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
393	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
394	0	0	1	0	0	0	0	0	0	Root;Bacteria
395	1	1	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
396	2	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
397	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
398	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
399	0	0	0	0	0	0	13	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
400	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
401	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
402	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
403	0	0	0	0	0	0	0	1	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Prevotellaceae
404	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
405	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
406	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
407	1	0	0	0	0	4	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
408	1	5	3	2	0	0	0	0	1	Root;Bacteria;Bacteroidetes
409	0	0	0	0	0	0	0	1	1	Root;Bacteria;Bacteroidetes
410	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
411	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
412	0	0	0	0	2	0	0	0	0	Root;Bacteria;Bacteroidetes
413	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
414	1	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
415	0	0	0	0	0	7	0	2	2	Root;Bacteria;Bacteroidetes
416	0	1	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales"""
        self.fasting_map_str = \
"""#SampleID	BarcodeSequence	LinkerPrimerSequence	Treatment	DOB	Description
#Example mapping file for the QIIME analysis package.  These 9 samples are from a study of the effects of exercise and diet on mouse cardiac physiology (Crawford, et al, PNAS, 2009).
PC.354	AGCACGAGCCTA	YATGCTGCCTCCCGTAGGAGT	Control	20061218	Control_mouse_I.D._354
PC.355	AACTCGTCGATG	YATGCTGCCTCCCGTAGGAGT	Control	20061218	Control_mouse_I.D._355
PC.356	ACAGACCACTCA	YATGCTGCCTCCCGTAGGAGT	Control	20061126	Control_mouse_I.D._356
PC.481	ACCAGCGACTAG	YATGCTGCCTCCCGTAGGAGT	Control	20070314	Control_mouse_I.D._481
PC.593	AGCAGCACTTGT	YATGCTGCCTCCCGTAGGAGT	Control	20071210	Control_mouse_I.D._593
PC.607	AACTGTGCGTAC	YATGCTGCCTCCCGTAGGAGT	Fast	20071112	Fasting_mouse_I.D._607
PC.634	ACAGAGTCGGCT	YATGCTGCCTCCCGTAGGAGT	Fast	20080116	Fasting_mouse_I.D._634
PC.635	ACCGCAGAGTCA	YATGCTGCCTCCCGTAGGAGT	Fast	20080116	Fasting_mouse_I.D._635
PC.636	ACGGTGAGTGTC	YATGCTGCCTCCCGTAGGAGT	Fast	20080116	Fasting_mouse_I.D._636
"""

    def test_pool_otu_table(self):
        """pool_otu_table should eliminate samples and add new one """
        def get_result():
            otu_infile = StringIO(self.fasting_otu_table_str)
            otu_outfile = StringIO('')
            pool_otu_table(otu_infile, otu_outfile,
                'pooledguy', ['PC.355', 'PC.481'])
            otu_outfile.seek(0)
            return otu_outfile.read()
        otu_res = get_result()

        assert('PC.355' not in otu_res)
        assert('pooledguy' in otu_res)

    def test_pool_map(self):
        """pool_map should eliminate samples and add new one """
        def get_result():
            map_infile = StringIO(self.fasting_map_str)
            map_outfile = StringIO('')
            pool_map(map_infile, map_outfile,
                'pooledguy', ['PC.355', 'PC.481'])
            map_outfile.seek(0)
            return map_outfile.read()
        map_res = get_result()

        assert('PC.355' not in map_res)
        assert('pooledguy' in map_res)



#run tests if called from command line
if __name__ == "__main__":
    main()
