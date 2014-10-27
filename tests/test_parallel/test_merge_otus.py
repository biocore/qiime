#!/usr/bin/env python

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Daniel McDonald", "Jai Ram Rideout",
               "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"

from unittest import TestCase, main
from os import close
from os.path import exists, join
from shutil import rmtree
from tempfile import mkdtemp, mkstemp

from biom.parse import parse_biom_table
from biom.util import biom_open
from skbio.util.misc import remove_files

from qiime.util import get_qiime_temp_dir
from qiime.test import initiate_timeout, disable_timeout
from qiime.parallel.merge_otus import ParallelMergeOtus


class MergeTests(TestCase):

    def setUp(self):
        self._files_to_remove = []
        self._dirs_to_remove = []

        tmp_dir = get_qiime_temp_dir()
        self.test_out = mkdtemp(dir=tmp_dir,
                                prefix='qiime_parallel_tests_',
                                suffix='')
        self._dirs_to_remove.append(self.test_out)

        fd, self.t1_fp = mkstemp(dir=self.test_out, prefix='t1_',
                                 suffix='.biom')
        close(fd)
        with open(self.t1_fp, 'w') as f:
            f.write(t1)
        self._files_to_remove.append(self.t1_fp)

        fd, self.t2_fp = mkstemp(dir=self.test_out, prefix='t2_',
                                 suffix='.biom')
        close(fd)
        with open(self.t2_fp, 'w') as f:
            f.write(t2)
        self._files_to_remove.append(self.t2_fp)

        fd, self.t3_fp = mkstemp(dir=self.test_out, prefix='t3_',
                                 suffix='.biom')
        close(fd)
        with open(self.t3_fp, 'w') as f:
            f.write(t3)
        self._files_to_remove.append(self.t3_fp)

        fd, self.t4_fp = mkstemp(dir=self.test_out, prefix='t4_',
                                 suffix='.biom')
        close(fd)
        with open(self.t4_fp, 'w') as f:
            f.write(t4)
        self._files_to_remove.append(self.t4_fp)

        fd, exp_even_fp = mkstemp(dir=self.test_out, prefix='exp_even_',
                                  suffix='.biom')
        close(fd)
        with open(exp_even_fp, 'w') as f:
            f.write(exp_even)
        self._files_to_remove.append(exp_even_fp)
        with biom_open(exp_even_fp) as f:
            self.exp_even = parse_biom_table(f)

        fd, exp_odd_fp = mkstemp(dir=self.test_out, prefix='exp_odd_',
                                 suffix='.biom')
        close(fd)
        with open(exp_odd_fp, 'w') as f:
            f.write(exp_odd)
        self._files_to_remove.append(exp_odd_fp)
        with biom_open(exp_odd_fp) as f:
            self.exp_odd = parse_biom_table(f)

        initiate_timeout(60)

    def tearDown(self):
        disable_timeout()
        remove_files(self._files_to_remove)
        # remove directories last, so we don't get errors trying to remove
        # files which may be in the directories
        for d in self._dirs_to_remove:
            if exists(d):
                rmtree(d)

    def test_parallel_merge_otus_even(self):
        infiles = [self.t1_fp, self.t2_fp, self.t3_fp, self.t4_fp]
        app = ParallelMergeOtus()
        app(infiles, self.test_out)

        # Confirm that the output OTU table have been merged correctly
        merged_fp = join(self.test_out, "merged.biom")
        self.assertTrue(exists(merged_fp))
        with biom_open(merged_fp) as f:
            obs = parse_biom_table(f)
        self.assertEqual(obs, self.exp_even)

    def test_parallel_merge_otus_odd(self):
        infiles = [self.t1_fp, self.t2_fp, self.t3_fp]
        app = ParallelMergeOtus()
        app(infiles, self.test_out)

        # Confirm that the output OTU table have been merged correctly
        merged_fp = join(self.test_out, "merged.biom")
        self.assertTrue(exists(merged_fp))
        with biom_open(merged_fp) as f:
            obs = parse_biom_table(f)
        self.assertEqual(obs, self.exp_odd)

t1 = """{
    "id":null,
    "format": "Biological Observation Matrix 1.0.0-dev",
    "format_url": "http://biom-format.org",
    "type": "OTU table",
    "generated_by": "QIIME revision XYZ",
    "date": "2011-12-19T19:00:00",
    "rows":[
            {"id":"GG_OTU_1", "metadata":null},
            {"id":"GG_OTU_2", "metadata":null},
            {"id":"GG_OTU_3", "metadata":null},
            {"id":"GG_OTU_4", "metadata":null},
            {"id":"GG_OTU_5", "metadata":null}
        ],
    "columns": [
            {"id":"Sample1", "metadata":null},
            {"id":"Sample2", "metadata":null},
            {"id":"Sample3", "metadata":null},
            {"id":"Sample4", "metadata":null},
            {"id":"Sample5", "metadata":null},
            {"id":"Sample6", "metadata":null}
        ],
    "matrix_type": "sparse",
    "matrix_element_type": "int",
    "shape": [5, 6],
    "data":[[0,2,1],
            [1,0,5],
            [1,1,1],
            [1,3,2],
            [1,4,3],
            [1,5,1],
            [2,2,1],
            [2,3,4],
            [2,5,2],
            [3,0,2],
            [3,1,1],
            [3,2,1],
            [3,5,1],
            [4,1,1],
            [4,2,1]
           ]
}"""

t2 = """{
     "id":null,
     "format": "Biological Observation Matrix 1.0.0-dev",
     "format_url": "http://biom-format.org",
     "type": "OTU table",
     "generated_by": "QIIME revision XYZ",
     "date": "2011-12-19T19:00:00",
     "rows":[
        {"id":"GG_OTU_1", "metadata":{"taxonomy":["k__Bacteria", "p__Proteobac\
teria", "c__Gammaproteobacteria", "o__Enterobacteriales", "f__Enterobacteriace\
ae", "g__Escherichia", "s__"]}},
        {"id":"GG_OTU_20", "metadata":{"taxonomy":["k__Bacteria", "p__Cyanobac\
teria", "c__Nostocophycideae", "o__Nostocales", "f__Nostocaceae", "g__Dolichos\
permum", "s__"]}},
        {"id":"GG_OTU_3", "metadata":{"taxonomy":["k__Archaea", "p__Euryarchae\
ota", "c__Methanomicrobia", "o__Methanosarcinales", "f__Methanosarcinaceae", "\
g__Methanosarcina", "s__"]}},
        {"id":"GG_OTU_44", "metadata":{"taxonomy":["k__Bacteria", "p__Firmicut\
es", "c__Clostridia", "o__Halanaerobiales", "f__Halanaerobiaceae", "g__Halanae\
robium", "s__Halanaerobiumsaccharolyticum"]}},
        {"id":"GG_OTU_5", "metadata":{"taxonomy":["k__Bacteria", "p__Proteobac\
teria", "c__Gammaproteobacteria", "o__Enterobacteriales", "f__Enterobacteriace\
ae", "g__Escherichia", "s__"]}}
        ],
     "columns":[
        {"id":"Sample1", "metadata":{
                                "BarcodeSequence":"CGCTTATCGAGA",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"gut",
                                "Description":"human gut"}},
        {"id":"Sample2", "metadata":{
                                "BarcodeSequence":"CATACCAGTAGC",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"gut",
                                "Description":"human gut"}},
        {"id":"Sample3", "metadata":{
                                "BarcodeSequence":"CTCTCTACCTGT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"gut",
                                "Description":"human gut"}},
        {"id":"Sample4", "metadata":{
                                "BarcodeSequence":"CTCTCGGCCTGT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"skin",
                                "Description":"human skin"}},
        {"id":"Sample5", "metadata":{
                                "BarcodeSequence":"CTCTCTACCAAT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"skin",
                                "Description":"human skin"}},
        {"id":"Sample6", "metadata":{
                                "BarcodeSequence":"CTAACTACCAAT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"skin",
                                "Description":"human skin"}}
                ],
     "matrix_type": "dense",
     "matrix_element_type": "int",
     "shape": [5,6],
     "data":  [[0,0,1,0,0,0],
               [5,1,0,2,3,1],
               [0,0,1,4,2,0],
               [2,1,1,0,0,1],
               [0,1,1,0,0,0]]
}"""

t3 = """{
        "id":null,
        "format": "Biological Observation Matrix 1.0.0-dev",
        "format_url": "http://biom-format.org",
        "type": "OTU table",
        "generated_by": "QIIME revision XYZ",
        "date": "2011-12-19T19:00:00",
        "rows":[
                {"id":"GG_OTU_1", "metadata":null},
                {"id":"GG_OTU_2", "metadata":null},
                {"id":"GG_OTU_3", "metadata":null},
                {"id":"GG_OTU_4", "metadata":null},
                {"id":"GG_OTU_5", "metadata":null}
            ],
        "columns": [
                {"id":"Sample1", "metadata":null},
                {"id":"Sample2", "metadata":null},
                {"id":"Sample3", "metadata":null},
                {"id":"Sample4", "metadata":null},
                {"id":"Sample5", "metadata":null},
                {"id":"Sample6", "metadata":null}
            ],
        "matrix_type": "sparse",
        "matrix_element_type": "int",
        "shape": [5, 6],
        "data":[[0,2,1],
                [1,0,5],
                [1,1,1],
                [1,3,2],
                [1,4,3],
                [1,5,1],
                [2,2,1],
                [2,3,4],
                [2,5,2],
                [3,0,2],
                [3,1,1],
                [3,2,1],
                [3,5,1],
                [4,1,1],
                [4,2,1]
               ]
}"""

t4 = """{
     "id":null,
     "format": "Biological Observation Matrix 1.0.0-dev",
     "format_url": "http://biom-format.org",
     "type": "OTU table",
     "generated_by": "QIIME revision XYZ",
     "date": "2011-12-19T19:00:00",
     "rows":[
        {"id":"GG_OTU_1", "metadata":{"taxonomy":["k__Bacteria", "p__Proteobac\
teria", "c__Gammaproteobacteria", "o__Enterobacteriales", "f__Enterobacteriace\
ae", "g__Escherichia", "s__"]}},
        {"id":"GG_OTU_20", "metadata":{"taxonomy":["k__Bacteria", "p__Cyanobac\
teria", "c__Nostocophycideae", "o__Nostocales", "f__Nostocaceae", "g__Dolichos\
permum", "s__"]}},
        {"id":"GG_OTU_3", "metadata":{"taxonomy":["k__Archaea", "p__Euryarchae\
ota", "c__Methanomicrobia", "o__Methanosarcinales", "f__Methanosarcinaceae", "\
g__Methanosarcina", "s__"]}},
        {"id":"GG_OTU_44", "metadata":{"taxonomy":["k__Bacteria", "p__Firmicut\
es", "c__Clostridia", "o__Halanaerobiales", "f__Halanaerobiaceae", "g__Halanae\
robium", "s__Halanaerobiumsaccharolyticum"]}},
        {"id":"GG_OTU_5", "metadata":{"taxonomy":["k__Bacteria", "p__Proteobac\
teria", "c__Gammaproteobacteria", "o__Enterobacteriales", "f__Enterobacteriace\
ae", "g__Escherichia", "s__"]}}
        ],
     "columns":[
        {"id":"Sample1", "metadata":{
                                "BarcodeSequence":"CGCTTATCGAGA",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"gut",
                                "Description":"human gut"}},
        {"id":"Sample2", "metadata":{
                                "BarcodeSequence":"CATACCAGTAGC",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"gut",
                                "Description":"human gut"}},
        {"id":"Sample3", "metadata":{
                                "BarcodeSequence":"CTCTCTACCTGT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"gut",
                                "Description":"human gut"}},
        {"id":"Sample4", "metadata":{
                                "BarcodeSequence":"CTCTCGGCCTGT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"skin",
                                "Description":"human skin"}},
        {"id":"Sample5", "metadata":{
                                "BarcodeSequence":"CTCTCTACCAAT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"skin",
                                "Description":"human skin"}},
        {"id":"Sample6", "metadata":{
                                "BarcodeSequence":"CTAACTACCAAT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"skin",
                                "Description":"human skin"}}
                ],
     "matrix_type": "dense",
     "matrix_element_type": "int",
     "shape": [5,6],
     "data":  [[0,0,1,0,0,0],
               [5,1,0,2,3,1],
               [0,0,1,4,2,0],
               [2,1,1,0,0,1],
               [0,1,1,0,0,0]]
}"""

exp_even = """{"id": "No Table ID","format": "Biological Observation Matrix 1.\
0.0","format_url": "http://biom-format.org","matrix_type": "sparse","generated\
_by": "BIOM-Format 2.0.1-dev","date": "2014-07-22T12:11:47.418442","type": nul\
l,"matrix_element_type": "float","shape": [7, 6],"data": [[0,2,4.0],[1,\
0,10.0],[1,1,2.0],[1,3,4.0],[1,4,6.0],[1,5,2.0],[2,2,4.0],[2,3,16.0],[2,4,4.0]\
,[2,5,4.0],[3,0,4.0],[3,1,2.0],[3,2,2.0],[3,5,2.0],[4,1,4.0],[4,2,4.0],[5,0,10\
.0],[5,1,2.0],[5,3,4.0],[5,4,6.0],[5,5,2.0],[6,0,4.0],[6,1,2.0],[6,2,2.0],[6,5\
,2.0]],"rows": [{"id": "GG_OTU_1", "metadata": {"taxonomy": ["k__Bacteria", "p\
__Proteobacteria", "c__Gammaproteobacteria", "o__Enterobacteriales", "f__Enter\
obacteriaceae", "g__Escherichia", "s__"]}},{"id": "GG_OTU_2", "metadata": {"ta\
xonomy": null}},{"id": "GG_OTU_3", "metadata": {"taxonomy": ["k__Archaea", "p_\
_Euryarchaeota", "c__Methanomicrobia", "o__Methanosarcinales", "f__Methanosarc\
inaceae", "g__Methanosarcina", "s__"]}},{"id": "GG_OTU_4", "metadata": {"taxon\
omy": null}},{"id": "GG_OTU_5", "metadata": {"taxonomy": ["k__Bacteria", "p__P\
roteobacteria", "c__Gammaproteobacteria", "o__Enterobacteriales", "f__Enteroba\
cteriaceae", "g__Escherichia", "s__"]}},{"id": "GG_OTU_20", "metadata": {"taxo\
nomy": ["k__Bacteria", "p__Cyanobacteria", "c__Nostocophycideae", "o__Nostocal\
es", "f__Nostocaceae", "g__Dolichospermum", "s__"]}},{"id": "GG_OTU_44", "meta\
data": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Hala\
naerobiales", "f__Halanaerobiaceae", "g__Halanaerobium", "s__Halanaerobiumsacc\
harolyticum"]}}],"columns": [{"id": "Sample1", "metadata": {"LinkerPrimerSeque\
nce": "CATGCTGCCTCCCGTAGGAGT", "BarcodeSequence": "CGCTTATCGAGA", "Description\
": "human gut", "BODY_SITE": "gut"}},{"id": "Sample2", "metadata": {"LinkerPri\
merSequence": "CATGCTGCCTCCCGTAGGAGT", "BarcodeSequence": "CATACCAGTAGC", "Des\
cription": "human gut", "BODY_SITE": "gut"}},{"id": "Sample3", "metadata": {"L\
inkerPrimerSequence": "CATGCTGCCTCCCGTAGGAGT", "BarcodeSequence": "CTCTCTACCTG\
T", "Description": "human gut", "BODY_SITE": "gut"}},{"id": "Sample4", "metada\
ta": {"LinkerPrimerSequence": "CATGCTGCCTCCCGTAGGAGT", "BarcodeSequence": "CTC\
TCGGCCTGT", "Description": "human skin", "BODY_SITE": "skin"}},{"id": "Sample5\
", "metadata": {"LinkerPrimerSequence": "CATGCTGCCTCCCGTAGGAGT", "BarcodeSeque\
nce": "CTCTCTACCAAT", "Description": "human skin", "BODY_SITE": "skin"}},{"id"\
: "Sample6", "metadata": {"LinkerPrimerSequence": "CATGCTGCCTCCCGTAGGAGT", "Ba\
rcodeSequence": "CTAACTACCAAT", "Description": "human skin", "BODY_SITE": "ski\
n"}}]}"""

exp_odd = """{"id": "No Table ID","format": "Biological Observation Matrix 1.0\
.0","format_url": "http://biom-format.org","matrix_type": "sparse","generated_\
by": "BIOM-Format 2.0.1-dev","date": "2014-07-22T12:22:36.723558","type": null\
,"matrix_element_type": "float","shape": [7, 6],"data": [[0,2,3.0],[1,0\
,10.0],[1,1,2.0],[1,3,4.0],[1,4,6.0],[1,5,2.0],[2,2,3.0],[2,3,12.0],[2,4,2.0],\
[2,5,4.0],[3,0,4.0],[3,1,2.0],[3,2,2.0],[3,5,2.0],[4,1,3.0],[4,2,3.0],[5,0,5.0\
],[5,1,1.0],[5,3,2.0],[5,4,3.0],[5,5,1.0],[6,0,2.0],[6,1,1.0],[6,2,1.0],[6,5,1\
.0]],"rows": [{"id": "GG_OTU_1", "metadata": {"taxonomy": ["k__Bacteria", "p__\
Proteobacteria", "c__Gammaproteobacteria", "o__Enterobacteriales", "f__Enterob\
acteriaceae", "g__Escherichia", "s__"]}},{"id": "GG_OTU_2", "metadata": {"taxo\
nomy": null}},{"id": "GG_OTU_3", "metadata": {"taxonomy": ["k__Archaea", "p__E\
uryarchaeota", "c__Methanomicrobia", "o__Methanosarcinales", "f__Methanosarcin\
aceae", "g__Methanosarcina", "s__"]}},{"id": "GG_OTU_4", "metadata": {"taxonom\
y": null}},{"id": "GG_OTU_5", "metadata": {"taxonomy": ["k__Bacteria", "p__Pro\
teobacteria", "c__Gammaproteobacteria", "o__Enterobacteriales", "f__Enterobact\
eriaceae", "g__Escherichia", "s__"]}},{"id": "GG_OTU_20", "metadata": {"taxono\
my": ["k__Bacteria", "p__Cyanobacteria", "c__Nostocophycideae", "o__Nostocales\
", "f__Nostocaceae", "g__Dolichospermum", "s__"]}},{"id": "GG_OTU_44", "metada\
ta": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Halana\
erobiales", "f__Halanaerobiaceae", "g__Halanaerobium", "s__Halanaerobiumsaccha\
rolyticum"]}}],"columns": [{"id": "Sample1", "metadata": {"LinkerPrimerSequenc\
e": "CATGCTGCCTCCCGTAGGAGT", "BarcodeSequence": "CGCTTATCGAGA", "Description":\
 "human gut", "BODY_SITE": "gut"}},{"id": "Sample2", "metadata": {"LinkerPrime\
rSequence": "CATGCTGCCTCCCGTAGGAGT", "BarcodeSequence": "CATACCAGTAGC", "Descr\
iption": "human gut", "BODY_SITE": "gut"}},{"id": "Sample3", "metadata": {"Lin\
kerPrimerSequence": "CATGCTGCCTCCCGTAGGAGT", "BarcodeSequence": "CTCTCTACCTGT"\
, "Description": "human gut", "BODY_SITE": "gut"}},{"id": "Sample4", "metadata\
": {"LinkerPrimerSequence": "CATGCTGCCTCCCGTAGGAGT", "BarcodeSequence": "CTCTC\
GGCCTGT", "Description": "human skin", "BODY_SITE": "skin"}},{"id": "Sample5",\
 "metadata": {"LinkerPrimerSequence": "CATGCTGCCTCCCGTAGGAGT", "BarcodeSequenc\
e": "CTCTCTACCAAT", "Description": "human skin", "BODY_SITE": "skin"}},{"id": \
"Sample6", "metadata": {"LinkerPrimerSequence": "CATGCTGCCTCCCGTAGGAGT", "Barc\
odeSequence": "CTAACTACCAAT", "Description": "human skin", "BODY_SITE": "skin"\
}}]}"""

if __name__ == '__main__':
    main()
