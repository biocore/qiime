#!/usr/bin/env python
#file test_filter_otu_table.py

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2011, The QIIME Project" #consider project name
__credits__ = ["Jesse Stombaugh"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Tony Walters"
__email__ = "William.A.Walters@colorado.edu"
__status__ = "Release"


from sys import argv
from string import strip
from shutil import rmtree

from cogent.util.unit_test import TestCase, main
from numpy import array
from qiime.util import get_tmp_filename
from cogent.util.misc import remove_files 

from qiime.filter_otu_table import (strip_quotes,split_tax,
                _filter_table_samples, filter_table, _filter_table_neg_control)

class TopLevelTests(TestCase):
    """Tests of top-level functions"""

    def setUp(self):
        """define some top-level data"""
        self.otu_fname='otu_table.txt'
        self.min_count=1
        self.min_samples=2
        self.include_taxonomy='Root;Bacteria'
        self.exclude_taxonomy='Root;Archaea'
        self.dir_path='./'
        self.taxon1='"Root;Bacteria"'
        self.taxon2='Root;Archaea'
        
        self.sample_unfiltered_otu_table = sample_unfiltered_otu_table
        self.sample_unfiltered_otu_table_no_taxa =\
         sample_unfiltered_otu_table_no_taxa
        
        self.sample_input_otu_table =\
         get_tmp_filename(prefix = "sample_otu_table", suffix = ".txt")
        seq_file = open(self.sample_input_otu_table, 'w')
        seq_file.write(self.sample_unfiltered_otu_table)
        seq_file.close()
        
        self.sample_input_otu_table_no_taxa =\
         get_tmp_filename(prefix = "sample_otu_table", suffix = ".txt")
        seq_file = open(self.sample_input_otu_table_no_taxa, 'w')
        seq_file.write(self.sample_unfiltered_otu_table_no_taxa)
        seq_file.close()
        
        self._files_to_remove =\
         [self.sample_input_otu_table, self.sample_input_otu_table_no_taxa]
         
        self.expected_otu_table_output_leniant =\
         expected_otu_table_output_leniant
        self.expected_otu_table_output_default =\
         expected_otu_table_output_default
        self.expected_otu_table_output_sequence_count_three =\
         expected_otu_table_output_sequence_count_three
        self.expected_otu_table_include_taxa =\
         expected_otu_table_include_taxa
        self.expected_otu_table_exclude_taxa =\
         expected_otu_table_exclude_taxa
        self.expected_otu_table_include_exclude_taxa =\
         expected_otu_table_include_exclude_taxa
        self.expected_otu_table_output_default_no_taxa =\
         expected_otu_table_output_default_no_taxa
         
    def tearDown(self):
        remove_files(self._files_to_remove)



    def test_strip_quotes(self):
        """strip_quotes: removes leading and trailing quotes from tax string"""

        self.sample_to_extract='SampleA,SampleB'
        exp1='Root;Bacteria'
        exp2='Root;Archaea'
        
        obs1=strip_quotes(self.taxon1)
        obs2=strip_quotes(self.taxon2)
        
        self.assertEqual(obs1,exp1)
        self.assertEqual(obs2,exp2)

    def test_split_tax(self):
        """split_tax: Splits the tax string on comma and semicolon"""

        exp1=['"Root','Bacteria"']
        exp2=['Root','Archaea']
        
        obs1=split_tax(self.taxon1)
        obs2=split_tax(self.taxon2)
        
        self.assertEqual(obs1,exp1)
        self.assertEqual(obs2,exp2)
       
    def test_filter_table_samples(self):
        """_filter_table_samples removes small samples from OTU table
        """

        otu_table = """# QIIME v%s OTU table\n#OTU ID\tsample1\tsample2\tsample3
0\t0\t2\t0
1\t1\t0\t0
2\t1\t1\t1""" % __version__
        otu_table = otu_table.split('\n')
        result = _filter_table_samples(otu_table, 2)
        self.assertEqual(result, "# QIIME v%s OTU table\n#OTU ID\tsample1\tsample2\n0\t0\t2\n1\t1\t0\n2\t1\t1" % __version__)
        result = _filter_table_samples(otu_table, 1)
        self.assertEqual(result, '\n'.join(otu_table))
        result = _filter_table_samples(otu_table, 3)
        self.assertEqual(result, "# QIIME v%s OTU table\n#OTU ID\tsample2\n0\t2\n1\t0\n2\t1" % __version__)
        
    def test_filter_table_lenient(self):
        """ filter_table does not remove any OTUs with lax settings """
        
        params = {'min_otu_count': 0, 'min_otu_samples': 0,
         'included_taxa': '', 'excluded_taxa': ''}
         
        otu_file = open(self.sample_input_otu_table, "U")
        
        filtered_otu_table_fp = get_tmp_filename(prefix = "filtered_otu_table_",
          suffix = ".txt")
        
        filtered_otu_table_f =\
         open(filtered_otu_table_fp, "w")
        
        
        filter_table(params, filtered_otu_table_f, otu_file)
        
        filtered_otu_table_f.close()
        
        # Output is the same as input otu file, except for header
        
        actual_result_f = open(filtered_otu_table_fp, "U")
        
        
        actual_results = "\n".join([line.strip() for line in actual_result_f])
        
        self.assertEqual(actual_results, self.expected_otu_table_output_leniant)
        
        self._files_to_remove.append(filtered_otu_table_fp)
        
    def test_filter_table_default(self):
        """ filter_table removes low sample count OTUs on default settings """
        
        params = {'min_otu_count': 1, 'min_otu_samples': 2,
         'included_taxa': '', 'excluded_taxa': ''}
         
        otu_file = open(self.sample_input_otu_table, "U")
        
        filtered_otu_table_fp = get_tmp_filename(prefix = "filtered_otu_table_",
          suffix = ".txt")
        
        filtered_otu_table_f =\
         open(filtered_otu_table_fp, "w")
        
        
        filter_table(params, filtered_otu_table_f, otu_file)
        
        filtered_otu_table_f.close()
        
        # Output has 3 of the OTUs removed, that only show up in 1 sample
        
        actual_result_f = open(filtered_otu_table_fp, "U")
        
        
        actual_results = "\n".join([line.strip() for line in actual_result_f])
        
        self.assertEqual(actual_results, self.expected_otu_table_output_default)
        
        self._files_to_remove.append(filtered_otu_table_fp)
        
    def test_filter_table_sequence_count(self):
        """ filter_table removes low sequence count OTUs on default settings """
        
        params = {'min_otu_count': 3, 'min_otu_samples': 1,
         'included_taxa': '', 'excluded_taxa': ''}
         
        otu_file = open(self.sample_input_otu_table, "U")
        
        filtered_otu_table_fp = get_tmp_filename(prefix = "filtered_otu_table_",
          suffix = ".txt")
        
        filtered_otu_table_f =\
         open(filtered_otu_table_fp, "w")
        
        
        filter_table(params, filtered_otu_table_f, otu_file)
        
        filtered_otu_table_f.close()
        
        # Output has 3 of the OTUs removed, that only show up in 1 sample
        
        actual_result_f = open(filtered_otu_table_fp, "U")
        
        
        actual_results = "\n".join([line.strip() for line in actual_result_f])
        
        self.assertEqual(actual_results,
         self.expected_otu_table_output_sequence_count_three)
        
        self._files_to_remove.append(filtered_otu_table_fp)
        
        
    def test_filter_table_include_taxa(self):
        """ filter_table removes OTUs lacking targetted taxa """
        
        params = {'min_otu_count': 0, 'min_otu_samples': 0,
         'included_taxa': set(['Clostridiales']), 'excluded_taxa': ''}
         
        otu_file = open(self.sample_input_otu_table, "U")
        
        filtered_otu_table_fp = get_tmp_filename(prefix = "filtered_otu_table_",
          suffix = ".txt")
        
        filtered_otu_table_f =\
         open(filtered_otu_table_fp, "w")
        
        
        filter_table(params, filtered_otu_table_f, otu_file)
        
        filtered_otu_table_f.close()
        
        # Output has 3 of the OTUs removed, that only show up in 1 sample
        
        actual_result_f = open(filtered_otu_table_fp, "U")
        
        
        actual_results = "\n".join([line.strip() for line in actual_result_f])
        
        self.assertEqual(actual_results,
         self.expected_otu_table_include_taxa)
        
        self._files_to_remove.append(filtered_otu_table_fp)
        
        
    def test_filter_table_exclude_taxa(self):
        """ filter_table removes OTUs labeled for exclusion """
        
        params = {'min_otu_count': 0, 'min_otu_samples': 0,
         'included_taxa': '', 'excluded_taxa': set(['Bacteroidetes'])}
         
        otu_file = open(self.sample_input_otu_table, "U")
        
        filtered_otu_table_fp = get_tmp_filename(prefix = "filtered_otu_table_",
          suffix = ".txt")
        
        filtered_otu_table_f =\
         open(filtered_otu_table_fp, "w")
        
        
        filter_table(params, filtered_otu_table_f, otu_file)
        
        filtered_otu_table_f.close()
        
        # Output has 3 of the OTUs removed, that only show up in 1 sample
        
        actual_result_f = open(filtered_otu_table_fp, "U")
        
        
        actual_results = "\n".join([line.strip() for line in actual_result_f])
        
        self.assertEqual(actual_results,
         self.expected_otu_table_exclude_taxa)
        
        self._files_to_remove.append(filtered_otu_table_fp)
        
    def test_filter_table_include_exclude_taxa(self):
        """ filter_table handles both include and exclude taxa filters """
        
        params = {'min_otu_count': 0, 'min_otu_samples': 0,
         'included_taxa': set(['Firmicutes']),
         'excluded_taxa': set(['Bacillales'])}
         
        otu_file = open(self.sample_input_otu_table, "U")
        
        filtered_otu_table_fp = get_tmp_filename(prefix = "filtered_otu_table_",
          suffix = ".txt")
        
        filtered_otu_table_f =\
         open(filtered_otu_table_fp, "w")
        
        
        filter_table(params, filtered_otu_table_f, otu_file)
        
        filtered_otu_table_f.close()
        
        # Output has 3 of the OTUs removed, that only show up in 1 sample
        
        actual_result_f = open(filtered_otu_table_fp, "U")
        
        
        actual_results = "\n".join([line.strip() for line in actual_result_f])
        
        self.assertEqual(actual_results,
         self.expected_otu_table_include_exclude_taxa)
        
        self._files_to_remove.append(filtered_otu_table_fp)
        
    def test_filter_table_no_taxa(self):
        """ filter_table handles input OTU table with no taxa """
        
        params = {'min_otu_count': 1, 'min_otu_samples': 2,
         'included_taxa': '', 'excluded_taxa': ''}
         
        otu_file = open(self.sample_input_otu_table_no_taxa, "U")
        
        filtered_otu_table_fp = get_tmp_filename(prefix = "filtered_otu_table_",
          suffix = ".txt")
        
        filtered_otu_table_f =\
         open(filtered_otu_table_fp, "w")
        
        
        filter_table(params, filtered_otu_table_f, otu_file)
        
        filtered_otu_table_f.close()
        
        # Output has 3 of the OTUs removed, that only show up in 1 sample
        
        actual_result_f = open(filtered_otu_table_fp, "U")
        
        
        actual_results = "\n".join([line.strip() for line in actual_result_f])
        
        self.assertEqual(actual_results,
         self.expected_otu_table_output_default_no_taxa)
        
        self._files_to_remove.append(filtered_otu_table_fp)

    def test_filter_table_neg_control(self):
        """_filter_table_neg_control removes otus found in neg control samples
        """
        otu_table = """# QIIME v%s OTU table\n#OTU ID\tsample1\tsample2\tsample3
0\t0\t2\t0
1\t1\t0\t0
2\t1\t1\t1""" % __version__
        otu_table = otu_table.split('\n')
        samples = ['sample2', 'sample3']
        result = _filter_table_neg_control(otu_table, samples)
        self.assertEqual(result, """# QIIME v1.3.0 OTU table
#OTU ID\tsample1
1\t1""")
        #works with lineages
        otu_table = """# QIIME v%s OTU table\n#OTU ID\tsample1\tsample2\tsample3\tConsensus Lineage
0\t0\t2\t0\ttaxon1
1\t1\t0\t0\ttaxon2
2\t1\t1\t1\ttaxon3""" % __version__
        otu_table = otu_table.split('\n')
        samples = ['sample2', 'sample3']
        result = _filter_table_neg_control(otu_table, samples)
        self.assertEqual(result, """# QIIME v1.3.0 OTU table
#OTU ID\tsample1\tConsensus Lineage
1\t1\ttaxon2""")
        samples = ['sample3']
        result = _filter_table_neg_control(otu_table, samples)
        self.assertEqual(result, """# QIIME v1.3.0 OTU table
#OTU ID\tsample1\tsample2\tConsensus Lineage
0\t0\t2\ttaxon1
1\t1\t0\ttaxon2""")

# Large strings at the end for better readability
sample_unfiltered_otu_table = """#Full OTU Counts
#OTU ID	PC.354	PC.355	PC.356	PC.481	PC.593	PC.607	PC.634	PC.635	PC.636	Consensus Lineage
0	0	0	0	0	0	0	0	0	1	Root;Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae
1	0	0	0	3	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
5	0	0	0	0	1	0	0	2	1	Root;Bacteria;Bacteroidetes
6	0	0	2	0	0	0	0	1	0	Root;Bacteria
12	0	1	0	0	0	3	1	1	0	Root;Bacteria;Bacteroidetes
13	0	0	0	1	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
14	0	0	0	0	0	0	0	0	2	Root;Bacteria;Firmicutes;"Bacilli";Bacillales;"Staphylococcaceae";Staphylococcus"""

sample_unfiltered_otu_table_no_taxa = """#Full OTU Counts
#OTU ID	PC.354	PC.355	PC.356	PC.481	PC.593	PC.607	PC.634	PC.635	PC.636
0	0	0	0	0	0	0	0	0	1
1	0	0	0	3	0	0	0	0	0
5	0	0	0	0	1	0	0	2	1
6	0	0	2	0	0	0	0	1	0
12	0	1	0	0	0	3	1	1	0
13	0	0	0	1	0	0	0	1	0
14	0	0	0	0	0	0	0	0	2"""


# identical to input file, except for the header
expected_otu_table_output_leniant = """# QIIME v1.3.0 OTU table
#OTU ID	PC.354	PC.355	PC.356	PC.481	PC.593	PC.607	PC.634	PC.635	PC.636	Consensus Lineage
0	0	0	0	0	0	0	0	0	1	Root;Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae
1	0	0	0	3	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
5	0	0	0	0	1	0	0	2	1	Root;Bacteria;Bacteroidetes
6	0	0	2	0	0	0	0	1	0	Root;Bacteria
12	0	1	0	0	0	3	1	1	0	Root;Bacteria;Bacteroidetes
13	0	0	0	1	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
14	0	0	0	0	0	0	0	0	2	Root;Bacteria;Firmicutes;"Bacilli";Bacillales;"Staphylococcaceae";Staphylococcus"""

# OTUs that only appear in one sample should be removed
expected_otu_table_output_default = """# QIIME v1.3.0 OTU table
#OTU ID	PC.354	PC.355	PC.356	PC.481	PC.593	PC.607	PC.634	PC.635	PC.636	Consensus Lineage
5	0	0	0	0	1	0	0	2	1	Root;Bacteria;Bacteroidetes
6	0	0	2	0	0	0	0	1	0	Root;Bacteria
12	0	1	0	0	0	3	1	1	0	Root;Bacteria;Bacteroidetes
13	0	0	0	1	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales"""

# Remove all OTUs with sequences less than 3
expected_otu_table_output_sequence_count_three = """# QIIME v1.3.0 OTU table
#OTU ID	PC.354	PC.355	PC.356	PC.481	PC.593	PC.607	PC.634	PC.635	PC.636	Consensus Lineage
1	0	0	0	3	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
5	0	0	0	0	1	0	0	2	1	Root;Bacteria;Bacteroidetes
6	0	0	2	0	0	0	0	1	0	Root;Bacteria
12	0	1	0	0	0	3	1	1	0	Root;Bacteria;Bacteroidetes"""

# Removes all OTUs lacking the specified taxa Clostridiales
expected_otu_table_include_taxa = """# QIIME v1.3.0 OTU table
#OTU ID	PC.354	PC.355	PC.356	PC.481	PC.593	PC.607	PC.634	PC.635	PC.636	Consensus Lineage
1	0	0	0	3	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
13	0	0	0	1	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales"""

# Removes otus flagged for removal by excluded taxa, Bacteroidetes
expected_otu_table_exclude_taxa = """# QIIME v1.3.0 OTU table
#OTU ID	PC.354	PC.355	PC.356	PC.481	PC.593	PC.607	PC.634	PC.635	PC.636	Consensus Lineage
0	0	0	0	0	0	0	0	0	1	Root;Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae
1	0	0	0	3	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
6	0	0	2	0	0	0	0	1	0	Root;Bacteria
13	0	0	0	1	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
14	0	0	0	0	0	0	0	0	2	Root;Bacteria;Firmicutes;"Bacilli";Bacillales;"Staphylococcaceae";Staphylococcus"""

# Filter to include Firmicutes but exclude Bacillales
expected_otu_table_include_exclude_taxa = """# QIIME v1.3.0 OTU table
#OTU ID	PC.354	PC.355	PC.356	PC.481	PC.593	PC.607	PC.634	PC.635	PC.636	Consensus Lineage
1	0	0	0	3	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
13	0	0	0	1	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales"""


expected_otu_table_output_default_no_taxa = """# QIIME v1.3.0 OTU table
#OTU ID	PC.354	PC.355	PC.356	PC.481	PC.593	PC.607	PC.634	PC.635	PC.636
5	0	0	0	0	1	0	0	2	1
6	0	0	2	0	0	0	0	1	0
12	0	1	0	0	0	3	1	1	0
13	0	0	0	1	0	0	0	1	0"""

#run tests if called from command line
if __name__ == "__main__":
    main()
