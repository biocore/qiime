#!/usr/bin/env python

__author__ = "David Soergel"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["David Soergel"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "David Soergel"
__email__ = "soergel@cs.umass.edu"
__status__ = "Development"

from cogent.util.misc import remove_files

from unittest import TestCase, main
from qiime.pycogent_backports.rtax import (Rtax, assign_taxonomy)
from cogent.app.util import get_tmp_filename


class RtaxClassifierTests(TestCase):
    """ Tests of the RTAX classifier module """

    def setUp(self):
        self.maxDiff = None

        self.id_to_taxonomy_fp = get_tmp_filename(\
         prefix='RtaxTaxonAssignerTests_',suffix='.txt')
        self.input_seqs_fp = get_tmp_filename(\
         prefix='RtaxTaxonAssignerTests_',suffix='.fasta')
        self.reference_seqs_fp = get_tmp_filename(\
         prefix='RtaxTaxonAssignerTests_',suffix='.fasta')
        self.read_1_seqs_fp = get_tmp_filename(\
         prefix='RtaxTaxonAssignerTests_',suffix='.fasta')
        self.read_2_seqs_fp = get_tmp_filename(\
         prefix='RtaxTaxonAssignerTests_',suffix='.fasta')

        self._paths_to_clean_up = [self.id_to_taxonomy_fp,self.input_seqs_fp,self.reference_seqs_fp, self.read_1_seqs_fp,self.read_2_seqs_fp]

        a = open(self.id_to_taxonomy_fp,'w')
        a.write(rtax_reference_taxonomy)
        a.close()
        b = open(self.reference_seqs_fp,'w')
        b.write(rtax_reference_fasta)
        b.close()
        c = open(self.input_seqs_fp,'w')
        c.write(rtax_test_repset_fasta)
        c.close()
        d = open(self.read_1_seqs_fp,'w')
        d.write(rtax_test_read1_fasta)
        d.close()
        e = open(self.read_2_seqs_fp,'w')
        e.write(rtax_test_read2_fasta)
        e.close()

    def tearDown(self):
        remove_files(set(self._paths_to_clean_up),error_on_missing=False)

    def test_paired_end_classification(self):
        self._paths_to_clean_up += cleanAll(self.read_1_seqs_fp)
        self._paths_to_clean_up += cleanAll(self.read_2_seqs_fp)
        result = assign_taxonomy(self.input_seqs_fp, self.reference_seqs_fp, self.id_to_taxonomy_fp, self.read_1_seqs_fp, self.read_2_seqs_fp,single_ok=False,header_id_regex="\\S+\\s+(\\S+?)\/")
        self.assertEqual(result, rtax_expected_result_paired)

    def test_paired_end_classification_with_fallback(self):
        self._paths_to_clean_up += cleanAll(self.read_1_seqs_fp)
        self._paths_to_clean_up += cleanAll(self.read_2_seqs_fp)
        result = assign_taxonomy(self.input_seqs_fp, self.reference_seqs_fp, self.id_to_taxonomy_fp, self.read_1_seqs_fp, self.read_2_seqs_fp,single_ok=True,header_id_regex="\\S+\\s+(\\S+?)\/")
        self.assertEqual(result, rtax_expected_result_paired_with_fallback)

    def test_single_end_classification(self):
        self._paths_to_clean_up += cleanAll(self.read_1_seqs_fp)
        result = assign_taxonomy(self.input_seqs_fp, self.reference_seqs_fp, self.id_to_taxonomy_fp, self.read_1_seqs_fp, None ,header_id_regex="\\S+\\s+(\\S+?)\/")
        self.assertEqual(result, rtax_expected_result_single)

    # I'd like to add tests here that involve the TOOMANYHITS case.  However, that requires either a reference
    # database with >16,000 sequences, which we don't have handy for tests, or adjusting the maxMaxAccepts parameter to rtaxSearch.pl.
    # However the "rtax" wrapper shell script currently doesn't allow setting that option, and I'd prefer to leave that as is
    # unless someone actually wants to use it.  Thus the TOOMANYHITS situation is not easily testable at the moment.


def cleanAll(path):
    return [path, path + ".pos.db",  path + ".pos.dir", path + ".pos.pag", path + ".lines.db", path + ".lines.dir", path + ".lines.pag"]


# sample data copied from GreenGenes


rtax_reference_taxonomy = """508720	99.0	k__Bacteria	 p__Actinobacteria	 c__Actinobacteria	 o__Actinomycetales	 f__Propionibacteriaceae	 g__Propionibacterium	 s__Propionibacterium acnes
508050	99.0	k__Bacteria	 p__Proteobacteria	 c__Betaproteobacteria	 o__Burkholderiales	 f__Comamonadaceae	 g__Diaphorobacter	 s__
502492	99.0	k__Bacteria	 p__Proteobacteria	 c__Betaproteobacteria	 o__Burkholderiales	 f__	 g__Aquabacterium	 s__
"""

rtax_reference_fasta = """>508720
GACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACGGAAAGGCCCTGCTTTTGTGGGGTGCTCGAGTGGCGAACG
GGTGAGTAACACGTGAGTAACCTGCCCTTGACTTTGGGATAACTTCAGGAAACTGGGGCTAATACCGGATAGGAGCTCCT
GCTGCATGGTGGGGGTTGGAAAGTTTCGGCGGTTGGGGATGGACTCGCGGCTTATCAGCTTGTTGGTGGGGTAGTGGCTT
ACCAAGGCTTTGACGGGTAGCCGGCCTGAGAGGGTGACCGGCCACATTGGGACTGAGATACGGCCCAGACTCCTACGGGA
GGCAGCAGTGGGGAATATTGCACAATGGGCGGAAGCCTGATGCAGCAACGCCGCGTGCGGGATGACGGCCTTCGGGTTGT
AAACCGCTTTCGCCTGTGACGAAGCGTGAGTGACGGTAATGGGTAAAGAAGCACCGGCTAACTACGTGCCAGCAGCCGCG
GTGATACGTAGGGTGCGAGCGTTGTCCGGATTTATTGGGCGTAAAGGGCTCGTAGGTGGTTGATCGCGTCGGAAGTGTAA
TCTTGGGGCTTAACCCTGAGCGTGCTTTCGATACGGGTTGACTTGAGGAAGGTAGGGGAGAATGGAATTCCTGGTGGAGC
GGTGGAATGCGCAGATATCAGGAGGAACACCAGTGGCGAAGGCGGTTCTCTGGGCCTTTCCTGACGCTGAGGAGCGAAAG
CGTGGGGAGCGAACAGGCTTAGATACCCTGGTAGTCCACGCTGTAAACGGTGGGTACTAGGTGTGGGGTCCATTCCACGG
GTTCCGTGCCGTAGCTAACGCTTTAAGTACCCCGCCTGGGGAGTACGGCCGCAAGGCTAAAACTCAAAGGAATTGACGGG
GCCCCGCACAAGCGGCGGAGCATGCGGATTAATTCGATGCAACGCGTAGAACCTTACCTGGGTTTGACATGGATCGGGAG
TGCTCAGAGATGGGTGTGCCTCTTTTGGGGTCGGTTCACAGGTGGTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTT
GGGTTAAGTCCCGCAACGAGCGCAACCCTTGTTCACTGTTGCCAGCACGTTATGGTGGGGACTCAGTGGAGACCGCCGGG
GTCAACTCGGAGGAAGGTGGGGATGACGTCAAGTCATCATGCCCCTTATGTCCAGGGCTTCACGCATGCTACAATGGCTG
GTACAGAGAGTGGCGAGCCTGTGAGGGTGAGCGAATCTCGGAAAGCCGGTCTCAGTTCGGATTGGGGTCTGCAACTCGAC
CTCATGAAGTCGGAGTCGCTAGTAATCGCAGATCAGCAACGCTGCGGTGAATACGTTCCCGGGGCT
>508050
ATTGAACGCTGGCGGCATGCCTTACACATGCAAGTCGAACGGTAACAGGTCTTCGGATGCTGACGAGTGGCGAACGGGTG
AGTAATACATCGGAACGTGCCCGATCGTGGGGGATAACGAGGCGAAAGCTTTGCTAATACCGCATACGATCTACGGATGA
AAGCGGGGGATCTTCGGACCTCGCGCGGACGGAGCGGCCGATGGCAGATTAGGTAGTTGGTGGGATAAAAGCTTACCAAG
CCGACGATCTGTAGCTGGTCTGAGAGGATGATCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGC
AGTGGGGAATTTTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGCAGGATGAAGGCCTTCGGGTTGTAAACTG
CTTTTGTACGGAACGAAAAGCCTCTTTCTAATAAAGAGGGGTCATGACGGTACCGTAAGAATAAGCACCGGCTAACTACG
TGCCAGCAGCCGCGGTAATACGTAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGTTTTGTA
AGACAGAGGTGAAATCCCCGGGCTCAACCTGGGAACTGCCTTTGTGACTGCAAGGCTGGAGTGCGGCAGAGGGGGATGGA
ATTCCGCGTGTAGCAGTGAAATGCGTAGATATGCGGAGGAACACCGATGGCGAAGGCAATCCCCTGGGCCTGCACTGACG
CTCATGCACGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCCTAAACGATGTCAACTGGTTGTTG
GGTCTTCACTGACTCAGTAACGAAGCTAACGCGTGAAGTTGACCGCCTGGGGAGTACGGCCGCAAGGTTGAAACTCAAAG
GAATTGACGGGGACCCGCACAAGCGGTGGATGATGTGGTTTAATTCGATGCAACGCGAAAAACCTTACCCACCTTTGACA
TGGCAGGAAGTTTCCAGAGATGGATTCGTGCCCGAAAGGGAACCTGCACACAGGTGCTGCATGGCTGTCGTCAGCTCGTG
TCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTGCCATTAGTTGCTACGAAAGGGCACTCTAATGGGACTG
CCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAAGTCCTCATGGCCCTTATAGGTGGGGCTACACACGTCATACAAT
GGCTGGTACAGAGGGTTGCCAACCCGCGAGGGGGAGCTAATCCCATAAAGCCAGTCGTAGTCCGGATCGCAGTCTGCAAC
TCGACTGCGTGAAGTCGGAATCGCTAGTAATCGCGGATCAGAATGTCGCGGTGAATACGTTCCCGGGTCT
>502492
ATTGAACGCTGGCGGCATGCCTTACACATGCAAGTCGAACGGTAACGGGTCCTTCGGGATGCCGACGAGTGGCGAACGGG
TGAGTAATATATCGGAACGTGCCCAGTAGTGGGGGATAACTGCTCGAAAGAGCAGCTAATACCGCATACGACCTGAGGGT
GAAAGGGGGGGATCGCAAGACCTCTCGCTATTGGAGCGGCCGATATCAGATTAGCTAGTTGGTGGGGTAAAGGCCTACCA
AGGCAACGATCTGTAGTTGGTCTGAGAGGACGACCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCA
GCAGTGGGGAATTTTGGACAATGGGCGCAAGCCTGATCCAGCAATGCCGCGTGCAGGAAGAAGGCCTTCGGGTTGTAAAC
TGCTTTTGTCAGGGAAGAAATCTTCTGGGCTAATACCCCGGGAGGATGACGGTACCTGAAGAATAAGCACCGGCTAACTA
CGTGCCAGCAGCCGCGGTAATACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGCTTTG
CAAGACAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTTGTGACTGCAAGGCTAGAGTACGGCAGAGGGGGATG
GAATTCCGCGTGTAGCAGTGAAATGCGTAGATATGCGGAGGAACACCAATGGCGAAGGCAATCCCCTGGGCCTGTACTGA
CGCTCATGCACGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCCTAAACGATGTCAACTGGTTGT
TGGACGGCTTGCTGTTCAGTAACGAAGCTAACGCGTGAAGTTGACCGCCTGGGGAGTACGGCCGCAAGGTTGAAACTCAA
AGGAATTGACGGGGACCCGCACAAGCGGTGGATGATGTGGTTTAATTCGATGCAACGCGAAAAACCTTACCTACCCTTGA
CATGTCAAGAATTCTGCAGAGATGTGGAAGTGCTCGAAAGAGAACTTGAACACAGGTGCTGCATGGCCGTCGTCAGCTCG
TGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTGTCATTAGTTGCTACGCAAGAGCACTCTAATGAGAC
TGCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAGGTCCTCATGGCCCTTATGGGTAGGGCTACACACGTCATACA
ATGGCCGGTACAGAGGGCTGCCAACCCGCGAGGGGGAGCCAATCCCAGAAAACCGGTCGTAGTCCGGATCGTAGTCTGCA
ACTCGACTGCGTGAAGTCGGAATCGCTAGTAATCGCGGATCAGCTTGCCGCGGTGAATACGTTCCCGGGTCT
"""


rtax_test_repset_fasta = """>clusterIdA splitRead1IdA
ACCAAGGCTTTGACGGGTAGCCGGCCTGAGTGGGTGACCGGCCACATTGGGACTGAGATACGGCCCAGACTCCTACGGGA
>clusterIdB splitRead1IdB
CCGACGATCTGTAGCTGGTCTGAGAGGATGTTCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGC
>clusterIdC splitRead1IdC
AGGCAACGATCTGTAGTTGGTCTGAGAGGAGGACCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCA
>clusterIdD splitRead1IdD
AGGCAACGATCTGTAGTTGGTCTGAGAGGAGGACCAGCCACACTGGGACGGGGGGGGGGCCCAGACTCCTACGGGAGGCA
"""

# these reads are the 4th and 14th lines from the reference seqs

#rtax_test_read1_fasta = """>splitRead1IdA ampliconId_34563456/1
#ACCAAGGCTTTGACGGGTAGCCGGCCTGAGAGGGTGACCGGCCACATTGGGACTGAGATACGGCCCAGACTCCTACGGGA
#>splitRead1IdB ampliconId_
#CCGACGATCTGTAGCTGGTCTGAGAGGATGATCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGC
#>splitRead1IdC ampliconId_
#AGGCAACGATCTGTAGTTGGTCTGAGAGGACGACCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCA
#"""
#
#rtax_test_read2_fasta = """>splitRead2IdA ampliconId_34563456/3
#GGGTTAAGTCCCGCAACGAGCGCAACCCTTGTTCACTGTTGCCAGCACGTTATGGTGGGGACTCAGTGGAGACCGCCGGG
#>splitRead2IdB ampliconId_
#TCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTGCCATTAGTTGCTACGAAAGGGCACTCTAATGGGACTG
#>splitRead2IdC ampliconId_
#TGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTGTCATTAGTTGCTACGCAAGAGCACTCTAATGAGAC
#"""


# these reads are the 4th and 14th lines from the reference seqs, with one nucleotide changed each
# except D and E, which are unique to one read or the other
# and F and G, which are just decoys

rtax_test_read1_fasta = """>splitRead1IdA ampliconId_34563456/1
ACCAAGGCTTTGACGGGTAGCCGGCCTGAGTGGGTGACCGGCCACATTGGGACTGAGATACGGCCCAGACTCCTACGGGA
>splitRead1IdB ampliconId_12341234/1
CCGACGATCTGTAGCTGGTCTGAGAGGATGTTCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGC
>splitRead1IdC ampliconId_23452345/1
AGGCAACGATCTGTAGTTGGTCTGAGAGGAGGACCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCA
>splitRead1IdD ampliconId_45674567/1
AGGCAACGATCTGTAGTTGGTCTGAGAGGAGGACCAAAAAAAAAAAGACTGAGACACGGCCCAGACTCCTACGGGAGGCA
>splitRead1IdF ampliconId_56785678/1
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
"""

rtax_test_read2_fasta = """>splitRead2IdA ampliconId_34563456/3
GGGTTAAGTCCCGCAACGAGCGCAACCCTTATTCACTGTTGCCAGCACGTTATGGTGGGGACTCAGTGGAGACCGCCGGG
>splitRead2IdB ampliconId_12341234/3
TCGTGAGATGTTGGGTTAAGTCCCGCAACGTGCGCAACCCTTGCCATTAGTTGCTACGAAAGGGCACTCTAATGGGACTG
>splitRead2IdC ampliconId_23452345/3
TGTCGTGAGATGTTGGGTTAAGTCCCGCAAAGAGCGCAACCCTTGTCATTAGTTGCTACGCAAGAGCACTCTAATGAGAC
>splitRead2IdE ampliconId_67896789/3
TGTCGTGAGATGTTGGGTTAAAAAAAAAAAAAAACGCAACCCTTGTCATTAGTTGCTACGCAAGAGCACTCTAATGAGAC
>splitRead2IdG ampliconId_78907890/3
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
"""


rtax_expected_result_paired = {
    'clusterIdA splitRead1IdA': ('k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales; f__Propionibacteriaceae; g__Propionibacterium; s__Propionibacterium acnes', 1.0),
    'clusterIdB splitRead1IdB': ('k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Comamonadaceae; g__Diaphorobacter; s__', 1.0),
    'clusterIdC splitRead1IdC': ('k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__; g__Aquabacterium; s__', 1.0),
    }

rtax_expected_result_paired_with_fallback = {
    'clusterIdA splitRead1IdA': ('k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales; f__Propionibacteriaceae; g__Propionibacterium; s__Propionibacterium acnes', 1.0),
    'clusterIdB splitRead1IdB': ('k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Comamonadaceae; g__Diaphorobacter; s__', 1.0),
    'clusterIdC splitRead1IdC': ('k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__; g__Aquabacterium; s__', 1.0),
    'clusterIdD splitRead1IdD': ('k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__; g__Aquabacterium; s__', 1.0),
    }

rtax_expected_result_single = {
    'clusterIdA splitRead1IdA': ('k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales; f__Propionibacteriaceae; g__Propionibacterium; s__Propionibacterium acnes', 1.0),
    'clusterIdB splitRead1IdB': ('k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Comamonadaceae; g__Diaphorobacter; s__', 1.0),
    'clusterIdC splitRead1IdC': ('k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__; g__Aquabacterium; s__', 1.0),
    'clusterIdD splitRead1IdD': ('k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__; g__Aquabacterium; s__', 1.0),
    }

if __name__ == "__main__":
    main()


