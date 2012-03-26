#!/usr/bin/env python
"""Application controller for RAxML (v7.3.0).

WARNING: Because of the use of the -x option, this version is no longer
compatible with RAxML version VI.
"""
from cogent.app.parameters import FlagParameter, ValuedParameter, FilePath
from cogent.app.util import CommandLineApplication, ResultPath, \
                            get_tmp_filename,ApplicationError
from cogent.core.tree import PhyloNode
from cogent.core.alignment import Alignment
from cogent.core.moltype import DNA, RNA, PROTEIN
from random import choice, randint
from os import walk,listdir
from os.path import isabs,join,split
from cogent.parse.tree import DndParser
import re
from cogent.app.guppy import build_tree_from_json_using_params

__author__ = "Micah Hamady"
__copyright__ = "Copyright 2007-2011, The Cogent Project"
__credits__ = ["Micah Hamady", "Catherine Lozupone", "Rob Knight", \
               "Daniel McDonald","Jesse Stombaugh"]
__license__ = "GPL"
__version__ = "1.6.0dev"
__maintainer__ = "Micah Hamady"
__email__ = "hamady@colorado.edu"
__status__ = "Prototype"

class Raxml(CommandLineApplication):
    """RAxML application controller"""

    _options ={

        # Specify a column weight file name to assign individual wieghts to 
        # each column of the alignment. Those weights must be integers 
        # separated by any number and type of whitespaces whithin a separate 
        # file, see file "example_weights" for an example.
        '-a':ValuedParameter('-',Name='a',Delimiter=' '),

        # Specify one of the secondary structure substitution models implemented
        # in RAxML. The same nomenclature as in the PHASE manual is used, 
        # available models:  S6A, S6B, S6C, S6D, S6E, S7A, S7B, S7C, S7D, S7E, 
        # S7F, S16, S16A, S16B
        # DEFAULT: 16-state GTR model (S16)
        '-A':ValuedParameter('-',Name='A',Delimiter=' '),
        
        #  Specify an integer number (random seed) for bootstrapping
        '-b':ValuedParameter('-',Name='b',Delimiter=' '),
        
        # specify a floating point number between 0.0 and 1.0 that will be used 
        # as cutoff threshold for the MR-based bootstopping criteria. The 
        # recommended setting is 0.03.
        '-B':ValuedParameter('-',Name='B',Delimiter=' '),
        
        # Specify number of distinct rate catgories for raxml when 
        # ModelOfEvolution is set to GTRCAT or HKY85CAT.
        # Individual per-site rates are categorized into numberOfCategories 
        # rate categories to accelerate computations. (Default = 50)
        '-c':ValuedParameter('-',Name='c',Delimiter=' '),

        # Conduct model parameter optimization on gappy, partitioned multi-gene 
        # alignments with per-partition branch length estimates (-M enabled) 
        # using the fast method with pointer meshes described in:
        # Stamatakis and Ott: "Efficient computation of the phylogenetic 
        # likelihood function on multi-gene alignments and multi-core 
        # processors"
        # WARNING: We can not conduct useful tree searches using this method 
        # yet! Does not work with Pthreads version.
        '-C':ValuedParameter('-',Name='C',Delimiter=' '),          

        # This option allows you to start the RAxML search with a complete 
        # random starting tree instead of the default Maximum Parsimony 
        # Starting tree. On smaller datasets (around 100-200 taxa) it has 
        # been observed that this might sometimes yield topologies of distinct 
        # local likelihood maxima which better correspond to empirical 
        # expectations. 
        '-d':FlagParameter('-',Name='d'),

        # ML search convergence criterion. This will break off ML searches if 
        # the relative Robinson-Foulds distance between the trees obtained from 
        # two consecutive lazy SPR cycles is smaller or equal to 1%. Usage 
        # recommended for very large datasets in terms of taxa. On trees with 
        # more than 500 taxa this will yield execution time improvements of 
        # approximately 50% While yielding only slightly worse trees.
        # DEFAULT: OFF
        '-D':ValuedParameter('-',Name='D'), 

        # This allows you to specify up to which likelihood difference.
        # Default is 0.1 log likelihood units, author recommends 1 or 2 to
        # rapidly evaluate different trees.
        '-e':ValuedParameter('-',Name='e',Delimiter=' '),
        
        # specify an exclude file name, that contains a specification of 
        # alignment positions you wish to exclude. Format is similar to Nexus, 
        # the file shall contain entries like "100-200 300-400", to exclude a 
        # single column write, e.g., "100-100", if you use a mixed model, an 
        # appropriatly adapted model file will be written.
        '-E':ValuedParameter('-',Name='E',Delimiter=' '),

        # select search algorithm: 
        #   a rapid Bootstrap analysis and search for best-scoring ML tree in 
        #       one program run
        #   A compute marginal ancestral states on a ROOTED reference tree
        #       provided with "t" - ONLY IN 7.3.0
        #   b draw bipartition information on a tree provided with "-t" based on 
        #       multiple trees (e.g., from a bootstrap) in a file specifed by 
        #       "-z"
        #   c check if the alignment can be properly read by RAxML
        #   d for normal hill-climbing search (Default)
        #     when -f option is omitted this algorithm will be used
        #   e optimize model+branch lengths for given input tree under 
        #       GAMMA/GAMMAI only
        #   E execute very fast experimental tree search, at present only for 
        #       testing
        #   F execute fast experimental tree search, at present only for testing
        #   g compute per site log Likelihoods for one ore more trees passed via
        #       "-z" and write them to a file that can be read by CONSEL
        #       WARNING: does not print likelihoods in the original column order
        #   h compute log likelihood test (SH-test) between best tree passed via 
        #       "-t" and a bunch of other trees passed via "-z" 
        #   i EXPERIMENTAL do not use for real tree inferences: conducts a 
        #       single cycle of fast lazy SPR moves on a given input tree, to be 
        #       used in combination with -C and -M 
        #   I EXPERIMENTAL do not use for real tree inferences: conducts a 
        #       single cycle of thorough lazy SPR moves on a given input tree,  
        #       to be used in combination with -C and -M 
        #   j generate a bunch of bootstrapped alignment files from an original 
        #       alignemnt file. You need to specify a seed with "-b" and the 
        #       number of replicates with "-#" 
        # following "J" is for version 7.2.8
        #   J Compute SH-like support values on a given tree passed via "-t".
        #   m compare bipartitions between two bunches of trees passed via "-t" 
        #       and "-z" respectively. This will return the Pearson correlation 
        #       between all bipartitions found in the two tree files. A file 
        #       called RAxML_bipartitionFrequencies.outpuFileName will be 
        #       printed that contains the pair-wise bipartition frequencies of 
        #       the two sets
        #   n compute the log likelihood score of all trees contained in a tree 
        #       file provided by "-z" under GAMMA or GAMMA+P-Invar
        #   o old (slower) algorithm from v. 2.1.3
        #   p perform pure stepwise MP addition of new sequences to an 
        #       incomplete starting tree and exit
        #   r compute pairwise Robinson-Foulds (RF) distances between all pairs 
        #       of trees in a tree file passed via "-z" if the trees have node 
        #       labales represented as integer support values the program will 
        #       also compute two flavors of the weighted Robinson-Foulds (WRF)
        #       distance
        # following "R" is for version 7.2.8
        #   R compute rogue taxa using new statistical method based on the
        #       evolutionary placement algorithm
        #       WARNING: this is experimental code - DEPRECATED IN 7.3.0
        #   s (split) splits into individual genes, provided with model file
        # following "S" is for version 7.2.8
        #   S compute site-specific placement bias using a leave one out test
        #       inspired by the evolutionary placement algorithm
        #   t do randomized tree searches on one fixed starting tree
        #   u execute morphological weight calibration using maximum likelihood, 
        #       this will return a weight vector. you need to provide a 
        #       morphological alignment and a reference tree via "-t" 
        #   U execute morphological wieght calibration using parsimony, this 
        #       will return a weight vector. you need to provide a morphological 
        #       alignment and a reference tree via "-t" - DEPRECATED IN 7.3.0
        #   v classify a bunch of environmental sequences into a reference tree 
        #       using the slow heuristics without dynamic alignment you will 
        #       need to start RAxML with a non-comprehensive reference tree and 
        #       an alignment containing all sequences (reference + query)
        #   w compute ELW test on a bunch of trees passed via "-z" 
        #   x compute pair-wise ML distances, ML model parameters will be 
        #       estimated on an MP starting tree or a user-defined tree passed 
        #       via "-t", only allowed for GAMMA-based models of rate 
        #       heterogeneity
        #   y classify a bunch of environmental sequences into a reference tree 
        #       using the fast heuristics without dynamic alignment you will 
        #       need to start RAxML with a non-comprehensive reference tree and 
        #       an alignment containing all sequences (reference + query)
        '-f':ValuedParameter('-',Name='f',Delimiter=' ', Value="d"),

        # enable ML tree searches under CAT model for very large trees without 
        # switching to GAMMA in the end (saves memory). This option can also be 
        # used with the GAMMA models in order to avoid the thorough optimization 
        # of the best-scoring ML tree in the end.
        # DEFAULT: OFF
        '-F':FlagParameter('-',Name='F'),
        
        # select grouping file name: allows incomplete multifurcating constraint
        # tree in newick format -- resolves multifurcations randomly, adds
        # other taxa using parsimony insertion
        '-g':ValuedParameter('-', Name='g',Delimiter=' '),

        # enable the ML-based evolutionary placement algorithm heuristics by 
        # specifiyng a threshold value (fraction of insertion branches to be 
        # evaluated using slow insertions under ML).
        '-G':FlagParameter('-', Name='G'),

        # prints help and exits
        '-h':FlagParameter('-', Name='h'),

        # enable the MP-based evolutionary placement algorithm heuristics
        # by specifiyng a threshold value (fraction of insertion branches to be 
        # evaluated using slow insertions under ML) - DEPRECATED IN 7.3.0
        #'-H':ValuedParameter('-', Name='H',Delimiter=' '),
        
        # allows initial rearrangement to be constrained, e.g. 10 means
        # insertion will not be more than 10 nodes away from original.
        # default is to pick a "good" setting.
        '-i':ValuedParameter('-', Name='i', Delimiter=' '),

        # a posteriori bootstopping analysis. Use:
        #   "-I autoFC" for the frequency-based criterion
        #   "-I autoMR" for the majority-rule consensus tree criterion
        #   "-I autoMRE" for the extended majority-rule consensus tree criterion
        #   "-I autoMRE_IGN" for metrics similar to MRE, but include 
        #       bipartitions under the threshold whether they are compatible
        #       or not. This emulates MRE but is faster to compute.
        #   You also need to pass a tree file containg several bootstrap 
        #   replicates via "-z"
        '-I':ValuedParameter('-', Name='I', Delimiter=' '),
        
        # writes checkpoints (off by default)
        '-j':FlagParameter('-', Name='j'),

        # Compute majority rule consensus tree with "-J MR" or extended majority 
        # rule consensus tree with "-J MRE" or strict consensus tree with "-J 
        # STRICT" You will need to provide a tree file containing several 
        # UNROOTED trees via "-z"
        '-J':ValuedParameter('-', Name='J', Delimiter=' '),
        
        #specifies that RAxML will optimize model parameters (for GTRMIX and
        # GTRGAMMA) as well as calculating likelihoods for bootstrapped trees.
        '-k':FlagParameter('-', Name='k'),

        # Specify one of the multi-state substitution models (max 32 states) 
        # implemented in RAxML. Available models are: ORDERED, MK, GTR
        '-K':ValuedParameter('-', Name='K', Delimiter=' '),
        
        # Model of Binary (Morphological), Nucleotide, Multi-State, or Amino 
        #   Acid Substitution::
        # BINARY:
        #   -m BINCAT : Optimization of site-specific evolutionary rates which 
        #       are categorized into numberOfCategories distinct rate categories 
        #       for greater computational efficiency. Final tree might be 
        #       evaluated automatically under BINGAMMA, depending on the tree 
        #       search option
        #   -m BINCATI : Optimization of site-specific evolutionary rates which 
        #       are categorized into numberOfCategories distinct rate categories    
        #       for greater computational efficiency. Final tree might be 
        #       evaluated automatically under BINGAMMAI, depending on the tree 
        #       search option 
        #   -m BINGAMMA : GAMMA model of rate heterogeneity (alpha parameter 
        #       will be estimated)
        #   -m BINGAMMAI : Same as BINGAMMA, but with estimate of proportion of 
        #       invariable sites
        # NUCLEOTIDES
        #   -m GTRCAT: GTR + Optimization of substitution rates +  Optimization 
        #       of site-specific evolutionary rates which are categorized into 
        #       numberOfCategories distinct rate categories for greater 
        #       computational efficiency
        #   -m GTRCAT_FLOAT : Same as above but uses single-precision floating 
        #       point arithemtics instead of double-precision Usage only 
        #       recommened for testing, the code will run slower, but can save 
        #       almost 50% of memory. If you have problems with phylogenomic 
        #       datasets and large memory requirements you may give it a shot. 
        #       Keep in mind that numerical stability seems to be okay but needs 
        #       further testing. - DEPRECATED IN 7.3.0
        #   -m GTRCATI : GTR + Optimization of substitution rates + Optimization 
        #       of site-specific evolutionary rates which are categorized into 
        #       numberOfCategories distinct rate categories for greater 
        #       computational efficiency.  Final tree might be evaluated under 
        #       GTRGAMMAI, depending on the tree search option
        #   -m GTRGAMMA: GTR + Optimization of substitution rates + Gamma
        #   -m GTRGAMMA_FLOAT : Same as GTRGAMMA, but also with 
        #       single-precision arithmetics, same cautionary notes as for  
        #       GTRCAT_FLOAT apply. - DEPRECATED IN 7.3.0
        #   -m GTRGAMMAI : Same as GTRGAMMA, but with estimate of proportion of 
        #       invariable sites 
        # MULTI-STATE:
        #   -m MULTICAT : Optimization of site-specific evolutionary rates which 
        #       are categorized into numberOfCategories distinct rate categories 
        #       for greater computational efficiency. Final tree might be 
        #       evaluated automatically under MULTIGAMMA, depending on the tree 
        #       search option
        #   -m MULTICATI : Optimization of site-specific evolutionary rates 
        #       which are categorized into numberOfCategories distinct rate 
        #       categories for greater computational efficiency. Final tree 
        #       might be evaluated automatically under MULTIGAMMAI, depending on 
        #       the tree search option 
        #   -m MULTIGAMMA : GAMMA model of rate heterogeneity (alpha parameter 
        #       will be estimated)
        #   -m MULTIGAMMAI : Same as MULTIGAMMA, but with estimate of proportion 
        #       of invariable sites
        # You can use up to 32 distinct character states to encode multi-state
        # regions, they must be used in the following order: 0, 1, 2, 3, 4, 5, 
        # 6, 7, 8, 9, A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, 
        # T, U, V i.e., if you have 6 distinct character states you would use 0, 
        # 1, 2, 3, 4, 5 to encode these. The substitution model for the
        # multi-state regions can be selected via the "-K" option
        # Amino Acid Models:
        #   -m PROTCATmatrixName[F] : specified AA matrix + Optimization of 
        #       substitution rates + Optimization of site-specific evolutionary 
        #       rates which are categorized into numberOfCategories distinct 
        #       rate categories for greater computational efficiency.   Final
        #       tree might be evaluated automatically under 
        #       PROTGAMMAmatrixName[f], depending on the tree search option
        #   -m PROTCATmatrixName[F]_FLOAT : PROTCAT with single precision 
        #       arithmetics, same cautionary notes as for GTRCAT_FLOAT apply
        #       - DEPRECATED IN 7.3.0
        #   -m PROTCATImatrixName[F] : specified AA matrix + Optimization of 
        #       substitution rates + Optimization of site-specific
        #       evolutionary rates which are categorized into numberOfCategories 
        #       distinct rate categories for greater computational efficiency.   
        #       Final tree might be evaluated automatically under 
        #       PROTGAMMAImatrixName[f], depending on the tree search option
        #   -m PROTGAMMAmatrixName[F] : specified AA matrix + Optimization of 
        #       substitution rates + GAMMA model of rate heterogeneity (alpha 
        #       parameter will be estimated)
        #   -m PROTGAMMAmatrixName[F]_FLOAT : PROTGAMMA with single precision 
        #       arithmetics, same cautionary notes as for GTRCAT_FLOAT apply
        #       - DEPRECATED IN 7.3.0
        #   -m PROTGAMMAImatrixName[F] : Same as PROTGAMMAmatrixName[F], but 
        #       with estimate of proportion of invariable sites 
        # Available AA substitution models: DAYHOFF, DCMUT, JTT, MTREV, WAG, 
        # RTREV, CPREV, VT, BLOSUM62, MTMAM, LG, GTR. With the optional "F" 
        # appendix you can specify if you want to use empirical base frequencies
        # Please note that for mixed models you can in addition specify the 
        # per-gene AA model in the mixed model file (see manual for details). 
        # Also note that if you estimate AA GTR parameters on a partitioned
        # dataset, they will be linked (estimated jointly) across all partitions 
        # to avoid over-parametrization
        '-m':ValuedParameter('-',Name='m',Delimiter=' '),

        # Switch on estimation of individual per-partition branch lengths. Only 
        # has effect when used in combination with "-q". Branch lengths for 
        # individual partitions will be printed to separate files. A weighted 
        # average of the branch lengths is computed by using the respective 
        # partition lengths. 
        # DEFAULT: OFF
        '-M':FlagParameter('-',Name='M'),
        
        # Specifies the name of the output file.
        '-n':ValuedParameter('-',Name='n',Delimiter=' '),

        # Specifies the name of the outgroup (or outgroups: comma-delimited,
        # no spaces, should be monophyletic).
        '-o':ValuedParameter('-',Name='o',Delimiter=' '),

        # Enable checkpointing using the dmtcp library available at 
        # http://dmtcp.sourceforge.net/. This only works if you call the program 
        # by preceded by the command "dmtcp_checkpoint" and if you compile a 
        # dedicated binary using the appropriate Makefile. With "-O" you can 
        # specify the interval between checkpoints in seconds.
        # DEFAULT: 3600.0 seconds - DEPRECATED IN 7.3.0
        #'-O':ValuedParameter('-',Name='O',Delimiter=' ',Value=3600.0),

        # Specify a random number seed for the parsimony inferences. This allows 
        # you to reproduce your results and will help me debug the program.
        '-p':ValuedParameter('-',Name='p',Delimiter=' '),
        
        # Specify the file name of a user-defined AA (Protein) substitution 
        # model. This file must contain 420 entries, the first 400 being the AA 
        # substitution rates (this must be a symmetric matrix) and the last 20 
        # are the empirical base frequencies
        '-P':ValuedParameter('-',Name='P',Delimiter=' '),

        # Specified MultipleModel file name, in format:
        #    gene1 = 1-500
        #    gene2 = 501-1000
        #    (note: ranges can also be discontiguous, e.g. 1-100, 200-300,
        #     or can specify codon ranges as e.g. 1-100/3, 2-100/3, 3-100/3))
        '-q':ValuedParameter('-', Name='q', Delimiter=' '),

        # THE FOLLOWING "Q" is DEPRECATED IN 7.2.8
        # Turn on computation of SH-like support values on tree.
        # DEFAULT: OFF
        '-Q':FlagParameter('-', Name='Q'),
        
        # Constraint file name: allows a bifurcating Newick tree to be passed
        # in as a constraint file, other taxa will be added by parsimony.
        '-r':ValuedParameter('-',Name='r',Delimiter=' '),
        
        # THE FOLLOWING "R" is IN 7.2.8 
        # Specify the file name of a binary model parameter file that has
        # previously been generated with RAxML using the -f e tree evaluation
        # option. The file name should be:  RAxML_binaryModelParameters.runID
        '-R':ValuedParameter('-',Name='R',Delimiter=' '),
        
        # specify the name of the alignment data file, in relaxed PHYLIP
        # format.
        '-s':ValuedParameter('-',Name='s',Delimiter=' '),

        # Specify the name of a secondary structure file. The file can contain 
        # "." for alignment columns that do not form part of a stem and 
        # characters "()<>[]{}" to define stem regions and pseudoknots
        '-S':ValuedParameter('-',Name='S',Delimiter=' '),
        
        # Specify a user starting tree file name in Newick format
        '-t':ValuedParameter('-',Name='t',Delimiter=' '),

        # PTHREADS VERSION ONLY! Specify the number of threads you want to run.
        # Make sure to set "-T" to at most the number of CPUs you have on your 
        # machine, otherwise, there will be a huge performance decrease! 
        '-T':ValuedParameter('-',Name='T',Delimiter=' '),
        
        # THE FOLLOWING "U" is IN 7.2.8 
        # Try to save memory by using SEV-based implementation for gap columns
        # on large gappy alignments
        # WARNING: this will only work for DNA under GTRGAMMA and is still in an
        # experimental state.
        '-U':ValuedParameter('-',Name='U',Delimiter=' '),
        
        # Print the version
        '-v':FlagParameter('-',Name='v'),

        # Name of the working directory where RAxML-V will write its output 
        # files.
        '-w':ValuedParameter('-',Name='w',Delimiter=' '),

        # THE FOLLOWING "W" is IN 7.2.8
        # Sliding window size for leave-one-out site-specific placement bias
        # algorithm only effective when used in combination with "-f S" 
        #   DEFAULT: 100 sites
        '-W':ValuedParameter('-',Name='W',Delimiter=' '),
        
        # Specify an integer number (random seed) and turn on rapid 
        # bootstrapping. CAUTION: unlike in version 7.0.4 RAxML will conduct 
        # rapid BS replicates under the model of rate heterogeneity you 
        # specified via "-m" and not by default under CAT
        '-x':ValuedParameter('-',Name='x',Delimiter=' '),
        
        # EXPERIMENTAL OPTION: This option will do a per-site estimate of
        # protein substitution models by looping over all given, fixed models
        # LG, WAG, JTT, etc and using their respective base frequencies to
        # independently assign a prot subst. model to each site via ML
        # optimization. At present this option only works with the GTR+GAMMA
        # model, unpartitioned datasets, and in the sequential version only.
        #   DEFAULT: OFF
        '-X':FlagParameter('-', Name='X'),

        # Compute only randomized starting parsimony tree with RAxML, do not
        # optimize an ML analysis of the tree
        '-y':FlagParameter('-', Name='y'),

        # Do a more thorough parsimony tree search using a parsimony ratchet and 
        # exit. Specify the number of ratchet searches via "-#" or "-N". This 
        # has just been implemented for completeness, if you want a fast MP 
        # implementation use TNT
        # DEFAULT: OFF - DEPRECATED IN 7.3.0
        #'-Y':FlagParameter('-', Name='Y'),

        # Multiple tree file, for use with -f b (to draw bipartitions onto the
        # common tree specified with -t)
        '-z':ValuedParameter('-', Name='z', Delimiter=' '),

        # Specifies number of runs on distinct starting trees.
        '-#':ValuedParameter('-', Name='#', Delimiter=' ',Value=1),

        # Specifies number of runs on distinct starting trees.
        '-N':ValuedParameter('-', Name='N', Delimiter=' '),

    }

    _parameters = {}
    _parameters.update(_options)
    _command = "raxmlHPC"
    _out_format = "RAxML_%s.%s"

    def _format_output(self, outfile_name, out_type):
        """ Prepend proper output prefix to output filename """

        outfile_name = self._absolute(outfile_name)
        outparts = outfile_name.split("/") 
        outparts[-1] = self._out_format % (out_type, outparts[-1] )

        return '/'.join(outparts)

    def _input_as_seqs(self,data):
        lines = []
        for i,s in enumerate(data):
            #will number the sequences 1,2,3,etc.
            lines.append(''.join(['>',str(i+1)]))
            lines.append(s)
        return self._input_as_lines(lines)

    def _input_as_lines(self,data):
        if data:
            self.Parameters['-s']\
                .on(super(Raxml,self)._input_as_lines(data))
        return ''

    def _input_as_string(self,data):
        """Makes data the value of a specific parameter
     
        This method returns the empty string. The parameter will be printed
        automatically once set.
        """
        if data:
            self.Parameters['-in'].on(str(data))
        return ''

    def _input_as_multiline_string(self, data):
        if data:
            self.Parameters['-s']\
                .on(super(Raxml,self)._input_as_multiline_string(data))
        return ''
   
    def _absolute(self,path):
        path = FilePath(path)
        if isabs(path):
            return path
        elif self.Parameters['-w'].isOn():
            return self.Parameters['-w'].Value + path
        else:
            return self.WorkingDir + path

    def _log_out_filename(self):
        if self.Parameters['-n'].isOn():
            return self._format_output(str(self.Parameters['-n'].Value), "log")
        else:
            raise ValueError, "No output file specified." 

    def _info_out_filename(self):
        if self.Parameters['-n'].isOn():
            return self._format_output(str(self.Parameters['-n'].Value), "info")
        else:
            raise ValueError, "No output file specified." 

    def _parsimony_tree_out_filename(self):
        if self.Parameters['-n'].isOn():
            return self._format_output(str(self.Parameters['-n'].Value), \
                                            "parsimonyTree")
        else:
            raise ValueError, "No output file specified." 
    
    # added for tree-insertion
    def _originallabelled_tree_out_filename(self):
        if self.Parameters['-n'].isOn():
            return self._format_output(str(self.Parameters['-n'].Value), \
                                            "originalLabelledTree")
        else:
            raise ValueError, "No output file specified."
    
    # added for tree-insertion
    def _labelled_tree_out_filename(self):
        if self.Parameters['-n'].isOn():
            return self._format_output(str(self.Parameters['-n'].Value), \
                                            "labelledTree")
        else:
            raise ValueError, "No output file specified."

    # added for tree-insertion
    def _classification_out_filename(self):
        if self.Parameters['-n'].isOn():
            return self._format_output(str(self.Parameters['-n'].Value), \
                                            "classification")
        else:
            raise ValueError, "No output file specified."
    
    # added for tree-insertion
    def _classificationlikelihoodweights_out_filename(self):
        if self.Parameters['-n'].isOn():
            return self._format_output(str(self.Parameters['-n'].Value), \
                                            "classificationLikelihoodWeights")
        else:
            raise ValueError, "No output file specified."
    
    # added for tree-insertion
    def _best_tree_out_filename(self):
        if self.Parameters['-n'].isOn():
            return self._format_output(str(self.Parameters['-n'].Value), \
                                            "bestTree")
        else:
            raise ValueError, "No output file specified."

    # added for tree-insertion
    def _entropy_out_filename(self):
        if self.Parameters['-n'].isOn():
            return self._format_output(str(self.Parameters['-n'].Value), \
                                            "entropy")
        else:
            raise ValueError, "No output file specified."

    # added for tree-insertion
    def _json_out_filename(self):
        if self.Parameters['-n'].isOn():
            return self._format_output(str(self.Parameters['-n'].Value), \
                                            "portableTree")
        else:
            raise ValueError, "No output file specified."
    
    # added for tree-insertion
    def _parsimony_out_filename(self):
        if self.Parameters['-n'].isOn():
            return self._format_output(str(self.Parameters['-n'].Value), \
                                            "equallyParsimoniousPlacements")
        else:
            raise ValueError, "No output file specified."
            
    def _result_tree_out_filename(self):
        if self.Parameters['-n'].isOn():
            return self._format_output(str(self.Parameters['-n'].Value), \
                                            "result")
        else:
            raise ValueError, "No output file specified." 

    def _result_bootstrap_out_filename(self):
        if self.Parameters['-n'].isOn():
            return self._format_output(str(self.Parameters['-n'].Value), \
                                            "bootstrap")
        else:
            raise ValueError, "No output file specified"

    def _checkpoint_out_filenames(self):
        """
        RAxML generates a crapload of checkpoint files so need to
        walk directory to collect names of all of them.
        """
        out_filenames = []
        if self.Parameters['-n'].isOn():
            out_name = str(self.Parameters['-n'].Value)
            walk_root = self.WorkingDir
            if self.Parameters['-w'].isOn(): 
                walk_root = str(self.Parameters['-w'].Value)
            for tup in walk(walk_root):
                dpath, dnames, dfiles = tup
                if dpath == walk_root:
                    for gen_file in dfiles:
                        if out_name in gen_file and "checkpoint" in gen_file:
                            out_filenames.append(walk_root + gen_file)
                    break

        else:
            raise ValueError, "No output file specified." 
        return out_filenames

    def _handle_app_result_build_failure(self,out,err,exit_status,result_paths):
        """ Catch the error when files are not produced """

        try:
            raise ApplicationError, \
             'RAxML failed to produce an output file due to the following error: \n\n%s ' \
             % err.read()
        except:
            raise ApplicationError,\
                'RAxML failed to run properly.'

    def _get_result_paths(self,data):

        result = {}
        result['Info'] = ResultPath(Path=self._info_out_filename(),
                                            IsWritten=True)
        if self.Parameters['-k'].isOn():
            result['Bootstrap'] = ResultPath(
                            Path=self._result_bootstrap_out_filename(),
                            IsWritten=True)
        elif self.Parameters["-f"].Value == 'v':
            #these were added to handle the results from tree-insertion
            result['Classification'] = ResultPath(
                Path=self._classification_out_filename(),
                IsWritten=True)
            result['ClassificationLikelihoodWeights'] = ResultPath(  
                Path=self._classificationlikelihoodweights_out_filename(),
                IsWritten=True)
            result['OriginalLabelledTree'] = ResultPath(  
                Path=self._originallabelled_tree_out_filename(),
                IsWritten=True)
            result['Result'] = ResultPath(
                Path=self._labelled_tree_out_filename(),IsWritten=True)
            result['entropy'] = ResultPath(
                Path=self._entropy_out_filename(),IsWritten=True)
            result['json'] = ResultPath(
                Path=self._json_out_filename()+'.jplace',IsWritten=True)
        elif self.Parameters["-f"].Value == 'y':
            #these were added to handle the results from tree-insertion
            
            result['Parsimony'] = ResultPath(  
                Path=self._parsimony_out_filename(),
                IsWritten=True)
            result['OriginalLabelledTree'] = ResultPath(  
                Path=self._originallabelled_tree_out_filename(),
                IsWritten=True)
            result['json'] = ResultPath(
                Path=self._json_out_filename()+'.jplace',IsWritten=True)
        else:
            result['Log'] = ResultPath(Path=self._log_out_filename(),
                                            IsWritten=True)
            result['ParsimonyTree'] = ResultPath(
                                      Path=self._parsimony_tree_out_filename(),
                                      IsWritten=True)
            result['Result'] = ResultPath(
                            Path=self._result_tree_out_filename(),
                            IsWritten=True)
            #
            result['besttree'] = ResultPath(
                            Path=self._best_tree_out_filename(),
                            IsWritten=True)
        
        for checkpoint_file in self._checkpoint_out_filenames():
            checkpoint_num = checkpoint_file.split(".")[-1]
            try:
                checkpoint_num = int(checkpoint_num)
            except Exception, e:
                raise ValueError, "%s does not appear to be a valid checkpoint file"
            result['Checkpoint%d' % checkpoint_num] = ResultPath(
                        Path=checkpoint_file,
                        IsWritten=True)
 
        return result


#SOME FUNCTIONS TO EXECUTE THE MOST COMMON TASKS
def raxml_alignment(align_obj,
                 raxml_model="GTRCAT",
                 params={},
                 SuppressStderr=True,
                 SuppressStdout=True):
    """Run raxml on alignment object 

    align_obj: Alignment object
    params: you can set any params except -w and -n

    returns: tuple (phylonode, 
                    parsimonyphylonode, 
                    log likelihood, 
                    total exec time)
    """

    # generate temp filename for output
    params["-w"] = "/tmp/"
    params["-n"] = get_tmp_filename().split("/")[-1]
    params["-m"] = raxml_model
    params["-p"] = randint(1,100000)
    ih = '_input_as_multiline_string'
    seqs, align_map = align_obj.toPhylip()
    
    #print params["-n"]

    # set up command
    raxml_app = Raxml(
                   params=params,
                   InputHandler=ih,
                   WorkingDir=None,
                   SuppressStderr=SuppressStderr,
                   SuppressStdout=SuppressStdout)

    # run raxml
    ra = raxml_app(seqs)

    # generate tree
    tree_node =  DndParser(ra["Result"])

    # generate parsimony tree
    parsimony_tree_node =  DndParser(ra["ParsimonyTree"])

    # extract log likelihood from log file
    log_file = ra["Log"]
    total_exec_time = exec_time = log_likelihood = 0.0
    for line in log_file:
        exec_time, log_likelihood = map(float, line.split())
        total_exec_time += exec_time

    # remove output files
    ra.cleanUp()

    return tree_node, parsimony_tree_node, log_likelihood, total_exec_time

def build_tree_from_alignment(aln, moltype, best_tree=False, params={}):
    """Returns a tree from Alignment object aln.
    
    aln: an xxx.Alignment object, or data that can be used to build one.
    
    moltype: cogent.core.moltype.MolType object

    best_tree: best_tree suppport is currently not implemented
    
    params: dict of parameters to pass in to the RAxML app controller.
    
    The result will be an xxx.Alignment object, or None if tree fails.
    """
    if best_tree:
        raise NotImplementedError

    if '-m' not in params:
        if moltype == DNA or moltype == RNA:
            #params["-m"] = 'GTRMIX'
            # in version 7.2.3, GTRMIX is no longer supported but says GTRCAT
            # behaves like GTRMIX (http://www.phylo.org/tools/raxmlhpc2.html)
            params["-m"] = 'GTRGAMMA'
        elif moltype == PROTEIN:
            params["-m"] = 'PROTGAMMAmatrixName'
        else:
            raise ValueError, "Moltype must be either DNA, RNA, or PROTEIN"

    if not hasattr(aln, 'toPhylip'):
        aln = Alignment(aln)
    seqs, align_map = aln.toPhylip()

    # generate temp filename for output    
    params["-w"] = "/tmp/"    
    params["-n"] = get_tmp_filename().split("/")[-1]
    params["-k"] = True
    params["-p"] = randint(1,100000)
    params["-x"] = randint(1,100000)
    
    ih = '_input_as_multiline_string'    

    raxml_app = Raxml(params=params,
                      InputHandler=ih,
                      WorkingDir=None,
                      SuppressStderr=True,
                      SuppressStdout=True)
                      
    raxml_result = raxml_app(seqs)
    
    tree = DndParser(raxml_result['Bootstrap'], constructor=PhyloNode)
    
    for node in tree.tips():
        node.Name = align_map[node.Name]

    raxml_result.cleanUp()

    return tree
    
    
def insert_sequences_into_tree(seqs, moltype, params={},
                                           write_log=True):
    """Insert sequences into Tree.
    
    aln: an xxx.Alignment object, or data that can be used to build one.
    
    moltype: cogent.core.moltype.MolType object
    
    params: dict of parameters to pass in to the RAxML app controller.
    
    The result will be a tree.
    """
    
    ih = '_input_as_multiline_string'    

    raxml_app = Raxml(params=params,
                      InputHandler=ih,
                      WorkingDir=None,
                      SuppressStderr=False,
                      SuppressStdout=False,
                      HALT_EXEC=False)
    
    raxml_result = raxml_app(seqs)
    
    # write a log file
    if write_log:
        log_fp = join(params["-w"],'log_raxml_'+split(get_tmp_filename())[-1])
        log_file=open(log_fp,'w')
        log_file.write(raxml_result['StdOut'].read())
        log_file.close()
    
    ''' 
    # getting setup since parsimony doesn't output tree..only jplace, however
    # it is currently corrupt
        
    # use guppy to convert json file into a placement tree
    guppy_params={'tog':None}

    new_tree=build_tree_from_json_using_params(raxml_result['json'].name, \
                                               output_dir=params["-w"], \
                                               params=guppy_params)
    '''
    
    # get tree from 'Result Names'
    new_tree=raxml_result['Result'].readlines()
    filtered_tree=re.sub('\[I\d+\]','',str(new_tree))
    tree = DndParser(filtered_tree, constructor=PhyloNode)

    raxml_result.cleanUp()

    return tree



