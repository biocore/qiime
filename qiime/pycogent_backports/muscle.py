#!/usr/bin/env python
"""Application controller for muscle 3.6
"""
from os import remove
from cogent.app.parameters import FlagParameter, ValuedParameter
from cogent.app.util import CommandLineApplication, ResultPath, \
    get_tmp_filename, guess_input_handler
from random import choice
from cogent.core.alignment import SequenceCollection, Alignment
from cogent.parse.tree import DndParser
from cogent.core.tree import PhyloNode
from cogent.parse.fasta import MinimalFastaParser

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2011, The Cogent Project"
__credits__ = ["Micah Hamady", "Zongzhi Liu", "Mike Robeson",
       "Catherine Lozupone", "Rob Knight", "Daniel McDonald", "Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.6.0dev"
__maintainer__ = "Micah Hamady"
__email__ = "hamady@colorado.edu"
__status__ = "Prototype"

class Muscle(CommandLineApplication):
    """Muscle application controller"""
    
    _options ={
        # Minimum spacing between anchor columns. [Integer]
        '-anchorspacing':ValuedParameter('-',Name='anchorspacing',Delimiter=' '),
        # Center parameter. Should be negative [Float]
        '-center':ValuedParameter('-',Name='center',Delimiter=' '),
        
        # Clustering method. cluster1 is used in iteration 1
        # and 2, cluster2 in later iterations
        '-cluster1':ValuedParameter('-',Name='cluster1',Delimiter=' '),
        '-cluster2':ValuedParameter('-',Name='cluster2',Delimiter=' '),
        
        # Minimum length of diagonal.
        '-diaglength':ValuedParameter('-',Name='diaglength',Delimiter=' '),
        
        # Discard this many positions at ends of diagonal.
        '-diagmargin':ValuedParameter('-',Name='diagmargin',Delimiter=' '),
        
        # Distance measure for iteration 1.
        '-distance1':ValuedParameter('-',Name='distance1',Delimiter=' '),
        
        # Distance measure for iterations 2, 3 ...
        '-distance2':ValuedParameter('-',Name='distance2',Delimiter=' '),
        
        # The gap open score. Must be negative.
        '-gapopen':ValuedParameter('-',Name='gapopen',Delimiter=' '),
        
        # Window size for determining whether a region is hydrophobic.
        '-hydro':ValuedParameter('-',Name='hydro',Delimiter=' '),
        
        # Multiplier for gap open/close penalties in hydrophobic regions.
        '-hydrofactor':ValuedParameter('-',Name='hydrofactor',Delimiter=' '),
        
        # Where to find the input sequences.
        '-in':ValuedParameter('-',Name='in',Delimiter=' ', Quote="\""),
        '-in1':ValuedParameter('-',Name='in1',Delimiter=' ', Quote="\""),
        '-in2':ValuedParameter('-',Name='in2',Delimiter=' ', Quote="\""),
        
        # Log file name (delete existing file).
        '-log':ValuedParameter('-',Name='log',Delimiter=' '),
        
        # Log file name (append to existing file).
        '-loga':ValuedParameter('-',Name='loga',Delimiter=' '),
        
        # Maximum distance between two diagonals that allows them to merge
        # into one diagonal.
        '-maxdiagbreak':ValuedParameter('-',Name='maxdiagbreak',Delimiter=' '),
        
        # Maximum time to run in hours. The actual time may exceed the
        # requested limit by a few minutes. Decimals are allowed, so 1.5
        # means one hour and 30 minutes.
        '-maxhours':ValuedParameter('-',Name='maxhours',Delimiter=' '),
        
        # Maximum number of iterations.
        '-maxiters':ValuedParameter('-',Name='maxiters',Delimiter=' '),

        # Maximum memory in Mb
        '-maxmb': ValuedParameter('-', Name='maxmb', Delimiter=' '),
        
        # Maximum number of new trees to build in iteration 2.
        '-maxtrees':ValuedParameter('-',Name='maxtrees',Delimiter=' '),
        
        # Minimum score a column must have to be an anchor.
        '-minbestcolscore':ValuedParameter('-',Name='minbestcolscore',Delimiter=' '),
        
        # Minimum smoothed score a column must have to be an anchor.
        '-minsmoothscore':ValuedParameter('-',Name='minsmoothscore',Delimiter=' '),
        
        # Objective score used by tree dependent refinement.
        # sp=sum-of-pairs score.
        # spf=sum-of-pairs score (dimer approximation)
        # spm=sp for < 100 seqs, otherwise spf
        # dp=dynamic programming score.
        # ps=average profile-sequence score.
        # xp=cross profile score.
        '-objscore':ValuedParameter('-',Name='objscore',Delimiter=' '),
        
        # Where to write the alignment.
        '-out':ValuedParameter('-',Name='out',Delimiter=' ', Quote="\""),
        
        # Where to write the file in phylip sequenctial format (v3.6 only).
        '-physout':ValuedParameter('-',Name='physout',Delimiter=' '),
        
        # Where to write the file in phylip interleaved format (v3.6 only).
        '-phyiout':ValuedParameter('-',Name='phyiout',Delimiter=' '),

        # Set to profile for aligning two alignments and adding seqs to an 
        # existing alignment
        '-profile':FlagParameter(Prefix='-',Name='profile'),

        # Method used to root tree; root1 is used in iteration 1 and 2, root2
        # in later iterations.
        '-root1':ValuedParameter('-',Name='root1',Delimiter=' '),
        '-root2':ValuedParameter('-',Name='root2',Delimiter=' '),
        
        # Sequence type.
        '-seqtype':ValuedParameter('-',Name='seqtype',Delimiter=' '),
        
        # Maximum value of column score for smoothing purposes.
        '-smoothscoreceil':ValuedParameter('-',Name='smoothscoreceil',Delimiter=' '),
        
        # Constant used in UPGMB clustering. Determines the relative fraction
        # of average linkage (SUEFF) vs. nearest-neighbor linkage (1 . SUEFF).
        '-SUEFF':ValuedParameter('-',Name='SUEFF',Delimiter=' '),
        
        # Save tree produced in first or second iteration to given file in
        # Newick (Phylip-compatible) format.
        '-tree1':ValuedParameter('-',Name='tree1',Delimiter=' ', Quote="\""),
        '-tree2':ValuedParameter('-',Name='tree2',Delimiter=' ', Quote="\""),
        
        # Sequence weighting scheme.
        # weight1 is used in iterations 1 and 2.
        # weight2 is used for tree-dependent refinement.
        # none=all sequences have equal weight.
        # henikoff=Henikoff & Henikoff weighting scheme.
        # henikoffpb=Modified Henikoff scheme as used in PSI-BLAST.
        # clustalw=CLUSTALW method.
        # threeway=Gotoh three-way method.
        '-weight1':ValuedParameter('-',Name='weight1',Delimiter=' '),
        '-weight2':ValuedParameter('-',Name='weight2',Delimiter=' '),
        
        # Use anchor optimization in tree dependent refinement iterations
        '-anchors':FlagParameter(Prefix='-',Name='anchors'),
        
        # Write output in CLUSTALW format (default is FASTA).
        '-clw':FlagParameter(Prefix='-',Name='clw'),
        
        # Cluster sequences
        '-cluster':FlagParameter(Prefix='-',Name='cluster'),
        # neighborjoining is "unrecognized"
        #'-neighborjoining':FlagParameter(Prefix='-',Name='neighborjoining'),

        
        # Write output in CLUSTALW format with the "CLUSTAL W (1.81)" header
        # rather than the MUSCLE version. This is useful when a post-processing
        # step is picky about the file header.
        '-clwstrict':FlagParameter(Prefix='-',Name='clwstrict'),
        
        # Do not catch exceptions.
        '-core':FlagParameter(Prefix='-',Name='core'),
        
        # Write output in FASTA format. Alternatives include .clw,
        # .clwstrict, .msf and .html.
        '-fasta':FlagParameter(Prefix='-',Name='fasta'),
        
        # Group similar sequences together in the output. This is the default.
        # See also .stable.
        '-group':FlagParameter(Prefix='-',Name='group'),
        
        # Write output in HTML format (default is FASTA).
        '-html':FlagParameter(Prefix='-',Name='html'),
        
        # Use log-expectation profile score (VTML240). Alternatives are to use
        # -sp or -sv. This is the default for amino acid sequences.
        '-le':FlagParameter(Prefix='-',Name='le'),
        
        # Write output in MSF format (default is FASTA).
        '-msf':FlagParameter(Prefix='-',Name='msf'),
        
        # Disable anchor optimization. Default is -anchors.
        '-noanchors':FlagParameter(Prefix='-',Name='noanchors'),
        
        # Catch exceptions and give an error message if possible.
        '-nocore':FlagParameter(Prefix='-',Name='nocore'),
        
        # Do not display progress messages.
        '-quiet':FlagParameter(Prefix='-',Name='quiet'),
        
        # Input file is already aligned, skip first two iterations and begin
        # tree dependent refinement.
        '-refine':FlagParameter(Prefix='-',Name='refine'),
        
        # Use sum-of-pairs protein profile score (PAM200). Default is -le.
        '-sp':FlagParameter(Prefix='-',Name='sp'),
        
        # Use sum-of-pairs nucleotide profile score (BLASTZ parameters). This
        # is the only option for nucleotides, and is therefore the default.
        '-spn':FlagParameter(Prefix='-',Name='spn'),
        
        # Preserve input order of sequences in output file. Default is to group
        # sequences by similarity (-group).
        '-stable':FlagParameter(Prefix='-',Name='stable'),
        
        # Use sum-of-pairs profile score (VTML240). Default is -le.
        '-sv':FlagParameter(Prefix='-',Name='sv'),
        
        # Diagonal optimization
        '-diags':FlagParameter(Prefix='-',Name='diags'),
        '-diags1':FlagParameter(Prefix='-',Name='diags1'),
        '-diags2':FlagParameter(Prefix='-',Name='diags2'),

        
        # Terminal gaps penalized with full penalty.
        # [1] Not fully supported in this version.
        '-termgapsfull':FlagParameter(Prefix='-',Name='termgapsfull'),
        
        # Terminal gaps penalized with half penalty.
        # [1] Not fully supported in this version.
        '-termgapshalf':FlagParameter(Prefix='-',Name='termgapshalf'),
        
        # Terminal gaps penalized with half penalty if gap relative to
        # longer sequence, otherwise with full penalty.
        # [1] Not fully supported in this version.
        '-termgapshalflonger':FlagParameter(Prefix='-',Name='termgapshalflonger'),
        
        # Write parameter settings and progress messages to log file.
        '-verbose':FlagParameter(Prefix='-',Name='verbose'),
        
        # Write version string to stdout and exit.
        '-version':FlagParameter(Prefix='-',Name='version'),
    }
    
    _parameters = {}
    _parameters.update(_options)
    _command = "muscle"
    
    def _input_as_seqs(self,data):
        lines = []
        for i,s in enumerate(data):
            #will number the sequences 1,2,3,etc.
            lines.append(''.join(['>',str(i+1)]))
            lines.append(s)
        return self._input_as_lines(lines)
    
    def _input_as_lines(self,data):
        if data:
            self.Parameters['-in']\
                .on(super(Muscle,self)._input_as_lines(data))
        
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
            self.Parameters['-in']\
                .on(super(Muscle,self)._input_as_multiline_string(data))
        return ''

    def _input_as_multifile(self, data):
        """For use with the -profile option

        This input handler expects data to be a tuple containing two
        filenames. Index 0 will be set to -in1 and index 1 to -in2
        """
        if data:
            try:
                filename1, filename2 = data
            except:
                raise ValueError, "Expected two filenames"

            self.Parameters['-in'].off()
            self.Parameters['-in1'].on(filename1)
            self.Parameters['-in2'].on(filename2)
        return ''

    def _align_out_filename(self):
        
        if self.Parameters['-out'].isOn():
            aln_filename = self._absolute(str(self.Parameters['-out'].Value))
        else:
            raise ValueError, "No output file specified."
        return aln_filename
    
    def _tree1_out_filename(self):
        
        if self.Parameters['-tree1'].isOn():
            aln_filename = self._absolute(str(self.Parameters['-tree1'].Value))
        else:
            raise ValueError, "No tree output file specified."
        return aln_filename
    
    def _tree2_out_filename(self):
        
        if self.Parameters['-tree2'].isOn():
            tree_filename = self._absolute(str(self.Parameters['-tree2'].Value))
        else:
            raise ValueError, "No tree output file specified."
        return tree_filename
    
    def _get_result_paths(self,data):
        
        result = {}
        if self.Parameters['-out'].isOn():
            out_name = self._align_out_filename()
            result['MuscleOut'] = ResultPath(Path=out_name,IsWritten=True)
        if self.Parameters['-tree1'].isOn():
            out_name = self._tree1_out_filename()
            result['Tree1Out'] = ResultPath(Path=out_name,IsWritten=True)
        if self.Parameters['-tree2'].isOn():
            out_name = self._tree2_out_filename()
            result['Tree2Out'] = ResultPath(Path=out_name,IsWritten=True)
        return result

    
    def getHelp(self):
        """Muscle help"""
        
        help_str = """
"""
        return help_str

#SOME FUNCTIONS TO EXECUTE THE MOST COMMON TASKS
def muscle_seqs(seqs,
                 add_seq_names=False,
                 out_filename=None,
                 input_handler=None,
                 params={},
                 WorkingDir=None,
                 SuppressStderr=None,
                 SuppressStdout=None):
    """Muscle align list of sequences.
    
    seqs: a list of sequences as strings or objects, you must set add_seq_names=True
    or sequences in a multiline string, as read() from a fasta file
    or sequences in a list of lines, as readlines() from a fasta file
    or a fasta seq filename.
    
    == for eg, testcode for guessing
        #guess_input_handler should correctly identify input
        gih = guess_input_handler
        self.assertEqual(gih('abc.txt'), '_input_as_string')
        self.assertEqual(gih('>ab\nTCAG'), '_input_as_multiline_string')
        self.assertEqual(gih(['ACC','TGA'], True), '_input_as_seqs')
        self.assertEqual(gih(['>a','ACC','>b','TGA']), '_input_as_lines')
    
    == docstring for blast_seqs, apply to muscle_seqs ==
    seqs: either file name or list of sequence objects or list of strings or
    single multiline string containing sequences.
    
    WARNING: DECISION RULES FOR INPUT HANDLING HAVE CHANGED. Decision rules
    for data are as follows. If it's s list, treat as lines, unless
    add_seq_names is true (in which case treat as list of seqs). If it's a
    string, test whether it has newlines. If it doesn't have newlines, assume
    it's a filename. If it does have newlines, it can't be a filename, so
    assume it's a multiline string containing sequences.
    
    If you want to skip the detection and force a specific type of input
    handler, use input_handler='your_favorite_handler'.
    
    add_seq_names: boolean. if True, sequence names are inserted in the list
        of sequences. if False, it assumes seqs is a list of lines of some
        proper format that the program can handle
    
    Addl docs coming soon
    """
    
    if out_filename:
        params["-out"] = out_filename
    #else:
    #    params["-out"] = get_tmp_filename(WorkingDir)
    
    ih = input_handler or guess_input_handler(seqs, add_seq_names)
    muscle_app = Muscle(
                   params=params,
                   InputHandler=ih,
                   WorkingDir=WorkingDir,
                   SuppressStderr=SuppressStderr,
                   SuppressStdout=SuppressStdout)
    return muscle_app(seqs)


def cluster_seqs(seqs,
                 neighbor_join=False,
                 params={},
                 add_seq_names=True,
                 WorkingDir=None,
                 SuppressStderr=None,
                 SuppressStdout=None,
                 max_chars=1000000,
                 max_hours=1.0,
                 constructor=PhyloNode,
                 clean_up=True
                 ):
    """Muscle cluster list of sequences.
    
    seqs: either file name or list of sequence objects or list of strings or
        single multiline string containing sequences.
    
    Addl docs coming soon
    """
    num_seqs = len(seqs)
    if num_seqs < 2:
        raise ValueError, "Muscle requres 2 or more sequences to cluster."

    
    num_chars = sum(map(len, seqs))
    if num_chars > max_chars:
        params["-maxiters"] = 2
        params["-diags1"] = True
        params["-sv"] = True
        #params["-distance1"] = "kmer6_6"
        #params["-distance1"] = "kmer20_3"
        #params["-distance1"] = "kbit20_3"
        print "lots of chars, using fast align", num_chars

    
    params["-maxhours"] = max_hours
    #params["-maxiters"] = 10
    
    #cluster_type = "upgmb"
    #if neighbor_join:
    #    cluster_type = "neighborjoining"
    
    params["-cluster"] = True
    params["-tree1"] = get_tmp_filename(WorkingDir)
    
    muscle_res = muscle_seqs(seqs,
                 params=params,
                 add_seq_names=add_seq_names,
                 WorkingDir=WorkingDir,
                 SuppressStderr=SuppressStderr,
                 SuppressStdout=SuppressStdout)
    
    tree = DndParser(muscle_res["Tree1Out"], constructor=constructor)
    
    if clean_up:
        muscle_res.cleanUp()
    return tree

def aln_tree_seqs(seqs,
                 input_handler=None,
                 tree_type='neighborjoining',
                 params={},
                 add_seq_names=True,
                 WorkingDir=None,
                 SuppressStderr=None,
                 SuppressStdout=None,
                 max_hours=5.0,
                 constructor=PhyloNode,
                 clean_up=True
                 ):
    """Muscle align sequences and report tree from iteration2.
    
    Unlike cluster_seqs, returns tree2 which is the tree made during the
    second muscle iteration (it should be more accurate that the cluster from
    the first iteration which is made fast based on  k-mer words)
    
    seqs: either file name or list of sequence objects or list of strings or
        single multiline string containing sequences.
    tree_type: can be either neighborjoining (default) or upgmb for UPGMA
    clean_up: When true, will clean up output files
    """
    
    params["-maxhours"] = max_hours
    if tree_type:
        params["-cluster2"] = tree_type
    params["-tree2"] = get_tmp_filename(WorkingDir)
    params["-out"] = get_tmp_filename(WorkingDir)
    
    muscle_res = muscle_seqs(seqs,
                 input_handler=input_handler,
                 params=params,
                 add_seq_names=add_seq_names,
                 WorkingDir=WorkingDir,
                 SuppressStderr=SuppressStderr,
                 SuppressStdout=SuppressStdout)
    tree = DndParser(muscle_res["Tree2Out"], constructor=constructor)
    aln = [line for line in muscle_res["MuscleOut"]]
    
    if clean_up:
        muscle_res.cleanUp()
    return tree, aln

def fastest_aln_seqs(seqs,
                 params={},
                 out_filename=None,
                 add_seq_names=True,
                 WorkingDir=None,
                 SuppressStderr=None,
                 SuppressStdout=None
                 ):
    """Fastest (and least accurate) version of muscle
    
    seqs: either file name or list of sequence objects or list of strings or
        single multiline string containing sequences.
    
    Addl docs coming soon
    """
    
    params["-maxiters"] = 1
    params["-diags1"] = True
    params["-sv"] = True
    params["-distance1"] = "kbit20_3"
    
    muscle_res = muscle_seqs(seqs,
                 params=params,
                 add_seq_names=add_seq_names,
                 out_filename=out_filename,
                 WorkingDir=WorkingDir,
                 SuppressStderr=SuppressStderr,
                 SuppressStdout=SuppressStdout)
    return muscle_res

def align_unaligned_seqs(seqs, moltype, params=None):
    """Returns an Alignment object from seqs.

    seqs: SequenceCollection object, or data that can be used to build one.
    
    moltype: a MolType object.  DNA, RNA, or PROTEIN.

    params: dict of parameters to pass in to the Muscle app controller.
    
    Result will be an Alignment object.
    """
    if not params:
        params = {}
    #create SequenceCollection object from seqs
    seq_collection = SequenceCollection(seqs,MolType=moltype)
    #Create mapping between abbreviated IDs and full IDs
    int_map, int_keys = seq_collection.getIntMap()
    #Create SequenceCollection from int_map.
    int_map = SequenceCollection(int_map,MolType=moltype)
    #get temporary filename
    params.update({'-out':get_tmp_filename()})
    #Create Muscle app.
    app = Muscle(InputHandler='_input_as_multiline_string',\
                 params=params)
    #Get results using int_map as input to app
    res = app(int_map.toFasta())
    #Get alignment as dict out of results
    alignment = dict(MinimalFastaParser(res['MuscleOut'].readlines()))
    #Make new dict mapping original IDs
    new_alignment = {}
    for k,v in alignment.items():
        new_alignment[int_keys[k]]=v
    #Create an Alignment object from alignment dict
    new_alignment = Alignment(new_alignment,MolType=moltype)
    #Clean up
    res.cleanUp()
    del(seq_collection,int_map,int_keys,app,res,alignment,params)

    return new_alignment


def align_and_build_tree(seqs, moltype, best_tree=False, params=None):
    """Returns an alignment and a tree from Sequences object seqs.
    
    seqs: a cogent.core.alignment.SequenceCollection object, or data that can 
    be used to build one.
    
    moltype: cogent.core.moltype.MolType object

    best_tree: if True (default:False), uses a slower but more accurate
    algorithm to build the tree.
    
    params: dict of parameters to pass in to the Muscle app controller.
    
    The result will be a tuple containing a cogent.core.alignment.Alignment 
    and a cogent.core.tree.PhyloNode object (or None for the alignment 
    and/or tree if either fails).
    """
    aln = align_unaligned_seqs(seqs, moltype=moltype, params=params)
    tree = build_tree_from_alignment(aln, moltype, best_tree, params)
    return {'Align':aln, 'Tree':tree}

def build_tree_from_alignment(aln, moltype, best_tree=False, params=None):
    """Returns a tree from Alignment object aln.
    
    aln: a cogent.core.alignment.Alignment object, or data that can be used 
    to build one.
    
    moltype: cogent.core.moltype.MolType object

    best_tree: unsupported
    
    params: dict of parameters to pass in to the Muscle app controller.
    
    The result will be an cogent.core.tree.PhyloNode object, or None if tree 
    fails.
    """
    # Create instance of app controller, enable tree, disable alignment
    app = Muscle(InputHandler='_input_as_multiline_string', params=params, \
                   WorkingDir='/tmp')

    app.Parameters['-cluster'].on()
    app.Parameters['-tree1'].on(get_tmp_filename(app.WorkingDir))
    app.Parameters['-seqtype'].on(moltype.label)

    seq_collection = SequenceCollection(aln, MolType=moltype)

    #Create mapping between abbreviated IDs and full IDs
    int_map, int_keys = seq_collection.getIntMap()
    #Create SequenceCollection from int_map.
    int_map = SequenceCollection(int_map,MolType=moltype)


    # Collect result
    result = app(int_map.toFasta())

    # Build tree
    tree = DndParser(result['Tree1Out'].read(), constructor=PhyloNode)
    
    for tip in tree.tips():
        tip.Name = int_keys[tip.Name]

    # Clean up
    result.cleanUp()
    del(seq_collection, app, result)

    return tree

def add_seqs_to_alignment(seqs, aln, params=None):
    """Returns an Alignment object from seqs and existing Alignment.
    
    seqs: a cogent.core.alignment.SequenceCollection object, or data that can 
    be used to build one.
    
    aln: a cogent.core.alignment.Alignment object, or data that can be used 
    to build one
    
    params: dict of parameters to pass in to the Muscle app controller.
    """
    if not params:
        params = {}

    #create SequenceCollection object from seqs
    seqs_collection = SequenceCollection(seqs)
    #Create mapping between abbreviated IDs and full IDs
    seqs_int_map, seqs_int_keys = seqs_collection.getIntMap(prefix='seq_')
    #Create SequenceCollection from int_map.
    seqs_int_map = SequenceCollection(seqs_int_map)

    #create SequenceCollection object from aln
    aln_collection = SequenceCollection(aln)
    #Create mapping between abbreviated IDs and full IDs
    aln_int_map, aln_int_keys = aln_collection.getIntMap(prefix='aln_')
    #Create SequenceCollection from int_map.
    aln_int_map = SequenceCollection(aln_int_map)

    #set output and profile options
    params.update({'-out':get_tmp_filename(), '-profile':True})

    #save seqs to tmp file
    seqs_filename = get_tmp_filename()
    seqs_out = open(seqs_filename,'w')
    seqs_out.write(seqs_int_map.toFasta())
    seqs_out.close()

    #save aln to tmp file
    aln_filename = get_tmp_filename()
    aln_out = open(aln_filename, 'w')
    aln_out.write(aln_int_map.toFasta())
    aln_out.close()

    #Create Muscle app and get results
    app = Muscle(InputHandler='_input_as_multifile', params=params)
    res = app((aln_filename, seqs_filename))

    #Get alignment as dict out of results
    alignment = dict(MinimalFastaParser(res['MuscleOut'].readlines()))
    #Make new dict mapping original IDs
    new_alignment = {}
    for k,v in alignment.items():
        if k in seqs_int_keys:
            new_alignment[seqs_int_keys[k]] = v
        else:
            new_alignment[aln_int_keys[k]] = v

    #Create an Alignment object from alignment dict
    new_alignment = Alignment(new_alignment)

    #Clean up
    res.cleanUp()
    del(seqs_collection, seqs_int_map, seqs_int_keys)
    del(aln_collection, aln_int_map, aln_int_keys)
    del(app, res, alignment, params)
    remove(seqs_filename)
    remove(aln_filename)

    return new_alignment

def align_two_alignments(aln1, aln2, params=None):
    """Returns an Alignment object from two existing Alignments.
    
    aln1, aln2: cogent.core.alignment.Alignment objects, or data that can be 
    used to build them.
    
    params: dict of parameters to pass in to the Muscle app controller.
    """
    if not params:
        params = {}

    #create SequenceCollection object from aln1
    aln1_collection = SequenceCollection(aln1)
    #Create mapping between abbreviated IDs and full IDs
    aln1_int_map, aln1_int_keys = aln1_collection.getIntMap(prefix='aln1_')
    #Create SequenceCollection from int_map.
    aln1_int_map = SequenceCollection(aln1_int_map)

    #create SequenceCollection object from aln2
    aln2_collection = SequenceCollection(aln2)
    #Create mapping between abbreviated IDs and full IDs
    aln2_int_map, aln2_int_keys = aln2_collection.getIntMap(prefix='aln2_')
    #Create SequenceCollection from int_map.
    aln2_int_map = SequenceCollection(aln2_int_map)

    #set output and profile options
    params.update({'-out':get_tmp_filename(), '-profile':True})

    #save aln1 to tmp file
    aln1_filename = get_tmp_filename()
    aln1_out = open(aln1_filename,'w')
    aln1_out.write(aln1_int_map.toFasta())
    aln1_out.close()

    #save aln2 to tmp file
    aln2_filename = get_tmp_filename()
    aln2_out = open(aln2_filename, 'w')
    aln2_out.write(aln2_int_map.toFasta())
    aln2_out.close()

    #Create Muscle app and get results
    app = Muscle(InputHandler='_input_as_multifile', params=params)
    res = app((aln1_filename, aln2_filename))

    #Get alignment as dict out of results
    alignment = dict(MinimalFastaParser(res['MuscleOut'].readlines()))

    #Make new dict mapping original IDs
    new_alignment = {}
    for k,v in alignment.items():
        if k in aln1_int_keys:
            new_alignment[aln1_int_keys[k]] = v
        else:
            new_alignment[aln2_int_keys[k]] = v

    #Create an Alignment object from alignment dict
    new_alignment = Alignment(new_alignment)

    #Clean up
    res.cleanUp()
    del(aln1_collection, aln1_int_map, aln1_int_keys)
    del(aln2_collection, aln2_int_map, aln2_int_keys)
    del(app, res, alignment, params)
    remove(aln1_filename)
    remove(aln2_filename)

    return new_alignment
