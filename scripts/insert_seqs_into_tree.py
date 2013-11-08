#!/usr/bin/env python
# File created on 11 Oct 2011
from __future__ import division

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Jesse Stombaugh", "Jai Ram Rideout", "Emily TerAvest"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Development"
 
from qiime.util import parse_command_line_parameters, make_option, \
                       get_options_lookup,load_qiime_config,create_dir
from cogent.parse.fasta import MinimalFastaParser
from cogent.core.alignment import DenseAlignment
from qiime.parse import parse_qiime_parameters
from cogent.core.moltype import DNA
from qiime.util import get_tmp_filename
from os.path import abspath,join,split,splitext
from qiime.insert_seqs_into_tree import convert_tree_tips, \
                                    write_updated_tree_file, \
                                    strip_and_rename_unwanted_labels_from_tree
import cogent.app.raxml_v730
import cogent.app.parsinsert
import cogent.app.pplacer
from StringIO import StringIO

options_lookup = get_options_lookup()

qiime_config = load_qiime_config()

insertion_method_choices = ['pplacer','raxml_v730','parsinsert']

script_info = {}
script_info['brief_description'] = "Tree Insertion"
script_info['script_description'] = "This script takes a set of aligned sequences (query) either in the same file as the aligned reference set or separated (depending on method) along with a starting tree and produces a new tree containing the query sequences. This script requires that the user is running Raxml v7.3.0, PPlacer git repository version and ParsInsert 1.0.4."
script_info['script_usage'] = []
script_info['script_usage'].append(("""RAxML Example (default):""","""If you just want to use the default options, you can supply an alignment files where the query and reference sequences are included, along with a starting tree as follows:""","""%prog -i aligned_query_seqs.fasta -r aligned_reference_seqs.fasta -t starting_tree.tre -o insertion_results"""))
script_info['script_usage'].append(("""ParsInsert Example:""","""If you want to insert sequences using pplacer, you can supply a fasta file containg query sequences (aligned to reference sequences) along with the reference alignment, a starting tree and the stats file produced when building the starting tree via pplacer as follows:""","""%prog -i aligned_query_seqs.fasta -r aligned_reference_seqs.fasta -t starting_tree.tre -o insertion_results -m parsinsert"""))
script_info['script_usage'].append(("""Pplacer Example:""","""If you want to insert sequences using pplacer, you can supply a fasta file containg query sequences (aligned to reference sequences) along with the reference alignment, a starting tree and the stats file produced when building the starting tree via pplacer as follows:""","""%prog -i aligned_query_seqs.fasta -r aligned_reference_seqs.fasta -t starting_tree.tre -o insertion_results -m pplacer"""))
script_info['script_usage'].append(("""Parameters file:""","""Additionally, users can supply a parameters file to change the options of the underlying tools as follows:""","""%prog -i aligned_query_seqs.fasta -r aligned_reference_seqs.fasta -t starting_tree.tre -o insertion_results -p raxml_parameters.txt"""))
script_info['output_description']= "The result of this script produces a tree file (in Newick format) along with a log file containing the output from the underlying tool used for tree insertion."
script_info['required_options'] = [\
    options_lookup['fasta_as_primary_input'],
    options_lookup['output_dir'],
    make_option('-t','--starting_tree_fp',\
             type='existing_filepath',help='Starting Tree which you would like to insert into.'),
    make_option('-r','--refseq_fp',\
          type='existing_filepath',dest='refseq_fp',help='Filepath for '+\
          'reference alignment'),
]
script_info['optional_options'] = [\
    make_option('-m','--insertion_method',\
          type='choice',help='Method for aligning'+\
          ' sequences. Valid choices are: ' +\
          ', '.join(insertion_method_choices) + ' [default: %default]',
          choices=insertion_method_choices,\
          default='raxml_v730'),
    make_option('-s','--stats_fp',\
          type='existing_filepath',help='Stats file produced by tree-building software. REQUIRED if -m pplacer [default: %default]'),
    make_option('-p','--method_params_fp',\
            type='existing_filepath',help='Parameters file containing method-specific parameters to use. Lines should be formatted as 'raxml:-m GTRCAT' (note this is not a standard QIIME parameters file, but a RAxML parameters file). [default: %default]'),

]
script_info['version'] = __version__

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
      
    parameters={}
    
    # get the tree insertion method to use
    module = opts.insertion_method
    
    # create output directory
    output_dir=opts.output_dir
    create_dir(output_dir)
    
    # list of tree insertion methods
    tree_insertion_module_names = \
                {'raxml_v730':cogent.app.raxml_v730, 
                 'parsinsert':cogent.app.parsinsert,
                 'pplacer':cogent.app.pplacer}

    # load input sequences and convert to phylip since the tools require 
    # the query sequences to phylip-compliant names
    load_aln = MinimalFastaParser(open(opts.input_fasta_fp,'U'))   
    aln = DenseAlignment(load_aln)
    seqs, align_map = aln.toPhylip()
    
    if opts.method_params_fp:
        param_dict = parse_qiime_parameters(open(opts.method_params_fp,'U'))
        
    if module=='raxml_v730':
        # load the reference sequences
        load_ref_aln = \
            DenseAlignment(MinimalFastaParser(open(opts.refseq_fp,'U'))) 
        
        # combine and load the reference plus query
        combined_aln = MinimalFastaParser(StringIO(load_ref_aln.toFasta() + \
                                                   '\n' + aln.toFasta()))
        # overwrite the alignment map
        aln = DenseAlignment(combined_aln)
        seqs, align_map = aln.toPhylip()
        
        try: 
            parameters = param_dict['raxml']
        except:
            parameters = {}
            
        tree = convert_tree_tips(align_map,opts.starting_tree_fp)
        
        # write out the tree with phylip labels
        updated_tree_fp = join(output_dir, \
                                '%s_phylip_named_tree.tre' % (module))
        write_updated_tree_file(updated_tree_fp,tree)
        
        # set the primary parameters for raxml
        parameters['-w'] = abspath(output_dir)+'/'
        parameters["-n"] = split(splitext(get_tmp_filename())[0])[-1]
        parameters["-t"] = updated_tree_fp
        
        if "-f" not in parameters:
            parameters["-f"] = 'v'
        if "-m" not in parameters:
            parameters["-m"] = 'GTRGAMMA'
        
    elif module=='pplacer':
        try: 
            parameters=param_dict['pplacer']
        except:
            parameters={}
            
        # make sure stats file is passed
        if not opts.stats_fp:
            raise IOError, \
                'When using pplacer, the RAxML produced info file is required.'
                
        # set the primary parameters for pplacer - allow for user-defined
        parameters['--out-dir'] = abspath(output_dir)+'/'
        parameters["-t"] = opts.starting_tree_fp
        parameters['-r'] = opts.refseq_fp
        parameters['-s'] = opts.stats_fp
        
    elif module=='parsinsert':
        try: 
            parameters=param_dict['parsinsert']
        except:
            parameters={}
            
        # define log fp
        log_fp=join(output_dir,'parsinsert.log')

        # define tax assignment values fp
        tax_assign_fp=join(output_dir,'parsinsert_assignments.log')
        parameters["-l"] = log_fp
        parameters["-o"] = tax_assign_fp
        parameters["-s"] = opts.refseq_fp
        parameters["-t"] = opts.starting_tree_fp
    
    # call the module and return a tree object
    result = \
        tree_insertion_module_names[module].insert_sequences_into_tree(seqs, 
                                                moltype=DNA, params=parameters)
    
    result_tree=strip_and_rename_unwanted_labels_from_tree(align_map,result)
    
    # write out the resulting tree
    final_tree=join(output_dir,'%s_final_placement.tre' % (module))
    write_updated_tree_file(final_tree,result)


if __name__ == "__main__":
    main()
