#!/usr/bin/env python
# File created on 13 Jul 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from qiime.pycogent_backports.usearch import Usearch, clusters_from_blast_uc_file
from qiime.format import format_otu_map

def _usearch_search(query_fp,
                    refseqs_fp,
                    output_fp,
                    blast6_fp,
                    evalue,
                    min_id,
                    queryalnfract,
                    targetalnfract,
                    maxaccepts,
                    maxrejects,
                    temp_dir,
                    HALT_EXEC=False):
    params = {}
    app = Usearch(params,
                  WorkingDir=temp_dir,
                  HALT_EXEC=HALT_EXEC)
    
    app.Parameters['--evalue'].on(evalue)
    app.Parameters['--id'].on(min_id)
    app.Parameters['--queryalnfract'].on(queryalnfract)
    app.Parameters['--targetalnfract'].on(targetalnfract)
    app.Parameters['--maxaccepts'].on(maxaccepts)
    app.Parameters['--maxrejects'].on(maxrejects)
    
    data = {'--query':query_fp,
            '--uc':output_fp,
            '--db':refseqs_fp,
            '--blast6out':blast6_fp,
            }
    
    app_result = app(data)

def usearch_function_assigner(query_fp,
                              refseqs_fp,
                              output_fp,
                              failure_fp,
                              usearch_fp,
                              blast6_fp,
                              log_fp,
                              evalue,
                              min_id,
                              queryalnfract,
                              targetalnfract,
                              maxaccepts,
                              maxrejects,
                              temp_dir,
                              HALT_EXEC=False):

        _usearch_search(query_fp=query_fp,
                        refseqs_fp=refseqs_fp,
                        output_fp=usearch_fp,
                        blast6_fp=blast6_fp,
                        evalue=evalue,
                        min_id=min_id,
                        queryalnfract=queryalnfract,
                        targetalnfract=targetalnfract,
                        maxaccepts=maxaccepts,
                        maxrejects=maxrejects,
                        temp_dir=temp_dir,
                        HALT_EXEC=HALT_EXEC)

        otus, failures = clusters_from_blast_uc_file(\
         open(usearch_fp,'U'),9)
        function_map_f = open(output_fp,'w')
        for line in format_otu_map(otus.items(),''):
            function_map_f.write(line)
        function_map_f.close()

        ## Write log and failure files - these are place holders right now
        ## while the code is in development
        log_f = open(log_fp,'w')
        log_f.write('place holder...')
        log_f.close()
        
        failure_f = open(failure_fp,'w')
        failure_f.write("## This is not the file you're looking for.\n")
        failure_f.write("## All failed pairwise matches will be listed (hence single ids\n")
        failure_f.write("## being listed multiple times), and many of the ids listed here\n")
        failure_f.write("## will have ultimately been assigned. Generating these anyway as\n")
        failure_f.write("## a placeholder.\n")
        for failure in failures:
            failure_f.write(failure)
            failure_f.write('\n')
        failure_f.close()
