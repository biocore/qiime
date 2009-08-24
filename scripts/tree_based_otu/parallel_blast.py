#!/usr/bin/env python
#parallel_blast.py: make and run parallel blast given file of seqs and db

from os import mkdir, system, path, environ
from os.path import exists
#from split_fasta import split_fasta

def make_blast_dir(dirname):
    if not exists(dirname):
        mkdir(dirname)

def format_blast_database(db_path):
    system('formatdb -p F -i %s' % db_path)

def make_blast_commands(base_infile_path, db_path,  max_job, blast_executible_path,blastmat_path, blast_dir_name):
    s = blast_executible_path + " -n T -p blastn -m 9 -e 1e-30 -b 1 -i " + base_infile_path + ".%s -d " + db_path + " > "+blast_dir_name+"/blast_out.%s"

    for i in range(max_job):
        yield s % (i,i)

def make_blast_jobs(base_infile_path, db_path,  max, blast_executible_path,blastmat_path, blast_dir_name):
    out = open('jobs.txt', 'w')
    for i in make_blast_commands(base_infile_path, db_path, max, blast_executible_path, blastmat_path, blast_dir_name): 
        out.write(i)
        out.write('\n')
    out.close()

def submit_blast_jobs(path_to_cluster_jobs, job_prefix):
    system('%s -ms jobs.txt  %s' % (path_to_cluster_jobs, job_prefix))
    
def split_fasta(input_path, n_seqs, base_split_path):
    tmp_output = []
    file_count = 0
    rec_count = 0
    for i,l in enumerate(input_path):
        if i % 2:
            tmp_output.append(l)
            rec_count += 1
        else:
            tmp_output.append(l)
        
        if rec_count >= n_seqs:
            f = open(base_split_path + '.%d' % file_count, 'w')
            f.write(''.join(tmp_output))
            f.close()
            tmp_output = []
            rec_count = 0
            file_count += 1
    if tmp_output:
        f = open(base_split_path + '.%d' % file_count,'w')
        f.write(''.join(tmp_output))
        f.close()
        file_count += 1
    return file_count
        
        

if __name__ == '__main__':
    from sys import argv
    base_infile_path = '$RCAC_SCRATCH/large_fasttree_datasets/blast_db_inputs/output_badremoved_goodids.fasta'
    db_path = argv[1]
    seqs_per_job = 2500
    blast_dir_name = 'large_fasttree_datasets/blast_results'
    db_file = path.split(argv[1])[-1]
    blast_dir_name = path.join(blast_dir_name, db_file)
    blast_dir_name_abs = path.join(environ['RCAC_SCRATCH'], blast_dir_name)
    blast_dir_name_env = path.join('$RCAC_SCRATCH', blast_dir_name)
    blast_executible_path = '/home/ba01/u116/dmcdonal/blast-2.2.20/bin/blastall'
    blastmat_path = 'export BLASTMAT=$HOME/blast-2.2.20/data'
    job_prefix = "gg_blast"
    cluster_jobs_path = '$HOME/bin/cluster_jobs.py'
    format_blast_database(db_path)
    #max_job = split_fasta(open(base_infile_path), seqs_per_job, base_infile_path)
    make_blast_dir(blast_dir_name_abs)
    make_blast_jobs(base_infile_path, db_path,  1881, blast_executible_path,blastmat_path, blast_dir_name_env)
    #submit_blast_jobs()
