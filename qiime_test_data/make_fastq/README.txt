STARTING FILES (from qiime_tutorial):
    Fasting_Map.txt
    Fasting_Example.fna
    Fasting_Example.qual

COMMANDS:
    split_libraries.py -m $PWD/Fasting_Map.txt -f $PWD/Fasting_Example.fna
        OUTPUTS: seqs.fna
                 histograms.txt
                 split_libraries_log.txt

    make_fastq.py -f $PWD/seqs.fna -q $PWD/Fasting_Example.qual -s
        OUTPUTS: seqs.fna.fastq/PC.354.fastq
                 seqs.fna.fastq/PC.355.fastq
                 seqs.fna.fastq/PC.356.fastq
                 seqs.fna.fastq/PC.481.fastq
                 seqs.fna.fastq/PC.593.fastq
                 seqs.fna.fastq/PC.607.fastq
                 seqs.fna.fastq/PC.634.fastq
                 seqs.fna.fastq/PC.635.fastq
                 seqs.fna.fastq/PC.636.fastq
        (the output directory, seqs.fna.fastq, was renamed to
        seqs.fna.fastq_example to work with script_usage_tests.py)
