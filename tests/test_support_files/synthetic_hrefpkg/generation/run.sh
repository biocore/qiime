#!/bin/sh
set -eu

for x in a b; do
    (cd in && bppseqgen param=refpkg_${x}.bpp && mv refpkg_${x}.fasta ../out)
    nw_labels in/refpkg_${x}.tre | grep ^q > out/refpkg_${x}.tre.queries
    seqmagick convert --include-from-file out/refpkg_${x}.tre.queries \
        out/refpkg_${x}.fasta out/refpkg_${x}.queries.fasta
    seqmagick mogrify --exclude-from-file out/refpkg_${x}.tre.queries \
        out/refpkg_${x}.fasta
    FastTree -nt -gtr <out/refpkg_${x}.fasta >out/refpkg_${x}.fasttree.tre -log out/refpkg_${x}.fasttree.log
    taxit update ../${x}.refpkg \
        tree=out/refpkg_${x}.fasttree.tre \
        tree_stats=out/refpkg_${x}.fasttree.log \
        aln_fasta=out/refpkg_${x}.fasta
done
cat out/refpkg_[ab].fasta > out/allrefs.fasta
cat out/refpkg_[ab].queries.fasta > out/allqueries.fasta
seqmagick convert out/allrefs.fasta out/allrefs.sto
taxit update ../index.refpkg \
    aln_fasta=out/allrefs.fasta \
    aln_sto=out/allrefs.sto
