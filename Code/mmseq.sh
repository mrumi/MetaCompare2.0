#!/bin/sh

DB='/home/lqzhang/mmseq/GTDB/gtdb'
DB='/home/lqzhang/mmseq/pdb/swissprot'

echo ${1}
echo ${2}
echo ${3}



/home/lqzhang/mmseq/MMseqs2/build/bin/mmseqs createdb ${1} ${2}sample.contigs
/home/lqzhang/mmseq/MMseqs2/build/bin/mmseqs taxonomy ${2}sample.contigs $DB ${2}sample.assignments ${2}sample.tmpFolder --tax-lineage 1 --majority 0.7 --vote-mode 1 --lca-mode 3 --orf-filter 1
/home/lqzhang/mmseq/MMseqs2/build/bin/mmseqs createtsv ${2}sample.contigs ${2}sample.assignments ${2}${3}
rm ${2}sample.contigs*
rm ${2}sample.assignments*
rm -r ${2}sample.tmpFolder



