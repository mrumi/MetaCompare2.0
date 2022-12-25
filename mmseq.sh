#!/bin/sh

DB='BlastDB/pdb/swissprot'

mmseqs createdb ${1} ${2}sample.contigs
mmseqs taxonomy ${2}sample.contigs $DB ${2}sample.assignments ${2}sample.tmpFolder --tax-lineage 1 --majority 0.7 --vote-mode 1 --lca-mode 3 --orf-filter 1
mmseqs createtsv ${2}sample.contigs ${2}sample.assignments ${2}${3}
rm ${2}sample.contigs*
rm ${2}sample.assignments*
rm -r ${2}sample.tmpFolder
