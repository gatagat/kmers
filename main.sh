#!/bin/bash
#
# Command-line interface to count kmers in arbitrary regions.
#
# Counts occurences on both strands, ie. every position is counted twice. All Ns are ignored.
#
# Usage:
#   ./main.sh fasta_file kmers_k bed_file
#
# Example:
#   ./main.sh "/groups/stark/kazmar/data/genomes/dm3.fa" "2" "/groups/stark/kazmar/data/arnold_science_2013/processed/S2_DHSseq_rep1.peaks.txt"
#   AA	1776764
#   CA	1103939
#   GA	896042
#   TA	1139864
#   AC	899529
#   CC	660535
#   GC	915978
#   TC	896042
#   AG	868988
#   CG	738698
#   GG	660535
#   TG	1103939
#   AT	1371644
#   CT	868988
#   GT	899529
#   TT	1776764
#

genome_fasta=$1
kmer_k=$2
bedfile=$3
source ~kazmar/pys/exps/bin/activate
PYTHONPATH=~kazmar/projs:$PYTHONPATH python main.py "$genome_fasta" "$kmer_k" "$bedfile"
