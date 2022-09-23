#!/usr/bin/env nextflow

ch_genes = Channel.from( ['SRR0001', 'SRR0001_mRNA.txt'], ['SRR0002', 'SRR0002_mRNA.txt'] )
                  .view()

ch_mirna = Channel.from( ['SRR0001', 'SRR0001_miRNA.txt'], ['SRR0002', 'SRR0002_miRNA.txt'] )
                  .view()

all_files = ch_genes.join(ch_mirna).view()
