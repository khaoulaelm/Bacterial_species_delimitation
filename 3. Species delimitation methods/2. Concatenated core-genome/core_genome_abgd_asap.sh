#!/bin/bash

#ABGD
./abgd -a -o ~/abgd_core_output/ -d JC69 ~/roary/core_gene_alignment.aln

#ASAP
./asap -a -o ~/asap_core_output ~/roary/core_gene_alignment.aln
