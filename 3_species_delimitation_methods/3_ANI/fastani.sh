#!/bin/bash

mkdir -p ani_results

# genome_path_list.txt must contain one genome file path per line

#run FastANI (all-vs-all comparison)
fastANI --ql genome_path_list.txt \
        --rl genome_l_path_ist.txt \
        -o fastani_results.csv
