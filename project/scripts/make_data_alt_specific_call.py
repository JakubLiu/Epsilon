#!/usr/bin/python

import extract_features_alt_specific_call as extf
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--bamlist", type=str)
parser.add_argument("--reference_genome", type=str)
parser.add_argument('--loci_list', type = str)
parser.add_argument('--output_file_prefix', type = str, default = 'output')
parser.add_argument('--windowsize', type = int, default = 100)
args = parser.parse_args()


extf.extract_features_alt_specific(
    bamlist = args.bamlist,
    reference_genome = args.reference_genome,
    loci_list = args.loci_list,
    output_file_prefix = args.output_file_prefix,
    window_size = args.windowsize
)