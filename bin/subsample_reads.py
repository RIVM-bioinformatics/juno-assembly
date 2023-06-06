#!/usr/bin/env python3

import subprocess
import tempfile
import logging

def estimate_depth(input: list) -> float:
    """
    Estimate depth of paired-end reads using mash sketch.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        cmd_string = f"mash sketch -o {tmpdir}/tmpfile.msh -k 32 -r -m 3 {input[0]} {input[1]}"
        result = subprocess.run(cmd_string, capture_output=True, shell=True, check=True)
    decoded_result = result.stderr.decode('utf-8')
    genome_size = int(float(decoded_result.split('\n')[0].split(' ')[-1]))
    coverage = float(decoded_result.split('\n')[1].split(' ')[-1])
    logging.info(f"Genome size is estimated to be {genome_size}")
    logging.info(f"Depth of coverage is estimated to be {coverage}")
    return coverage

def calculate_fraction(estimated_depth: float, target_depth: int) -> float:
    """
    Calculate fraction of subsampling.
    """
    subsample_fraction = 1.0
    logging.info(f"Target depth is {args.depth}")
    subsample_fraction = args.depth / estimated_depth
    if subsample_fraction >= 1:
        logging.info(f"Estimated depth is lower than target depth, will not subsample")
    else:
        logging.info(f"Estimated depth is higher than target depth, will subsample")
        logging.info(f"Subsampling using fraction {subsample_fraction}")
    return subsample_fraction

def subsample_reads(input: list, output: list, fraction: float):
    """
    Subsample reads based on calculated fraction.
    """
    if fraction < 1.1: # can overestimate depth a bit
        cmd_string_r1 = f"seqtk seq -f {fraction} -s 1704 {input[0]} > {output[0]}"
        cmd_string_r2 = f"seqtk seq -f {fraction} -s 1704 {input[1]} > {output[1]}"
        logging.info(f"Subsampling files {input} to {output}, resp.")
        subprocess.run(cmd_string_r1, shell=True, check=True)
        subprocess.run(cmd_string_r2, shell=True, check=True)
        logging.info("Finished subsampling reads")
    else:
        cmd_cp_string_r1 = f"cp {input[0]} {output[0]}"
        cmd_cp_string_r2 = f"cp {input[1]} {output[1]}"
        logging.info(f"Copying files {input} to {output}, resp.")
        subprocess.run(cmd_cp_string_r1, shell=True, check=True)
        subprocess.run(cmd_cp_string_r2, shell=True, check=True)
        logging.info(f"Finished copying files")

def main(args):
    coverage = estimate_depth(args.input)
    fraction = calculate_fraction(coverage, args.depth)
    subsample_reads(args.input, args.output, fraction)

    
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input",
                        help="Paired input FASTQ files",
                        type=str,
                        nargs=2,
                        metavar="INPUT_FILE",
                        required=True)
    parser.add_argument("-d", "--depth",
                        help="Target depth [100]",
                        default=100,
                        metavar="INT",
                        type=int)
    parser.add_argument("-o", "--output",
                        help="Paired output FASTQ files",
                        type=str,
                        nargs=2,
                        metavar="OUTPUT_FILE")
    
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO,
                            format='[%(asctime)s] %(message)s',
                            datefmt='%Y/%m/%d %H:%M:%S')
    
    main(args)