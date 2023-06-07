#!/usr/bin/env python3

import subprocess
import tempfile
import logging
from pathlib import Path
import math

def estimate_genome_size(input: list) -> float:
    """
    Estimate genome size of paired-end read set using mash sketch.
    """
    logging.info(f"Estimating genome size for {[str(x) for x in input]}")
    with tempfile.TemporaryDirectory() as tmpdir:
        cmd_string = f"mash sketch -o {tmpdir}/tmpfile.msh -k 21 -r -m 3 {input[0]}"
        result = subprocess.run(cmd_string, capture_output=True, shell=True, check=True)
    decoded_result = result.stderr.decode('utf-8')
    genome_size = int(float(decoded_result.split('\n')[0].split(' ')[-1]))
    # mash depth estimation underestimates actual cov often
    # coverage = float(decoded_result.split('\n')[1].split(' ')[-1]) * 2
    logging.info(f"Genome size is estimated to be {genome_size}")
    return genome_size

def estimate_depth(input: list, genome_size) -> float:
    """
    Estimate depth of paired-end reads using seqtk size and genome size.
    """
    logging.info(f"Estimating depth for {[str(x) for x in input]}")
    with tempfile.TemporaryDirectory() as tmpdir:
        cmd_string = f"seqtk size {input[0]}"
        result = subprocess.run(cmd_string, capture_output=True, shell=True, check=True)
    decoded_result = result.stdout.decode('utf-8')
    total_nt_fw = int(decoded_result.split('\t')[1].rstrip('\n'))
    coverage = (total_nt_fw * 2) / genome_size
    logging.info(f"Depth of coverage is estimated to be {coverage}")
    return coverage

def calculate_fraction(estimated_depth: float, target_depth: int) -> float:
    """
    Calculate fraction of subsampling.
    """
    logging.info(f"Target depth is {target_depth}")
    subsample_fraction = target_depth / estimated_depth
    if subsample_fraction < 1:
        logging.info(f"Estimated depth is higher than target depth, will subsample")
        logging.info(f"Subsampling using fraction {subsample_fraction}")
    else:
        logging.info(f"Estimated depth is lower than or approx. equal to target depth, will not subsample")
    return subsample_fraction

def subsample_reads(input: list, output: list, fraction: float, n_threads: int):
    """
    Subsample reads based on calculated fraction.
    """
    if fraction < 1:
        cmd_string_r1 = f"seqtk seq -f {fraction} -s 1704 {input[0]} | pigz -p {n_threads} > {output[0]}"
        cmd_string_r2 = f"seqtk seq -f {fraction} -s 1704 {input[1]} | pigz -p {n_threads} > {output[1]}"
        logging.info(f"Subsampling files {[str(x) for x in input]} to {[str(x) for x in output]}, resp.")
        subprocess.run(cmd_string_r1, shell=True, check=True)
        subprocess.run(cmd_string_r2, shell=True, check=True)
        logging.info("Finished subsampling reads")
    else:
        cmd_cp_string_r1 = f"cp {input[0]} {output[0]}"
        cmd_cp_string_r2 = f"cp {input[1]} {output[1]}"
        logging.info(f"Copying files {str(input)} to {str(output)}, resp.")
        subprocess.run(cmd_cp_string_r1, shell=True, check=True)
        subprocess.run(cmd_cp_string_r2, shell=True, check=True)
        logging.info(f"Finished copying files")

def calculate_coverage_cutoff(estimated_depth: float, target_depth: int):
    if (target_depth / estimated_depth) < 1:
        logging.info(f"Basing coverage cutoff on target depth ({target_depth})")
        cov_cutoff = math.ceil(target_depth / 5)
    else:
        logging.info(f"Basing coverage cutoff on estimated depth ({estimated_depth})")
        cov_cutoff = math.ceil(estimated_depth / 5)
    logging.info(f"Coverage cutoff set to {cov_cutoff}")
    return cov_cutoff

def output_cov_cutoff(cov_cutoff, cov_cutoff_out):
    with open(cov_cutoff_out, "w") as file:
        file.write(str(cov_cutoff))

def main(args):
    genome_size = estimate_genome_size(args.input)
    coverage = estimate_depth(args.input, genome_size)
    fraction = calculate_fraction(coverage, args.depth)
    subsample_reads(args.input, args.output, fraction, args.threads)
    if args.cov_cutoff_in == "calculate":
        logging.info(f"Coverage cutoff was not specified on command line, will calculate value to use")
        cov_cutoff = calculate_coverage_cutoff(coverage, args.depth)
    else:
        logging.info(f"Coverage cutoff was set to \"{str(args.cov_cutoff_in)}\" on command line, will pass on this value")
        cov_cutoff = args.cov_cutoff_in
    output_cov_cutoff(cov_cutoff, args.cov_cutoff_out)
    
    
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input",
                        help="Paired input FASTQ files",
                        type=Path,
                        nargs=2,
                        metavar="INPUT_FILE",
                        required=True)
    parser.add_argument("-d", "--depth",
                        help="Target depth [100]",
                        default=100,
                        metavar="INT",
                        type=int)
    parser.add_argument("--cov-cutoff-in",
                        help="Input argument for coverage cutoff setting in SPAdes")
    parser.add_argument("--cov-cutoff-out",
                        help="Output file for new coverage cutoff",
                        type=Path,
                        metavar="STR")
    parser.add_argument("-o", "--output",
                        help="Paired output FASTQ files",
                        type=Path,
                        nargs=2,
                        metavar="STR")
    parser.add_argument("-t", "--threads",
                        help="Number of threads to use for pigz [1]",
                        default=1,
                        type=int)
    
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO,
                            format='[%(asctime)s] %(message)s',
                            datefmt='%Y/%m/%d %H:%M:%S')
    
    main(args)