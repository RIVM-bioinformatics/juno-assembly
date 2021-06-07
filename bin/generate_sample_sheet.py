"""Generate sample sheet for Juno pipeline.

Usage:
  generate_sample_sheet.py <source_dir>

<source_dir> is a directory containing input fastq files with typical
filenames as used in the legacy (non-automated) process. Output will
be a sample sheet in YAML, for example:

  sample1_id:
    R1:
      path/to/sample_R1.fq.gz
    R2:
      path/to/sample_R2.fq.gz
  sample2_id:
  ...
"""

import argparse
import pathlib
import re
import yaml
import warnings

fq_pattern = re.compile("(.*?)(?:_S\d+_|_S\d+.|_|\.)R?(1|2)(?:_.*\.|\..*\.|\.)f(ast)?q(\.gz)?")

def main(args):
    assert args.dir.is_dir(), "Argument must be a directory."
    
    samples = {}
    small_files = []
    for file_ in args.dir.iterdir():
        if file_.is_dir():
            continue
        if file_.stat().st_size > 3000:
            match = fq_pattern.fullmatch(file_.name)
            if match:
                sample = samples.setdefault(match.group(1), {})
                sample["R{}".format(match.group(2))] = str(file_)
        else:
            small_files.append(str(file_))

    if len(small_files) > 0:   
        warnings.warn("\n\n\033[91mThe following files are too small (<3000 bytes) and were not included in the analysis: \n-{}\033[0m\n\n".format('\n-'.join(small_files)))
    print(yaml.dump(samples, default_flow_style=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("dir", type=pathlib.Path, 
                       help="Directory where input files are located")
    main(parser.parse_args())
