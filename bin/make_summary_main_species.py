import argparse
import pandas as pd
import pathlib
import sys


class BrackenResult():
    '''
    Class that reads the bracken (species_content.txt) result file per sample 
    and can extract the top species hit
    '''
    def __init__(self,
                filepath):
        self.filepath = pathlib.Path(filepath)
        assert self.filepath.exists(), f"The provided file {filepath} does not exist."

    def update_filepath(self, filepath):
        self.filepath = filepath

    def read_bracken_result(self):
        '''Read the <sample>_species_content.txt produced by Bracken'''
        return pd.read_csv(self.filepath, sep = '\t')

    def find_top_hit(self):
        print(f"Finding top species for {self.filepath}")
        bracken_result = self.read_bracken_result()
        idx_top_species = bracken_result['fraction_total_reads'].argmax()
        return bracken_result.loc[[idx_top_species],['name', 'taxonomy_id', 'fraction_total_reads']]



class BrackenMultireport():
    '''Class that creates the multireport'''
    def __init__(self,
                input_dir = None,
                input_files = [],
                output_multireport = 'top1_species_multireport.csv'):
        if input_dir is not None:
            self.input_dir = pathlib.Path(input_dir)
            assert self.input_dir.exists(), f"The provided input directory {input_dir} does not exist."
        else:
            assert len(input_files) > 0, f"You need to provide either an input_dir or a list of input_files to make a Bracken multireport."
            self.input_dir = None
        self.input_files = input_files
        assert output_multireport.endswith('.csv'), "The output_multireport must have a csv extension (now provided {output_multireport})."
        self.output_multireport = pathlib.Path(output_multireport)

    def enlist_bracken_reports(self):
        '''
        Searches for all the <sample>_species_content.txt files in input_dir
        '''
        assert self.input_dir is not None, f"You need to provide an input dir to find the bracken result files in."
        return [file_ for file_ in self.input_dir.glob('*_species_content.txt')]

    def make_multireport(self):
        '''Read and concatenate the top 'n' result for multiple multireports'''
        print(f"Creating multireport...")
        if len(self.input_files) == 0:
            self.input_files = self.enlist_bracken_reports()
        top1_per_sample = [BrackenResult(file_).find_top_hit() for file_ in self.input_files]
        return pd.concat(top1_per_sample)

    def write_multireport_to_file(self):
        '''Make and write multireport to csv file'''
        report = self.make_multireport()
        print(f"Writing multireport to file {self.output_multireport}...")
        report.to_csv(self.output_multireport, index = False)
        return True

if __name__ == '__main__':
    argument_parser = argparse.ArgumentParser(description='Make bracken results multireport with top 1 species (species with higher score per sample).')
    argument_parser.add_argument('-i', '--input-dir', type=pathlib.Path, 
                        default=None,
                        required=not '-f' in sys.argv and not '--input-files' in sys.argv,
                        help='Directory where to find the individual bracken results (<sample>_species_content.txt files).')
    argument_parser.add_argument('-f', '--input-files', type=pathlib.Path, 
                        default=[], nargs='+',
                        required=not '-i' in sys.argv and not '--input-dir' in sys.argv,
                        help='List of bracken results files (<sample>_species_content.txt files).')
    argument_parser.add_argument('-o', '--output-multireport', type=str, 
                        default='top1_species_multireport.csv',
                        help='Path and name of the output file (must have .csv extension) for the desired multireport.')
    args = argument_parser.parse_args()
    multireport = BrackenMultireport(input_dir = args.input_dir,
                                    input_files = args.input_files,
                                    output_multireport = args.output_multireport)
    multireport.write_multireport_to_file()