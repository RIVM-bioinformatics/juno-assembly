import argparse

import csv
import os
from pathlib import Path
from sys import path
import unittest

main_script_path = str(Path(Path(__file__).parent.absolute()).parent.absolute())
path.insert(0, main_script_path)
from juno_assembly import JunoAssembly


def make_non_empty_file(file_path: Path, num_lines: int = 1000) -> None:
    content = "a\n" * num_lines
    with open(file_path, "w") as file_:
        file_.write(content)


class TestJunoAssemblyDryRun(unittest.TestCase):
    """Testing the junoassembly class (code specific for this pipeline)"""

    @classmethod
    def setUpClass(cls) -> None:
        fake_dirs = ["fake_dir_wsamples", "fake_empty_dir"]

        fake_files = [
            Path("fake_dir_wsamples/sample1_R1.fastq"),
            Path("fake_dir_wsamples/sample1_R2.fastq.gz"),
            Path("fake_dir_wsamples/sample2_R1_filt.fq"),
            Path("fake_dir_wsamples/sample2_R2_filt.fq.gz"),
            Path("fake_dir_wsamples/1234_R1.fastq.gz"),
            Path("fake_dir_wsamples/1234_R2.fastq.gz"),
        ]

        for folder in fake_dirs:
            Path(folder).mkdir(exist_ok=True)
        for file_ in fake_files:
            make_non_empty_file(file_)

        with open("fake_dir_wsamples/fake_metadata.csv", mode="w") as metadata_file:
            metadata_writer = csv.writer(
                metadata_file, delimiter=",", quotechar='"', quoting=csv.QUOTE_MINIMAL
            )
            metadata_writer.writerow(["sample", "genus", "species"])
            metadata_writer.writerow(["sample1", "Salmonella", "enterica"])
            metadata_writer.writerow(["sample2", "Escherichia", "coli"])
            metadata_writer.writerow(["1234", "campylobacter", "jejuni"])

    @classmethod
    def tearDownClass(cls) -> None:
        fake_dirs = ["fake_dir_wsamples", "fake_empty_dir", "test_output"]
        for folder in fake_dirs:
            os.system("rm -rf {}".format(str(folder)))

    def test_fails_with_empty_dir(self) -> None:
        """Testing the pipeline fails if input_dir is empty"""
        with self.assertRaisesRegex(
            ValueError,
            "does not contain any files with the expected format/naming",
        ):
            argv = ["-i", "fake_empty_dir", "-o", "test_output", "-n"]
            pipeline = JunoAssembly(argv=argv)
            pipeline.setup()

    def test_dryrun(self) -> None:
        """Testing the pipeline runs properly as a dry run"""
        full_input_dir = Path("fake_dir_wsamples").resolve()
        pipeline = JunoAssembly(
            argv=["-i", "fake_dir_wsamples", "-o", "test_output", "-n"]
        )
        pipeline.run()
        expected_sample_sheet = {
            "sample1": {
                "R1": str(full_input_dir.joinpath("sample1_R1.fastq")),
                "R2": str(full_input_dir.joinpath("sample1_R2.fastq.gz")),
                "genus": None,
            },
            "sample2": {
                "R1": str(full_input_dir.joinpath("sample2_R1_filt.fq")),
                "R2": str(full_input_dir.joinpath("sample2_R2_filt.fq.gz")),
                "genus": None,
            },
            "1234": {
                "R1": str(full_input_dir.joinpath("1234_R1.fastq.gz")),
                "R2": str(full_input_dir.joinpath("1234_R2.fastq.gz")),
                "genus": None,
            },
        }
        self.assertEqual(pipeline.sample_dict, expected_sample_sheet)

    def test_junoassembly_dryrun_if_genus_provided(self) -> None:
        """Testing the pipeline runs properly as a dry run"""
        pipeline = JunoAssembly(
            argv=[
                "-i",
                "fake_dir_wsamples",
                "-o",
                "test_output",
                "-n",
                "-g",
                "salmonella",
            ]
        )
        full_input_dir = Path("fake_dir_wsamples").resolve()
        pipeline.run()
        expected_sample_sheet = {
            "sample1": {
                "R1": str(full_input_dir.joinpath("sample1_R1.fastq")),
                "R2": str(full_input_dir.joinpath("sample1_R2.fastq.gz")),
                "genus": "salmonella",
            },
            "sample2": {
                "R1": str(full_input_dir.joinpath("sample2_R1_filt.fq")),
                "R2": str(full_input_dir.joinpath("sample2_R2_filt.fq.gz")),
                "genus": "salmonella",
            },
            "1234": {
                "R1": str(full_input_dir.joinpath("1234_R1.fastq.gz")),
                "R2": str(full_input_dir.joinpath("1234_R2.fastq.gz")),
                "genus": "salmonella",
            },
        }
        self.assertEqual(
            pipeline.sample_dict,
            expected_sample_sheet,
            {
                "pipeline": pipeline.sample_dict,
                "expected": expected_sample_sheet,
            },
        )

    def test_junoassembly_dryrun_wMetadata(self) -> None:
        """Testing the pipeline runs properly as a dry run when providing
        a metadata file
        """
        Path("fake_dir_wsamples/missingsamp_1.fastq").unlink()
        Path("fake_dir_wsamples/missingsamp_2.fastq").unlink()

        full_input_dir = Path("fake_dir_wsamples").resolve()
        pipeline = JunoAssembly(
            argv=[
                "-i",
                "fake_dir_wsamples",
                "-o",
                "test_output",
                "-n",
                "-g",
                "salmonella",
                "-m",
                str(full_input_dir.joinpath("fake_metadata.csv")),
            ]
        )
        pipeline.run()
        expected_sample_sheet = {
            "sample1": {
                "R1": str(full_input_dir.joinpath("sample1_R1.fastq")),
                "R2": str(full_input_dir.joinpath("sample1_R2.fastq.gz")),
                "genus": "salmonella",
            },
            "sample2": {
                "R1": str(full_input_dir.joinpath("sample2_R1_filt.fq")),
                "R2": str(full_input_dir.joinpath("sample2_R2_filt.fq.gz")),
                "genus": "escherichia",
            },
            "1234": {
                "R1": str(full_input_dir.joinpath("1234_R1.fastq.gz")),
                "R2": str(full_input_dir.joinpath("1234_R2.fastq.gz")),
                "genus": "campylobacter",
            },
        }
        self.assertEqual(
            pipeline.sample_dict,
            expected_sample_sheet,
            pipeline.sample_dict,
        )

    def test_junoassembly_dryrun_wrong_metadata_colnames(self) -> None:
        """
        Tests whether a good error message is given if the metadata does not have the
        expected column names
        """
        with open("fake_dir_wsamples/fake_metadata2.csv", mode="w") as metadata_file:
            metadata_writer = csv.writer(
                metadata_file, delimiter=",", quotechar='"', quoting=csv.QUOTE_MINIMAL
            )
            metadata_writer.writerow(["Sample", "Genus", "Species"])
            metadata_writer.writerow(["sample1", "Salmonella", "enterica"])
            metadata_writer.writerow(["sample2", "Escherichia", "coli"])
            metadata_writer.writerow(["1234", "campylobacter", "jejuni"])
        full_input_dir = Path("fake_dir_wsamples").resolve()
        with self.assertRaisesRegex(
            AssertionError, "does not contain one or more of the expected column names"
        ):
            pipeline = JunoAssembly(
                argv=[
                    "-i",
                    "fake_dir_wsamples",
                    "-o",
                    "test_output",
                    "-n",
                    "-m",
                    str(full_input_dir.joinpath("fake_metadata2.csv")),
                ]
            )
            pipeline.run()

    def test_metadata_overwrites_genus(self) -> None:
        """Testing the pipeline runs properly as a dry run when providing
        a metadata file. If both a genus and metadata are provided, the
        metadata should overwrite the genus (unless sample not present)
        """
        make_non_empty_file(Path("fake_dir_wsamples/missingsamp_1.fastq"))
        make_non_empty_file(Path("fake_dir_wsamples/missingsamp_2.fastq"))

        input_dir = "fake_dir_wsamples"
        full_input_dir = Path(input_dir).resolve()
        pipeline = JunoAssembly(
            argv=[
                "-i",
                "fake_dir_wsamples",
                "-o",
                "test_output",
                "-n",
                "-g",
                "salmonella",
                "-m",
                str(full_input_dir.joinpath("fake_metadata.csv")),
            ]
        )
        pipeline.run()
        expected_sample_sheet = {
            "sample1": {
                "R1": str(full_input_dir.joinpath("sample1_R1.fastq")),
                "R2": str(full_input_dir.joinpath("sample1_R2.fastq.gz")),
                "genus": "salmonella",
            },
            "sample2": {
                "R1": str(full_input_dir.joinpath("sample2_R1_filt.fq")),
                "R2": str(full_input_dir.joinpath("sample2_R2_filt.fq.gz")),
                "genus": "escherichia",
            },
            "1234": {
                "R1": str(full_input_dir.joinpath("1234_R1.fastq.gz")),
                "R2": str(full_input_dir.joinpath("1234_R2.fastq.gz")),
                "genus": "campylobacter",
            },
            "missingsamp": {
                "R1": str(full_input_dir.joinpath("missingsamp_1.fastq")),
                "R2": str(full_input_dir.joinpath("missingsamp_2.fastq")),
                "genus": "salmonella",
            },
        }
        self.assertEqual(
            pipeline.sample_dict,
            expected_sample_sheet,
            {
                "pipeline": pipeline.sample_dict,
                "expected": expected_sample_sheet,
            },
        )

    def test_junoassembly_dryrun_if_metadata_incomplete(self) -> None:
        """Testing the pipeline runs properly as a dry run when providing a
        metadata file and if a sample is not present in the metadata, then no genus
        is assigned
        """
        make_non_empty_file(Path("fake_dir_wsamples/missingsamp_1.fastq"))
        make_non_empty_file(Path("fake_dir_wsamples/missingsamp_2.fastq"))
        input_dir = "fake_dir_wsamples"
        full_input_dir = Path(input_dir).resolve()

        pipeline = JunoAssembly(
            argv=[
                "-i",
                "fake_dir_wsamples",
                "-o",
                "test_output",
                "-n",
                "-m",
                str(full_input_dir.joinpath("fake_metadata.csv")),
            ]
        )
        pipeline.run()
        expected_sample_sheet = {
            "sample1": {
                "R1": str(full_input_dir.joinpath("sample1_R1.fastq")),
                "R2": str(full_input_dir.joinpath("sample1_R2.fastq.gz")),
                "genus": "salmonella",
            },
            "sample2": {
                "R1": str(full_input_dir.joinpath("sample2_R1_filt.fq")),
                "R2": str(full_input_dir.joinpath("sample2_R2_filt.fq.gz")),
                "genus": "escherichia",
            },
            "1234": {
                "R1": str(full_input_dir.joinpath("1234_R1.fastq.gz")),
                "R2": str(full_input_dir.joinpath("1234_R2.fastq.gz")),
                "genus": "campylobacter",
            },
            "missingsamp": {
                "R1": str(full_input_dir.joinpath("missingsamp_1.fastq")),
                "R2": str(full_input_dir.joinpath("missingsamp_2.fastq")),
                "genus": None,
            },
        }
        self.assertEqual(
            pipeline.sample_dict,
            expected_sample_sheet,
            {
                "pipeline": pipeline.sample_dict,
                "expected": expected_sample_sheet,
            },
        )

    def test_junoassembly_fails_with_unsupported_genus(self) -> None:
        """Testing the pipeline runs properly as a dry run when providing a metadata file"""
        with self.assertRaises(argparse.ArgumentError):
            pipeline = JunoAssembly(
                argv=[
                    "-i",
                    "fake_dir_wsamples",
                    "-o",
                    "test_output",
                    "-g",
                    "fakegenus",
                ]
            )
            pipeline.parser.exit_on_error = False  # type: ignore
            pipeline.setup()


@unittest.skipIf(
    not Path("/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly").exists(),
    "Skipped in non-RIVM environments (because test data is needed)",
)
class TestJunoAssemblyPipeline(unittest.TestCase):
    """Testing the junoassembly class (code specific for this pipeline)"""

    @classmethod
    def setUpClass(cls) -> None:
        os.system("rm -rf test_output")

    @classmethod
    def tearDownClass(cls) -> None:
        os.system("rm -rf test_output")
        os.system("rm -rf test_output_sing")
        os.system("rm -rf test_output_sing_prefix")
        os.system("rm -rf sing_containers")

    def test_junoassembly_run_wMetadata_in_conda(self) -> None:
        """
        Testing the pipeline runs properly with real samples when providing
        a metadata file
        """
        output_dir = Path("test_output")
        input_dir = "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly"
        metadata_file = (
            "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/metadata.csv"
        )
        pipeline = JunoAssembly(
            argv=[
                "-i",
                input_dir,
                "-o",
                str(output_dir),
                "-m",
                metadata_file,
                "--no-containers",
            ]
        )

        pipeline.run()
        expected_sample_sheet = {
            "sample1": {
                "R1": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample1_S14_R1_001.fastq.gz",
                "R2": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample1_S14_R2_001.fastq.gz",
                "genus": "salmonella",
            },
            "sample2": {
                "R1": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample2_R1.fastq.gz",
                "R2": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample2_R2.fastq.gz",
                "genus": "escherichia",
            },
            "sample3": {
                "R1": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample3_R1_001.fastq.gz",
                "R2": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample3_R2_001.fastq.gz",
                "genus": "streptococcus",
            },
            "sample4": {
                "R1": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample4_R1.fastq.gz",
                "R2": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample4_R2.fastq.gz",
                "genus": "shigella",
            },
        }

        self.assertDictEqual(
            pipeline.sample_dict,
            expected_sample_sheet,
            pipeline.sample_dict,
        )
        self.assertTrue(output_dir.joinpath("multiqc", "multiqc.html").exists())
        self.assertTrue(
            output_dir.joinpath(
                "identify_species", "top1_species_multireport.csv"
            ).exists()
        )
        self.assertTrue(output_dir.joinpath("audit_trail", "log_git.yaml").exists())
        self.assertTrue(
            output_dir.joinpath("audit_trail", "log_pipeline.yaml").exists()
        )
        self.assertTrue(output_dir.joinpath("audit_trail", "log_conda.txt").exists())
        self.assertTrue(
            output_dir.joinpath("audit_trail", "juno_assembly_report.html").exists()
        )
        self.assertTrue(
            output_dir.joinpath("audit_trail", "sample_sheet.yaml").exists()
        )
        self.assertTrue(
            output_dir.joinpath("audit_trail", "user_parameters.yaml").exists()
        )

    def test_junoassembly_run_in_singularity(self) -> None:
        """Testing the pipeline runs properly with real samples when providing
        a metadata file
        """
        output_dir = Path("test_output_sing")
        input_dir = "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly"

        pipeline = JunoAssembly(
            argv=[
                "-i",
                input_dir,
                "-o",
                str(output_dir),
                "-m",
            ]
        )
        pipeline.run()

        expected_sample_sheet = {
            "sample1": {
                "R1": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample1_S14_R1_001.fastq.gz",
                "R2": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample1_S14_R2_001.fastq.gz",
                "genus": None,
            },
            "sample2": {
                "R1": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample2_R1.fastq.gz",
                "R2": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample2_R2.fastq.gz",
                "genus": None,
            },
            "sample3": {
                "R1": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample3_R1_001.fastq.gz",
                "R2": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample3_R2_001.fastq.gz",
                "genus": None,
            },
            "sample4": {
                "R1": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample4_R1.fastq.gz",
                "R2": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample4_R2.fastq.gz",
                "genus": None,
            },
        }

        self.assertEqual(
            pipeline.sample_dict,
            expected_sample_sheet,
            pipeline.sample_dict,
        )
        self.assertTrue(output_dir.joinpath("multiqc", "multiqc.html").exists())
        self.assertTrue(
            output_dir.joinpath(
                "identify_species", "top1_species_multireport.csv"
            ).exists()
        )
        self.assertTrue(output_dir.joinpath("audit_trail", "log_git.yaml").exists())
        self.assertTrue(
            output_dir.joinpath("audit_trail", "log_pipeline.yaml").exists()
        )
        self.assertTrue(output_dir.joinpath("audit_trail", "log_conda.txt").exists())
        self.assertTrue(
            output_dir.joinpath("audit_trail", "juno_assembly_report.html").exists()
        )
        self.assertTrue(
            output_dir.joinpath("audit_trail", "sample_sheet.yaml").exists()
        )
        self.assertTrue(
            output_dir.joinpath("audit_trail", "user_parameters.yaml").exists()
        )

    def test_junoassembly_wsingularity_prefix(self) -> None:
        """Testing the pipeline runs properly with real samples when providing
        a metadata file
        """
        output_dir = Path("test_output_sing_prefix")
        input_dir = "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly"

        pipeline = JunoAssembly(
            argv=[
                "-i",
                input_dir,
                "-o",
                str(output_dir),
                "-m",
                "--prefix",
                "sing_containers",
            ]
        )
        pipeline.run()

        expected_sample_sheet = {
            "sample1": {
                "R1": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample1_S14_R1_001.fastq.gz",
                "R2": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample1_S14_R2_001.fastq.gz",
                "genus": None,
            },
            "sample2": {
                "R1": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample2_R1.fastq.gz",
                "R2": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample2_R2.fastq.gz",
                "genus": None,
            },
            "sample3": {
                "R1": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample3_R1_001.fastq.gz",
                "R2": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample3_R2_001.fastq.gz",
                "genus": None,
            },
            "sample4": {
                "R1": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample4_R1.fastq.gz",
                "R2": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample4_R2.fastq.gz",
                "genus": None,
            },
        }

        self.assertEqual(
            pipeline.sample_dict,
            expected_sample_sheet,
            pipeline.sample_dict,
        )
        self.assertTrue(output_dir.joinpath("multiqc", "multiqc.html").exists())
        self.assertTrue(
            output_dir.joinpath(
                "identify_species", "top1_species_multireport.csv"
            ).exists()
        )
        self.assertTrue(output_dir.joinpath("audit_trail", "log_git.yaml").exists())
        self.assertTrue(
            output_dir.joinpath("audit_trail", "log_pipeline.yaml").exists()
        )
        self.assertTrue(output_dir.joinpath("audit_trail", "log_conda.txt").exists())
        self.assertTrue(
            output_dir.joinpath("audit_trail", "juno_assembly_report.html").exists()
        )
        self.assertTrue(
            output_dir.joinpath("audit_trail", "sample_sheet.yaml").exists()
        )
        self.assertTrue(
            output_dir.joinpath("audit_trail", "user_parameters.yaml").exists()
        )

    def test_junoassembly_wdifferent_kmersize(self) -> None:
        """Testing the pipeline runs properly with real samples when providing
        a metadata file
        """
        output_dir = Path("test_output_sing_prefix")
        input_dir = "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly"

        pipeline = JunoAssembly(
            argv=[
                "-i",
                input_dir,
                "-o",
                str(output_dir),
                "-m",
                "--prefix",
                "sing_containers",
                "--kmer-size",
                "21",
                "33",
                "55",
                "77",
            ]
        )
        pipeline.run()


if __name__ == "__main__":
    unittest.main()
