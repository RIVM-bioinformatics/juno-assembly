## Assembly pipeline

(example)Make samplesheet:
```
[verhager@rivm-biohn-l01p assemble_spades]$ source activate default

(default) [verhager@rivm-biohn-l01p assemble_spades]$ python --version
Python 3.6.8 :: Anaconda, Inc.

(default) [verhager@rivm-biohn-l01p assemble_spades]$ python scripts/generate_sample_sheet.py /data/BioGrid/verhager/Gastro/2018_test/2018-11_106934/
'1121800847':
  R1: /data/BioGrid/verhager/Gastro/2018_test/2018-11_106934/1121800847_R1.fastq.gz
  R2: /data/BioGrid/verhager/Gastro/2018_test/2018-11_106934/1121800847_R2.fastq.gz
'1121800852':
  R1: /data/BioGrid/verhager/Gastro/2018_test/2018-11_106934/1121800852_R1.fastq.gz
  R2: /data/BioGrid/verhager/Gastro/2018_test/2018-11_106934/1121800852_R2.fastq.gz
'1121800853':
  R1: /data/BioGrid/verhager/Gastro/2018_test/2018-11_106934/1121800853_R1.fastq.gz
  R2: /data/BioGrid/verhager/Gastro/2018_test/2018-11_106934/1121800853_R2.fastq.gz

(default) [verhager@rivm-biohn-l01p assemble_spades]$ python scripts/generate_sample_sheet.py /data/BioGrid/verhager/Gastro/2018_test/2018-11_106934/ > sample_sheet.yaml

(default) [verhager@rivm-biohn-l01p assemble_spades]$ cat profile/config.yaml| grep sample_sheet.yaml
  - "sample_sheet=sample_sheet.yaml"
```

Dry run:
```
snakemake --profile profile -n
```

Start pipeline:
```
snakemake --profile profile
```
