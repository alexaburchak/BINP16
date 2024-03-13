This repository contains three Python scripts designed for genetic data processing. 

**dna2aa.py**: Translates DNA sequences stored in a FASTA file into amino acid sequences using the standard genetic code. Provides error handling for incorrect file formats.

Usage:
```sh
python dna2aa.py DNA.fna output_file.txt[optional]
```
    
**amino_count.py**: Counts occurrences of standard and non-standard amino acids in a given FASTA file, producing counts for each type and writing them to an output file.

Usage:
```sh
python amino_count.py amino.faa output_file.txt[optional]
```
    
**barcode.py**: Processes a FASTQ file, separating sequences with predefined barcodes from those without, and outputs two separate FASTQ files accordingly. Incorporates error handling for data integrity and allows specification of output filenames.

Usage:
```sh
python barcode.py input_filename output_file1.txt[optional] output_file2.txt[optional]
```
