#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
barcode.py

Description: This program is designed to process a fastq file, separating sequences that contain specific barcodes from those 
that do not. The barcodes are predefined, and the program outputs two separate fastq files: one for sequences with barcodes and 
another for sequences without barcodes. The program also performs error handling, including checking the correctness of the input 
file's format and existence.
    
User defined functions:
    - check_files() checks the correctness of the input file's format and existence
    - contains_barcodes() uses a predefined list of barcodes to extract sequences containing one or more of these barcodes 
    
Procedure:
    1. Error handling
        - Check if an output file path is provided
        - Define a custom error class called WrongFormat to handle incorrect file formats
        - Define a function check_files(input_file) to validate the input file
    2. Determine output filenames
        - If none are provided by the user, default names will be assigned
    3. Define a function to extract barcodes based on predefined list 
    4. Call the main function for sequence processing
        - Initialize two output files for barcodes and no barcodes
        - For each sequence entry (consisting of header, sequence, separator, and quality score lines):
            - If the line starts with "@" (indicating a header line), it assigns this line to current_sequence
            - Read the sequence line and assign it to seq
            - Read the next line (usually a separator line) and assign it to plus
            - Read the quality score line and assign it to current_quality
            - Check if seq contains any of the predefined barcodes using the contains_barcode() function
            - If a barcode is found in seq, write the sequence, removing the barcode, to the output_file. It also includes 
                the header, separator, and quality score
            - If no barcode is found in seq, iwrite the entire sequence entry to the undetermined_file 
    5. Close output files 
    
Input:
    - Expected- input fastq file
    - Optional- output file paths for sequences with barcodes and sequences without barcodes
    
Output:
    - two Fastq output files: one for sequences with barcodes and another for sequences without barcodes
        - file names are set to default values if not specified by the user 



Usage: python barcode.py input_filename output_file1.txt[optional] output_file2.txt[optional]
Date: 12/10/2023
@author: alexaburchak
"""

#%% Error handling 
import sys

#Build custom error for incorrect file format
class WrongFormat(Exception):
    pass

#Check that the fastq file exists and is in the correct format 
def check_files(input_file):
    try:
        with open(input_file, 'r') as fastq: #check that fasta file is correctly formatted 
            if not '@' in fastq.readline():
                raise WrongFormat
        return True #return True if the fasta file is correctly formattted 
    except FileNotFoundError:  #raise error if the file does not exist 
         print('File does not exist.')
         sys.exit()
    except WrongFormat: #raise error if the file is incorrectly formatted  
         print('File is not in the correct format.')
         sys.exit()

#%% Define function to extract a matching barcode from a sequence

# Define the barcode sequences
barcodes = ["TATCCTCT", "GTAAGGAG", "TCTCTCCG"]

def contains_barcode(sequence, barcodes):
    contains_barcode = False # create boolean that is False by default 
    for barcode in barcodes:
        if barcode in sequence[:8] or barcode in sequence[-8:]: # check for barcodes only at the beginning and end of the sequence
            contains_barcode = True # change boolean to True if a barcode is found 
    return contains_barcode

#%% Main function to filter and separate sequences

def main(input_file):
    # Initialize output files for barcodes and no barcodes (undetermined)
    output_file = f"{output_filename}"
    undetermined_file = f"{undetermined_filename}"
    
    output = open(output_file, "w")
    undetermined = open(undetermined_file, "w")
    
    # Process the input fastq file
    with open(input_file, "r") as input_fastq:
        current_sequence = None
        current_quality = None

        for line in input_fastq:
            if line.startswith("@"):  # Header line
                current_sequence = line.strip() # assign header line to current_sequence
                seq = input_fastq.readline() # assign sequence line to seq 
                plus = input_fastq.readline()
                current_quality = input_fastq.readline() # assign quality header to current_quality
                if contains_barcode(seq, barcodes):
                    # Write sequence and quality information to the barcode file
                    output.write(current_sequence + "\n")
                    # Remove barcodes then add sequence to output file 
                    for barcode in barcodes:
                        seq = seq.replace(barcode, "") 
                    output.write(seq)  
                    output.write(plus)
                    output.write(current_quality)
                else:
                    # Write sequence to undetermined file
                    undetermined.write(current_sequence + "\n")
                    undetermined.write(seq)
                    undetermined.write(plus)
                    undetermined.write(current_quality)
                # Reset current_sequence and current_quality 
                current_sequence = None
                current_quality = None

    # Close all output files
    output.close()
    undetermined.close()


input_file = sys.argv[1]
# Define default filenames
output_filename = "sample.fastq" 
undetermined_filename = "undetermined.fastq"

# Warn the user if no output filepaths are provided, then make default output filenames
if len(sys.argv) == 2:
    print("No output filename provided. Output assigned to {} and {}".format(output_filename, undetermined_filename))
# If only one output filepath is provided, check if this is for barcodes or no barcodes 
elif len(sys.argv) == 3:
    x = input("You have only provided one output filename. Is this filepath for sequences with or without barcodes? [with/without]")
    if x == 'with':
        output_filename = sys.argv[2]
    else:
        undetermined_filename = sys.argv[3]
# If two output filepaths are provided then use these names 
elif len(sys.argv) == 4:
    output_filename = sys.argv[2]
    undetermined_filename = sys.argv[3]

# Call the main function to process the input file
main(input_file)
