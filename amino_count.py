#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
amino_count.py

Description: This program counts the occurrences of standard and non-standard amino acids in a given input file and returns 
two dictionaries: one for counts of standard amino acids and another for counts of non-standard amino acids. The main function 
then processes the results and writes the counts to an output file based on user input or default settings defined in the program.
    
User defined functions:
    - check_files() checks the correctness of the input file's format and existence
    - count_amino_acids() counts the occurrences of standard and non-standard amino acids in a given input file
    
Procedure:
    1. Error handling
        - Check if an output file path is provided
        - Define a custom error class called WrongFormat to handle incorrect file formats
        - Define a function check_files(input_file) to validate the input file
    2. Count amino acids
        - Define a string containing the standard amino acids.
        - Define a function count_amino_acids(input_file) to count amino acids in the input file. This function uses regular 
            expressions to clean sequences and counts both standard and non-standard amino acids.
    3. Main function
        - Call count_amino_acids(input_file) to obtain counts
        - store non-standard amino acids as 'X'
    
Input:
    Amino acid sequences stored in a FASTA file
    An optional output filename 
Output:
    A file containing counts of standard and non-standard amino acids


Usage: python amino_count.py amino.faa output_file.txt[optional]
Date: 12/10/2023
@author: alexaburchak
"""

#%% Error handling 
import sys

input_file = sys.argv[1]

# Warn the user if no output filepath is provided, then create a default output filename
if len(sys.argv) <3:
    output_file = "count_output"
    print("No output filename provided. Output assigned to default 'count_output'.")
else:
    output_file = sys.argv[2]
    
# Build custom error for incorrect file format
class WrongFormat(Exception):
    pass

# Check that the fasta file exists and is in the correct format 
def check_files(input_file):
    try:
        with open(input_file, 'r') as fasta: #check that fasta file is correctly formatted 
            if not '>' in fasta.readline():
                raise WrongFormat
        return True #return True if the fasta file is correctly formattted 
    except FileNotFoundError:  #raise error if the file does not exist 
         print('File does not exist.')
         sys.exit()
    except WrongFormat: #raise error if the file is incorrectly formatted  
         print('File is not in the correct format.')
         sys.exit()
         
#%% Count amino acids
import re  

# Define a string containing the standard amino acids
standard_amino_acids = "ACDEFGHIKLMNPQRSTVWY"

# Define a function to count amino acids in a given input file
def count_amino_acids(input_file):
    # Create a dictionary to store counts of standard amino acids and initialize them to 0
    amino_acid_counts = {aa: 0 for aa in standard_amino_acids}
    
    # Create a dictionary to store counts of non-standard amino acids
    non_standard_counts = {}

    with open(input_file, "r") as infile:
        current_sequence = ""  # Initialize empty string to store the current sequence
        for line in infile:
            if line.startswith(">"):
                if current_sequence:
                    # If there is a current sequence, clean it by removing non-alphabetical characters
                    # and convert it to uppercase
                    current_sequence = re.sub(r'[^a-zA-Z]', '', current_sequence.upper())
                    for aa in current_sequence:
                        if aa in standard_amino_acids:
                            # Update the counts of standard amino acids
                            amino_acid_counts[aa] += 1
                        else:
                            # Update the counts of non-standard amino acids
                            non_standard_counts[aa] = non_standard_counts.get(aa, 0) + 1
                current_sequence = ""  # Reset the current sequence
            else:
                current_sequence += line.strip()  # Append the current line to the current sequence

        # Process the last sequence in the file
        if current_sequence:
            current_sequence = re.sub(r'[^a-zA-Z]', '', current_sequence.upper())
            for aa in current_sequence:
                if aa in standard_amino_acids:
                    amino_acid_counts[aa] += 1
                else:
                    non_standard_counts[aa] = non_standard_counts.get(aa, 0) + 1
                    
    # Calculate the count of non-standard amino acids
    non_standard_counts = sum(non_standard_counts.values())
    
    # Return the counts
    return amino_acid_counts, non_standard_counts  

#%% Apply count_amino_acids to an input file 

def main(input_file, output_file):
    amino_acid_counts, non_standard_counts = count_amino_acids(input_file)

    if output_file:
        with open(output_file, "w") as outfile:
            for amino_acid in standard_amino_acids:
                count = amino_acid_counts[amino_acid]
                # Write the counts of standard amino acids to the output file
                outfile.write(f"{amino_acid}\t{count}\n")
            # Write the count of non-standard amino acids as 'X'
            outfile.write(f"X\t{non_standard_counts}\n")

        print(f"Amino acid abundances written to {output_file}")

    else:
        # If no output file is provided, print the counts of standard amino acids to the console
        for amino_acid in standard_amino_acids:
            count = amino_acid_counts[amino_acid]
            print(f"{amino_acid}\t{count}")
            
        # Print the count of non-standard amino acids as 'X'
        print(f"X\t{non_standard_counts}")

if __name__ == "__main__":
    if check_files(input_file):  # main function will only run if check_files returns True
        main(input_file, output_file)
