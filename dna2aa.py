#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dna2aa.py

Description: This script accepts a DNA file in fasta format, with one or more sequences, converts T->U, and converts the 
resulting RNA sequence to amino acids using the standard genetic code with stop codons denoted as "*". The resulting protein 
sequence is then saved to an output file. This program also includes error handling for incorrect file formats and provides 
the option to specify an output filename, with a default assigned if none is provided.

User defined functions:
    - check_files() checks the correctness of the input file's format and existence
    - translate_dna() translates a DNA sequence into an amino acid sequence
    
Procedure:
    1. Error handling
        - Check if an output file path is provided
        - Define a custom error class called WrongFormat to handle incorrect file formats
        - Define a function check_files(input_file) to validate the input file
    2. Translate DNA to Amino Acids 
        - Define a standard genetic code dictionary to map codons to amino acids, including stop codons as '*'
        - Define a function translate_dna(dna_sequence) to:
            - translate DNA to RNA by replacing 'T' with 'U'
            - translate RNA to amino acids 
    3. Main translation function
        - Open the input file in read mode and the output file in write mode.
        - Initialize variables to store sequence information.
        - Loop through the input file:
            - If a line starts with '>', parse it to get sequence ID and description.
            - If it's a sequence line, append it to the DNA sequence.
            - Upon encountering a new header, translate the previous DNA sequence, write the result to the output file, and 
                reset the variables.
        - Translate and write the last sequence.
    4. Execution 
        - Check if the input file is correctly formatted using the check_files function.
        - If the check is successful (returns True), call the main function to perform the translation.
        - The translated sequences are saved in the specified output file, and the program concludes.
    
Input:
    A DNA sequence stored in a FASTA file
    An optional output filename 
Output:
    A translated amino acid sequence saved in an output file
    Error messages for incorrect file formats and file absence


Usage: python dna2aa.py DNA.fna output_file.txt[optional]
Date: 12/10/2023
@author: alexaburchak
"""

#%% Error handling 
import sys

input_file = sys.argv[1]

# Warn the user if no output filepath is provided, then create a default output filename
if len(sys.argv) <3:
    output_file = "AA_output"
    print("No output filename provided. Output assigned to default 'AA_output'.")
else:
    output_file = sys.argv[2]


#Build custom error for incorrect file format
class WrongFormat(Exception):
    pass

#Check that the fasta file exists and is in the correct format 
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

#%% Translate DNA to Amino Acid Sequence
import re

# Standard genetic code dictionary, stop codons are denoted by "*"
genetic_code = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}

def translate_dna(dna_sequence):
    # Translate DNA sequence to RNA 
    rna_sequence = dna_sequence.replace("T", "U") # Replace T with U 
    protein_sequence = "" # Initiate an empty protein sequence 

    # Translate the RNA sequence into amino acids
    for i in range(0, len(rna_sequence), 3):
        codon = rna_sequence[i:i + 3] 
        if codon in genetic_code: # add codons to protein sequence by their abbreviations
            protein_sequence += genetic_code[codon]
        else:
            protein_sequence += "X"  # If the codon is invalid, use "X" for unknown amino acid

    return protein_sequence

#%% Apply translate_dna function to provided FASTA file 
def main(input_file, output_file):
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        # Initialize variables to store sequence information
        sequence_id = ""
        description = ""
        dna_sequence = ""
        for line in infile:
            if line.startswith(">"):
                if dna_sequence:
                    # Translate and write the previous sequence to the output file
                    protein_sequence = translate_dna(dna_sequence)
                    outfile.write(f"{sequence_id} {description}\n{protein_sequence}\n")
                
                # Reset for the next sequence
                header = line.strip() # Remove leading/trailing whitespaces from the header
                sequence_id, description = re.split(r'\s', header, 1) # Split the header into sequence_id and description using the first space as a delimiter
                dna_sequence = ""
            else:
                dna_sequence += line.strip() # Remove leading/trailing whitespaces from the sequence lines
        
        # Translate and write the last sequence
        protein_sequence = translate_dna(dna_sequence)
        outfile.write(f"{sequence_id} {description}\n{protein_sequence}\n")

    print("Translation completed and saved to", output_file)


if __name__ == "__main__":
    if check_files(input_file): #main will only run if check_files returns True 
        main(input_file, output_file)
    
