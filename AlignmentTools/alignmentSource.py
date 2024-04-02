"""
Copyright: Dalton J. Nelson, Vanderbilt University
First Created: April 2023
Last Modified: March 28, 2024
"""

# This module contains functions to perform (1) mafft and (2) clustalo alignments
# For use aligning .fasta and GenBank Accession provided sequences as part of OLAgen

import subprocess

def clustAlign(fastaCompiledFile): # Works as desired for a singular input file as of Jan 31, 2024
    
    clustalo_cmd = ['./AlignmentTools/clustalo', '-i', fastaCompiledFile, '-o', 'outClustalo.fasta', '--outfmt=fasta', '--force']
    n = 0 
    # Execute the ClustalO command
    try:
        # Run the ClustalO command
        # Here, add a loading indicator rather than a progress bar
        subprocess.run(clustalo_cmd, check=True)
        print("ClustalO alignment completed successfully.")
        n = 1
        return n
    except subprocess.CalledProcessError as e:
        # Handle any errors that occur during command execution
        print("Error executing ClustalO:", e)
        
#clustAlign('seqs_trial.fasta')

def mafftAlign(fastaCompiledFile): # works as desired for a singular input file as of Jan 31, 2024

    try:
        cmd = ['mafft', '--genafpair', fastaCompiledFile]  # Command to execute
        #cmd = ['mafft', fastaCompiledFile]  # Command to execute
        
        # Open a file for writing the output
        with open('outMafft.fasta', 'w') as output_file:
            # Execute the command and redirect output to the file
            subprocess.run(cmd, stdout=output_file, check=True)
        
        print("Command executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        
# mafftAlign('seqs_trial.fasta')
