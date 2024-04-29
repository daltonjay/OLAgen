# This module contains functions to perform (1) mafft and (2) clustalo alignments
# For use aligning .fasta provided sequences as part of OLAgen
# Dalton J. Nelson, copyright 2024

import subprocess
import platform

# def clustAlign(fastaCompiledFile): # Works as desired for a singular input file as of Jan 31, 2024
    
#     clustalo_cmd = ['AlignmentTools\\clustalo.exe', '-i', fastaCompiledFile, '-o', 'outClustalo.fasta', '--outfmt=fasta', '--force']

#     #clustalo_cmd = ['./AlignmentTools/clustalo', '-i', fastaCompiledFile, '-o', 'outClustalo.fasta', '--outfmt=fasta', '--force']
#     n = 0 
#     # Execute the ClustalO command
#     try:
#         # Run the ClustalO command
#         # Here, add a loading indicator rather than a progress bar
#         subprocess.run(clustalo_cmd, check=True)
#         print("ClustalO alignment completed successfully.")
#         n = 1
#         return n
#     except subprocess.CalledProcessError as e:
#         # Handle any errors that occur during command execution
#         print("Error executing ClustalO:", e)
def clustAlign(fastaCompiledFile): # Works as desired for a singular input file as of Jan 31, 2024
    system = platform.system()
    if system == 'Darwin' or system == 'Linux':
        # Adjust the path to the clustalo.exe executable as needed
        clustalo_cmd = ['./utils/alignments/clustalo', '-i', fastaCompiledFile, '-o', 'outClustalo.fasta', '--outfmt=fasta', '--force']
       
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
    else:
        print('Unsupported operating system')
        
#clustAlign('seqs_trial.fasta')

def mafftAlign(fastaCompiledFile): # works as desired for a singular input file as of Jan 31, 2024
    system = platform.system()
    print(system)
    
    try:
        # Construct the command
        if system == 'Windows': 
            # Convert Windows path to WSL path
            wsl_path = fastaCompiledFile.replace('\\', '/').replace('C:', '/mnt/c')
            
            cmd = ['wsl.exe', 'bash', '-c', f"mafft --genafpair --out outMafft.fasta {wsl_path}"]
        
        if system == 'Darwin' or system == 'Linux':
            cmd = ['mafft', '--genafpair', fastaCompiledFile]  # Command to execute
        
            with open('outMafft.fasta', 'w') as output_file:
                # Execute the command and redirect output to the file
                subprocess.run(cmd, stdout=output_file, check=True)
        
        # Execute the command
        #subprocess.run(cmd, check=True)
        print("Command executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

    #try:
        #cmd = ['wsl.exe',  'mafft', '--genafpair', fastaCompiledFile]  # Command to execute
        #cmd = ['mafft', fastaCompiledFile]  # Command to execute
        #cmd = ['wsl.exe', 'mafft', '--genafpair', fastaCompiledFile] #forWindows
   
        # Open a file for writing the output
        #with open('outMafft.fasta', 'w') as output_file:
            # Execute the command and redirect output to the file
            #subprocess.run(cmd, stdout=output_file, check=True)
            #subprocess.run(' '.join(cmd), shell=True, check=True)
        
        #print("Command executed successfully.")
    #except subprocess.CalledProcessError as e:
        #print(f"Error executing command: {e}")
    #except Exception as e:
        #print(f"An unexpected error occurred: {e}")
        
# mafftAlign('seqs_trial.fasta')
