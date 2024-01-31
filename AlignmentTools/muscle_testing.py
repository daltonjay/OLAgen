import subprocess

muscle_exe = r"/Users/daltonjaynelson/Documents/Research/OLAgen/AlignmentTools/muscle"
in_file = r"/Users/daltonjaynelson/Documents/Research/OLAgen/spikes.fasta"
out_file = r"/Users/daltonjaynelson/Documents/Research/OLAgen/outputTesting.fasta"

import subprocess

# Define the Muscle command and arguments
muscle_cmd = ['/Users/daltonjaynelson/Documents/Research/OLAgen/AlignmentTools/muscle', '-in', 'seqs_trial.fasta', '-out', 'output_muscle.fasta']

# Execute the Muscle command
try:
    # Run the Muscle command
    subprocess.run(muscle_cmd, check=True)
    print("Muscle alignment completed successfully.")
except subprocess.CalledProcessError as e:
    # Handle any errors that occur during command execution
    print("Error executing Muscle:", e)
