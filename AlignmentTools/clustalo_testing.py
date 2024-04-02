import subprocess
from Bio import SeqIO
import matplotlib.pyplot as plt

# Define the ClustalO command and arguments
clustalo_cmd = ['/Users/daltonjaynelson/Documents/Research/OLAgen/AlignmentTools/clustalo', '-i', 'seqs_trial.fasta', '-o', 'outputClustalo.fasta', '--outfmt=fasta']

# Execute the ClustalO command
try:
    # Run the ClustalO command
    subprocess.run(clustalo_cmd, check=True)
    print("ClustalO alignment completed successfully.")
except subprocess.CalledProcessError as e:
    # Handle any errors that occur during command execution
    print("Error executing ClustalO:", e)
    
    
fastafiles = list(SeqIO.parse('outputClustalo.fasta', format = 'fasta'))

# Confirm the file loaded properly
for entry in fastafiles:
    print(entry.id)
    
# Make a dictionary for easy access moving forward.
sequences = {}
for entry in fastafiles:
    sequences[entry.id] = entry

def gapped_pos(sequ, pos):
    non_gap = 0
    gaps = 0
    for nt in sequ:
        if nt != '-':
            non_gap += 1
        else:
            gaps += 1
        if non_gap == pos:
            return pos + gaps

spike_start = gapped_pos(sequences['Wuhan_strain'].seq, 21563)
spike_end = gapped_pos(sequences['Wuhan_strain'].seq, 25384)

# Now let's make a dictionary of the spike protein sequences to work with
spikes = {}
for s in fastafiles:
    spikes[s.id] = s.seq[21563-1:25393]
    
# Find the mutations in the spike protein!    
def get_mutations(initial, variant):
    seq_list = list(zip(initial, variant))
    mut_out = []
    for pos, nt in enumerate(seq_list):
        if nt[0] != nt[1]:
            # print(nt[0].upper() + str(pos + 1) + nt[1].upper())
            mut_out.append(nt[0].upper() + str(pos + 1) + nt[1].upper())
    return mut_out
            
            
# ref_sequence = input('What is the name of the reference sequence?: ')
# var_sequence = input('What is the name of the variant sequence?: ')


# Yet, these may be synonymous. Let's look for non-synonymous changes!
# When we aligned originally, our program was unaware. The gaps may cause errors.
# Take spike nt sequences generated, remove dashes, add to a file, translate them
# with biopython's .translate() feature, then re-align them once more.

with open('spikes.fasta', 'w') as f:
    for spike in spikes:
        out = spikes[spike].replace('-', '').translate() # change to AAs!
        f.write('>' + spike + '\n')
        f.write(str(out).upper()+'\n')

# now, align the spikes!
clustalo_cmd = ['/Users/daltonjaynelson/Documents/Research/OLAgen/AlignmentTools/clustalo', '-i', 'spikes.fasta', '-o', 'outputSpikeClustalo.fasta', '--outfmt=fasta']

# Execute the ClustalO command
try:
    # Run the ClustalO command
    subprocess.run(clustalo_cmd, check=True)
    print("ClustalO alignment completed successfully.")
except subprocess.CalledProcessError as e:
    # Handle any errors that occur during command execution
    print("Error executing ClustalO:", e)

spikes_aa = list(SeqIO.parse('outputSpikeClustalo.fasta', format = 'fasta'))

# and create a dictionary with the amino acid sequence we have generated
spike_aa = {}
for entry in spikes_aa:
    spike_aa[entry.id] = entry
    
mutations = {}
for item in spike_aa:
    mutations[item] = get_mutations(spike_aa['Wuhan_strain'].seq, spike_aa[item])

plt.figure(figsize = (15, 5))
for y, item in enumerate(sequences):
    plt.plot((0, len(sequences['Wuhan_strain'])), (y,y), color = 'lightgrey')
    plt.text(-160, y +.35, item, va = 'center', ha = 'left')
    
    for mutation in mutations[item]:
        pos = int(mutation[1:-1])
        aa_change = mutation[-1]
        plt.text(pos, y, aa_change, va = 'center', ha = 'center')
    
    plt.xlim(-175, len(spike_aa['Wuhan_strain']) + 100)
    plt.ylim(-.75, len(sequences) - 0.25)
plt.show()