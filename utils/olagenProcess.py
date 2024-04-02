# Copyright Dalton J. Nelson, 2024

from Bio import SeqIO
from Bio.Seq import Seq

from .alignments.alignmentSource import *
        
def runMafftAlignment(fastasOLA): # here, I will separate the standard alignment -- perhaps call alignmentSource instead
        
        print(fastasOLA)
        
        fastafiles = list(SeqIO.parse(fastasOLA, format = 'fasta'))

        # Confirm the file loaded properly
        for entry in fastafiles:
            print(entry.id)
            
        # Make a dictionary for easy access moving forward.
        sequences = {}
        for entry in fastafiles:
            sequences[entry.id] = entry.seq
            
        mafftAlign(fastasOLA)
        
        fastasAlign = list(SeqIO.parse('outMafft.fasta', format = 'fasta'))
        
        aligned_seq_nts = {}
        aligned_seqEx_nts = {}
        for entry in fastasAlign:
            aligned_seq_nts[entry.id] = entry
            aligned_seqEx_nts[entry.id] = entry.seq
            
        # We also want to have the mutations in amino acids for easy reference
        with open('target_AAs.fasta', 'w') as f:
            for seq in aligned_seqEx_nts:
                out = aligned_seqEx_nts[seq].replace('-', '').translate() # change to AAs!
                f.write('>' + seq + '\n')
                f.write(str(out).upper()+'\n')
                
        mafftAlign('target_AAs.fasta')

        # Now let's check to make sure all of the sequences are aligned by looking at 
        # entry lengths.

        aligned_AAs = list(SeqIO.parse('outMafft.fasta', format = 'fasta'))
        
        # Confirm the file loaded properly
        for entry in aligned_AAs:
            print(entry.id)
            print(len(entry.seq))

        # Make a dictionary for easy access moving forward.
        seq_AAs_aligned = {}
        for entry in aligned_AAs:
            seq_AAs_aligned[entry.id] = entry
        
        return seq_AAs_aligned, aligned_seq_nts, fastasAlign

# alignedAAs, alignedNTs, seqIO_Seqs = runMafftAlignment('spikes_nts.fasta')
def runClustAlignment(fastasOLA): # here, I will separate the standard alignment -- perhaps call alignmentSource instead
        
        print(fastasOLA)
        
        fastafiles = list(SeqIO.parse(fastasOLA, format = 'fasta'))

        # Confirm the file loaded properly
        for entry in fastafiles:
            print(entry.id)
            
        # Make a dictionary for easy access moving forward.
        sequences = {}
        for entry in fastafiles:
            sequences[entry.id] = entry.seq
            
        wait_var = clustAlign(fastasOLA)
    
        if wait_var == 1:
            fastasAlign = list(SeqIO.parse('outClustalo.fasta', format = 'fasta'))
            
            aligned_seq_nts = {}
            aligned_seqEx_nts = {}
            for entry in fastasAlign:
                aligned_seq_nts[entry.id] = entry
                aligned_seqEx_nts[entry.id] = entry.seq
                
            # We also want to have the mutations in amino acids for easy reference
            with open('target_AAs.fasta', 'w') as f:
                for seq in aligned_seqEx_nts:
                    out = aligned_seqEx_nts[seq].replace('-', '').translate() # change to AAs!
                    f.write('>' + seq + '\n')
                    f.write(str(out).upper()+'\n')
                    
            wait_var2 = clustAlign('target_AAs.fasta')
            
            if wait_var2 == 1:
            
                # Now let's check to make sure all of the sequences are aligned by looking at 
                # entry lengths.

                aligned_AAs = list(SeqIO.parse('outClustalo.fasta', format = 'fasta'))
                
                # Confirm the file loaded properly
                for entry in aligned_AAs:
                    print(entry.id)
                    print(len(entry.seq))

                # Make a dictionary for easy access moving forward.
                seq_AAs_aligned = {}
                for entry in aligned_AAs:
                    seq_AAs_aligned[entry.id] = entry
                
                return seq_AAs_aligned, aligned_seqEx_nts, fastasAlign

# aas, nts, fastas = runClustAlignment('spikes_nts.fasta')
# print(aas)
# print(nts)
# print(fastas)

def determine_mutations(initial, variant):
    seq_list = list(zip(initial, variant))
    mut_out = []
    for pos, nt in enumerate(seq_list):
        if nt[0] != nt[1]:
            # print(nt[0].upper() + str(pos + 1) + nt[1].upper())
            mut_out.append(nt[0].upper() + str(pos + 1) + nt[1].upper())
    return mut_out

# mutations = {}
# for item in alignedAAs:
#     mutations[item] = determine_mutations(alignedAAs['Wuhan_strain'].seq, alignedAAs[item])


def viability_test(mutationsDict, targetInterest):
    testList = mutationsDict[targetInterest]
    viableSNPs = []
    
    for element in testList:
        if '-' not in element:
            viableSNPs.append(element)
    return viableSNPs

# viableTry = viability_test(mutations, 'B.1.1.7|Alpha')

def viableSNP_sequences(globalSeqs, target, SNPs, orientation = 'normal'):
    '''This is a function that produces the sequences for viable SNP options to be
    targeted in the ligation reaction being designed. It also provides both variable
    ligation probe (VPs) and common ligation probe (CPs) outputs for each SNP.
    In total, three dicts are provided as an output to this function.'''
    
    base_in_codon = 0
    SNP_sequences = {}
    SNP_CPs = {}
    SNP_VPs = {}
    SNP_WTs = {}
    
    extracted_Seqs = {}
    for entry in globalSeqs:
        extracted_Seqs[entry.id] = entry.seq
    
    for SNP in SNPs:
        
        # Extract the amino acid location from viable SNPs unique to our target
        num = ''
        location_AA = [num + i for i in SNP if i.isdigit()]
        location_AA = ''.join(location_AA)
        location_AA = int(location_AA)
        
        # Convert the location into a specific base location.
        location_seq = (3 * location_AA) - 3
        
        # Identify the target codon based on the location determined.
        target_codon = extracted_Seqs[target][location_seq:location_seq+ 3]
        
        # Recognize the exact base change to center the sequence around it for
        # downstream ligation probe generation.
        # Do this by comparing to the first sequence (likely reference)
        seq_names = list(extracted_Seqs.keys())
        comparator_codon = extracted_Seqs[seq_names[0]][location_seq:location_seq + 3]
        
        for indx, codon in enumerate(target_codon):
            if comparator_codon[indx] != codon:
                base_in_codon = indx
        
        seqBasis = extracted_Seqs[target]
        wtBasis = extracted_Seqs[seq_names[0]] 
        
        base_location = base_in_codon + location_seq
              
        # if orientation == 'normal':
        #     # Variable Probe Sequence
        #     ligation_hybrid_VP = seqBasis[base_location - 24: base_location + 1] 
        #     ligation_WT_VP = wtBasis[base_location - 24: base_location + 1]
        #     # Common Probe Sequence
        #     ligation_hybrid_CP = seqBasis[base_location + 1: base_location + 27]
        # elif orientation == 'reorient':
        #     ligation_hybrid_VP = seqBasis[base_location - 26: base_location].reverse_complement()
        #     ligation_WT_VP = wtBasis[base_location - 26: base_location].reverse_complement()
        #     # Common Probe Sequence
        #     ligation_hybrid_CP = seqBasis[base_location: base_location + 24].reverse_complement()
        #     print(ligation_hybrid_VP)
        # else:
        #     print('provide a valid orientation')
                
        
        
        # Variable Probe Sequence
        ligation_hybrid_VP = extracted_Seqs[target][base_location - 24: base_location + 1] 
        ligation_WT_VP = extracted_Seqs[seq_names[0]][base_location - 24: base_location + 1]
        # Common Probe Sequence
        VP_excess = extracted_Seqs[target][base_location - 26: base_location - 24] # useful when doing reverse orientation
        ligation_hybrid_CP = extracted_Seqs[target][base_location + 1: base_location + 27]
        # Entire Hybridization Region
        ligation_hybridization_region = VP_excess + ligation_hybrid_VP + ligation_hybrid_CP
        
        if orientation == 'reorient':
            ligation_hybridization_region = ligation_hybridization_region.reverse_complement()
            ligation_hybrid_VP = ligation_hybridization_region[2:27]
            ligation_hybrid_CP = ligation_hybridization_region[27:]
            ligation_WT_VP = ligation_hybrid_VP[:-1] + str(Seq(ligation_WT_VP[-1]).reverse_complement())
        
        SNP_sequences[SNP] = ligation_hybridization_region
        SNP_CPs[SNP] = ligation_hybrid_CP
        SNP_VPs[SNP] = ligation_hybrid_VP
        SNP_WTs[SNP] = ligation_WT_VP
        
    return SNP_sequences, SNP_VPs, SNP_CPs, SNP_WTs


# full_region, VP_region, CP_region = viableSNP_sequences(seqIO_Seqs, 'B.1.1.7|Alpha', viableTry)
# print('The targeted ligation probe hybridization regions and their corresponding SNP are: ')
# print(full_region)

