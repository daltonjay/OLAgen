
# -*- coding: utf-8 -*-
"""
Copyright: Dalton J. Nelson, Vanderbilt University
First Created: April 2023
Last Modified: January 29, 2024
"""
import sys
import os
import subprocess
from PyQt5 import QtWidgets
from PyQt5.QtCore import (Qt, pyqtSignal, QStringListModel)
from PyQt5.QtWidgets import *
from PyQt5 import uic
from Bio import SeqIO
from PyQt5.QtWidgets import QWidget
import matplotlib.pyplot as plt

class HelpWindow(QDialog):
    def __init__(self):
        super().__init__()
        
        self.initUI()
        
    def initUI(self):
        self.setWindowTitle('Help Window')
        self.setGeometry(400, 400, 300, 200)
        
        layout = QVBoxLayout()
        
        label = QLabel('This is a help window. Yet, there is no help. Good luck.')
        layout.addWidget(label)
        
        self.setLayout(layout)
        
class genbankWindow(QDialog):
    
    def __init__(self):
        super(genbankWindow, self).__init__()
        uic.loadUi("/Users/daltonjaynelson/Documents/Research/OLAgen/GUI/genbank_window.ui", self)
        self.show()
        
        self.genbuttonBox.setStandardButtons(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        self.genbuttonBox.accepted.connect(self.on_ok_clicked)
        self.genbuttonBox.rejected.connect(self.on_cancel_clicked)
        
    def on_ok_clicked(self):
        output_window = outputWindow()
        widget.addWidget(output_window)
        widget.setCurrentIndex(widget.currentIndex()+1)
        
    def on_cancel_clicked(self):
        main_window = olaGUI()
        widget.addWidget(main_window)
        widget.setCurrentIndex(widget.currentIndex()+1)
        
class fastaWindow(QDialog):
    
    def __init__(self):
        super(fastaWindow, self).__init__()
        uic.loadUi('GUI/fasta_window.ui', self)
        
        self.setWindowTitle('OLAgen - .fasta run')
        self.fileUplBtn.clicked.connect(self.fastaUpload)
        self.mafftCheck.stateChanged.connect(self.checkboxClicked)
        self.clustalCheck.stateChanged.connect(self.checkboxClicked)
         
    def fastaUpload(self):
        options = QFileDialog.Options()
        options |= QFileDialog.ReadOnly 
        
        file_dialog = QFileDialog()
        files, _ = file_dialog.getOpenFileNames(self, 'Select File(s)', '', 'FASTA sequence file(s) (*.fasta or *.fas)', options = options)
        
        if files:
            fasta_file_path = files[0]
            fasta_file_name = self.get_file_name(fasta_file_path)
            
            self.update_fasta_name(fasta_file_name)
            self.update_entryList(fasta_file_path)
            #     self.runOlagen(self.fasta_file_name) 
            
    def get_file_name(self, file_path):
        return os.path.basename(file_path)    
    
    def update_fasta_name(self, text):
        print(text)
        self.userFileLabel.setText("<i>" + text + "</i>")
        
    def update_entryList(self, file_path):  
        fastafiles = list(SeqIO.parse(file_path, format = 'fasta'))
        
        # Confirm the file loaded properly
        elements = [entry.id for entry in fastafiles]
        
        model = self.entryIDList.model()
        if model is None:
            model = QStringListModel()
            self.entryIDList.setModel(model)
        model.setStringList(elements)
        

class outputWindow(QDialog):
    
    def __init__(self):
        super(outputWindow, self).__init__()
        uic.loadUi('GUI/output_window.ui', self)
        
        self.csvExportBtn.clicked.connect(self.promptMainWindow)
        
    def promptMainWindow(self):
        main_window = olaGUI()
        widget.addWidget(main_window)
        widget.setCurrentIndex(widget.currentIndex()+1)

    

class olaGUI(QDialog):
    
    def __init__(self):
        super(olaGUI, self).__init__()
        uic.loadUi("GUI/main_window.ui", self)
        self.show()
        self.setWindowTitle('OLAgen')
        
        self.fasta_window = fastaWindow()
        
        self.fastaInputButton.clicked.connect(self.fastaInit)
        self.helpButton.clicked.connect(self.helpPrompt)
        self.genbankInputButton.clicked.connect(self.genbankPrompt)
        
    def fastaInit(self):       
        fasta_window = fastaWindow()
        widget.addWidget(fasta_window)
        widget.setCurrentIndex(widget.currentIndex()+1)
        

            
    def runOlagen(self, fastasOLA):
        print(fastasOLA)
        fastafiles = list(SeqIO.parse(fastasOLA, format = 'fasta'))

        # Confirm the file loaded properly
        for entry in fastafiles:
            print(entry.id)
            
        # Make a dictionary for easy access moving forward.
        sequences = {}
        for entry in fastafiles:
            sequences[entry.id] = entry
            
        # Use mafft for simple alignment - we first must make sure mafft is installed
        # use this: conda install -c biocore mafft
        # update 1.24.2024: conda install bioconda::mafft

        os.system('mafft ' + fastasOLA +  '> seqs_ali.fasta') # assign aligned sequences to seqs_ali
        # try --leavegappyregion

        # Now let's check to make sure all of the sequences are aligned by looking at 
        # entry lengths.

        ali_files = list(SeqIO.parse('seqs_ali.fasta', format = 'fasta'))

        # Confirm the file loaded properly
        for entry in ali_files:
            print(entry.id)
            print(len(entry.seq))

        # Make a dictionary for easy access moving forward.
        sequences_ali = {}
        for entry in ali_files:
            sequences_ali[entry.id] = entry

        # For now, we are just going to look at the spike protein. 21563...25384
        # But we need to re-align this range with the new alignments generated by mafft

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
        for s in ali_files:
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
        os.system('mafft --genafpair spikes.fasta > spikes_ali.fasta')

        spikes_aa = list(SeqIO.parse('spikes_ali.fasta', format = 'fasta'))

        # and create a dictionary with the amino acid sequence we have generated
        spike_aa = {}
        for entry in spikes_aa:
            spike_aa[entry.id] = entry

                
        # print('And the following are the Amino Acid changes between sequences: ')
        # get_mutations(spike_aa[ref_sequence], spike_aa[var_sequence])

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
            
        target_seq = input('Which sequence is being targeted?: ')
        print('The following are the Amino Acid changes between the target and reference: ')
        print(mutations[target_seq])

        # We would like to compare the target sequence to all other input sequences.
        # There are seemingly a few options: (1) use the generated mutations as comparators
        # or (2) do pairwise comparisons using the get_mutations feature.

        # First, I will attempt (1) by list comparison.

        #for target_mut in target_muts:
          
        def find_target_SNPs(target):
            testlist = mutations[target]
            viable_muts = []
            
            comparator_list = []
            for sequence_element in sequences:
                if sequence_element != target:    
                    for mutation_element in mutations[sequence_element]:
                        comparator_list.append(mutation_element)
                        
            comparator_set = [*set(comparator_list)] # Remove duplicates from the list

            for target_SNP in testlist:
                if target_SNP not in comparator_set:
                    viable_muts.append(target_SNP)
                        
            return viable_muts

        viable_muts = find_target_SNPs(target_seq)
        viableSNPs = []

        for element in viable_muts:
            if '-' not in element:
                viableSNPs.append(element)
                
        print('The following are viable SNPs characteristic to your target for the OLA: ')  
        print(viableSNPs)

        def viableSNP_sequences(target, SNPs):
            '''This is a function that produces the sequences for viable SNP options to be
            targeted in the ligation reaction being designed. It also provides both variable
            ligation probe (VPs) and common ligation probe (CPs) outputs for each SNP.
            In total, three dicts are provided as an output to this function.'''
            
            SNP_sequences = {}
            SNP_CPs = {}
            SNP_VPs = {}
            
            for SNP in SNPs:
                
                # Extract the amino acid location from viable SNPs unique to our target
                num = ''
                location_AA = [num + i for i in SNP if i.isdigit()]
                location_AA = ''.join(location_AA)
                location_AA = int(location_AA)
                
                # Convert the location into a specific base location.
                location_seq = (3 * location_AA) - 3
                
                # Identify the target codon based on the location determined.
                target_codon = spikes[target][location_seq:location_seq+ 3]
                
                # Recognize the exact base change to center the sequence around it for
                # downstream ligation probe generation.
                # Do this by comparing to the first sequence (likely reference)
                seq_names = list(sequences.keys())
                comparator_codon = spikes[seq_names[0]][location_seq:location_seq + 3]
                
                for indx, codon in enumerate(target_codon):
                    if comparator_codon[indx] != codon:
                        base_in_codon = indx
                        
                base_location = base_in_codon + location_seq
                # Variable Probe Sequence
                ligation_hybrid_VP = spikes[target][base_location - 24: base_location + 1]
                # Common Probe Sequence
                ligation_hybrid_CP = spikes[target][base_location + 1: base_location + 27]
                # Entire Hybridization Region
                ligation_hybridization_region = ligation_hybrid_VP + ligation_hybrid_CP
                
                SNP_sequences[SNP] = ligation_hybridization_region
                SNP_CPs[SNP] = ligation_hybrid_CP
                SNP_VPs[SNP] = ligation_hybrid_VP
                
            return SNP_sequences, SNP_VPs, SNP_CPs

        full_region, VP_region, CP_region = viableSNP_sequences(target_seq, viableSNPs)
        print('The targeted ligation probe hybridization regions and their corresponding SNP are: ')
        print(full_region)

        # Now that we've identified the targeted region, we need to know whether the 
        # target sequence corresponds to RNA or DNA. If RNA without PCR, the ligation
        # probes explicitly have to match the cDNA. Otherwise, either works.
        na_type = input('Is your target RNA or DNA?: ')
        if 'r' in na_type.lower():
            pre_rtpcr = input('Are you performing PCR or RT-PCR before the ligation assay? (Y/N): ')

            if 'n' in pre_rtpcr.lower():
                print('If you are performing OLA directly after RT, OLA must match cDNA strand.')
                SNP_user_interest = input('Which of the listed SNPs would you like to target? You may request "all": ')
                
                if SNP_user_interest.lower().strip() == 'all':
                    for snp_seq in full_region:
                        print('SNP: {0} - Rev Complement Seq: {1}'.format(snp_seq, full_region[snp_seq].reverse_complement()))
                elif SNP_user_interest.upper().strip() in full_region:
                    print('SNP: {0} - Rev Complement Seq: {1}'.format(SNP_user_interest, full_region[SNP_user_interest].reverse_complement()))
                    print('Variable Probe: {0}'.format(full_region[SNP_user_interest].reverse_complement()[0:27]))
                    print('Common Probe: /5PHOS/{0}'.format(full_region[SNP_user_interest].reverse_complement()[27:]))
                else:
                    print('Did not register as a possible SNP target.')
           
            elif 'y' in pre_rtpcr.lower():
                SNP_user_interest = input('Which of the listed SNPs would you like to target? You may request "all": ')
                
                if SNP_user_interest.lower().strip() == 'all':
                    for snp_seq in full_region:
                        print('SNP: {0} - Seq: {1}'.format(snp_seq, full_region[snp_seq]))
                elif SNP_user_interest.upper().strip() in full_region:
                    print('SNP: {0} - Seq: {1}'.format(SNP_user_interest, full_region[SNP_user_interest.upper()]))
                    print('Variable Probe: {0}'.format(VP_region[SNP_user_interest.upper()]))
                    print('Common Probe: /5PHOS/{0}'.format(CP_region[SNP_user_interest.upper()]))
                else:
                    print('Did not register as a possible SNP target.')
            
        else:
            SNP_user_interest = input('Which of the listed SNPs do you wish to target? You may request "all": ')
            if SNP_user_interest.lower().strip() == 'all':
                for snp_seq in VP_region:
                    print('Target SNP: {}'.format(snp_seq))
                    print('Variable Probe: {}'.format(VP_region[snp_seq]))
                    print('Common Probe: /5PHOS/{}'.format(CP_region[snp_seq]))
            elif SNP_user_interest.upper().strip() in full_region:
                print('Target SNP: {}'.format(SNP_user_interest.upper()))
                print('Variable Probe: {}'.format(VP_region[SNP_user_interest.upper()]))
                print('Common Probe: /5PHOS/{}'.format(CP_region[SNP_user_interest.upper()]))  


        # Cool, now we have provided the variable and common probes of interest. Next,
        # we should provide a few primers sequence options that don't hybridize with the
        # ligation probes to be used in the full sequence design.
        primer_decision = input('Would you also like some suggested primer regions? (Y/N): ')
        if primer_decision.lower().strip() == 'y':
            multiplex_decision = input('Will your reaction be multiplexed? (Y/N): ')
            if multiplex_decision.lower().strip() == 'y':
                print('Sorry, there is no multiplexing built into this tool just yet! Coming Soon!')
            elif SNP_user_interest.lower().strip() == 'all':
                snp_official = input('Please select one SNP to provide primers for: ')
                official_seq = full_region[snp_official]
                print('These are the potential primer BINDING regions!')
            else:
                official_seq = full_region[SNP_user_interest.upper().strip()]
        else:
            print('Nice! Good luck and thanks for using this tool! Pls dont forget to cite us :)')

        # NOTE 4.20.23: NEED to add a way to not get kicked out of the program if you want to
        # try multiple inputs for a given prompt! 
                    
    
    def helpPrompt(self):
        help_window = HelpWindow()
        help_window.exec_()
        
    def genbankPrompt(self):
        genbank_window = genbankWindow()
        widget.addWidget(genbank_window)
        widget.setCurrentIndex(widget.currentIndex()+1)
       
    def outputPrompt(self):
        output_window = outputWindow()
        widget.addWidget(output_window)
        widget.setCurrentIndex(widget.currentIndex()+1)
        
app = QApplication([])
widget=QtWidgets.QStackedWidget()
main_window = olaGUI()
widget.addWidget(main_window)
widget.setFixedHeight(600)
widget.setFixedWidth(800)
widget.show()
app.exec_()

        
    
    