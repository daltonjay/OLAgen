import os
from Bio import SeqIO
from PyQt5.QtCore import pyqtSignal, QStringListModel
from PyQt5.QtWidgets import QDialog, QFileDialog, QVBoxLayout, QLabel
from PyQt5 import uic
import sys

current_dir = os.path.dirname(os.path.abspath(__file__))
project_dir = os.path.dirname(current_dir)
if project_dir not in sys.path:
    sys.path.insert(0, project_dir)

from utils.olagenProcess import *

class FastaWindow(QDialog):
    
    return_global_dict = pyqtSignal(dict)
    return_global_list = pyqtSignal(list)
    switch_view = pyqtSignal(str)  # Signal to indicate view switch
    
    def __init__(self):
        super(FastaWindow, self).__init__()
        ui_path = os.path.join(os.path.dirname(__file__), '..', 'views', 'fasta_window.ui')
        uic.loadUi(ui_path, self)
        
        self.setWindowTitle('OLAgen - .fasta run')
        self.fileUplBtn.clicked.connect(self.fastaUpload)
        self.mafftCheck.stateChanged.connect(self.mafftBoxClicked)
        self.clustalCheck.stateChanged.connect(self.clustalBoxClicked)
        self.alignBtnBox.accepted.connect(self.alignByCheck)
        self.alignBtnBox.rejected.connect(self.returnHome)
        self.alignBtnBox.accepted.connect(self.openOutput)
        
        self.user_fasta_file = ''
         
    def fastaUpload(self):
        options = QFileDialog.Options()
        options |= QFileDialog.ReadOnly 
        
        file_dialog = QFileDialog()
        files, _ = file_dialog.getOpenFileNames(self, 'Select File(s)', '', 'FASTA sequence file(s) (*.fasta or *.fas)', options = options)
        
        if files:
            fasta_file_path = files[0]
            self.user_fasta_file = fasta_file_path
            fasta_file_name = self.get_file_name(fasta_file_path)
            
            self.update_fasta_name(fasta_file_name)
            self.update_entryList(fasta_file_path)
        
            
    def get_file_name(self, file_path):
        return os.path.basename(file_path)    
    
    def update_fasta_name(self, text):
        print(text)
        self.userFileLabel.setText("<i>" + text + "</i>")
        
    def update_entryList(self, file_path):  
        global target_names
        fastafiles = list(SeqIO.parse(file_path, format = 'fasta'))
        
        # Confirm the file loaded properly
        elements = [entry.id for entry in fastafiles]
        target_names = elements
        self.return_global_list.emit(target_names)
        
        model = self.entryIDList.model()
        if model is None:
            model = QStringListModel()
            self.entryIDList.setModel(model)
        model.setStringList(elements)
        
    def mafftBoxClicked(self, state):
        if state == 2: # Checked State
            self.clustalCheck.setChecked(False)
    
    def clustalBoxClicked(self, state):
        if state == 2: # Checked State
            self.mafftCheck.setChecked(False)
            
    def alignByCheck(self):
        global global_muts
        global global_AA_seqs
        global global_SeqIO_seqs
        
        if self.user_fasta_file:
            if self.mafftCheck.isChecked():
                self.statusLabel.setText('')
                
                alignedAAs, alignedNTs, seqIO_data = runMafftAlignment(self.user_fasta_file)
                global_AA_seqs = alignedAAs
                selected_indexes = self.entryIDList.selectedIndexes()
                if selected_indexes:
                    selected_index = selected_indexes[0]
                    reference_sequence = selected_index.data()
                    
                    mut_out_AAs = self.getMutations(alignedAAs, reference_sequence)
                    global_muts = mut_out_AAs
                    global_SeqIO_seqs = seqIO_data
                    #mut_out_NTs = self.getMutations(alignedNTs)
                    self.return_global_dict.emit(global_muts)
                
            elif self.clustalCheck.isChecked():
                self.statusLabel.setText('')
                
                alignedAAs, alignedNTs, seqIO_data = runClustAlignment(self.user_fasta_file)
                global_AA_seqs = alignedAAs
                selected_indexes = self.entryIDList.selectedIndexes()
                if selected_indexes:
                    selected_index = selected_indexes[0]
                    reference_sequence = selected_index.data()
                    
                    mut_out_AAs = self.getMutations(alignedAAs, reference_sequence)
                    global_muts = mut_out_AAs
                    global_SeqIO_seqs = seqIO_data
                    #mut_out_NTs = self.getMutations(alignedNTs)
                    self.return_global_dict.emit(global_muts)
                
            else: 
                self.statusLabel.setText('Please select an alignment method.')
        else:
            self.statusLabel.setText('Please select a .fasta file.')
    
    def getMutations(self, pre_aligned_input, ref_ID):
        mutations = {}
        for item in pre_aligned_input:
            mutations[item] = determine_mutations(pre_aligned_input[ref_ID].seq, pre_aligned_input[item])
        return mutations
    
    def returnHome(self):
        self.switch_view.emit('home')

    def openOutput(self):
        if self.user_fasta_file and (self.clustalCheck.isChecked() or self.mafftCheck.isChecked()):
            self.switch_view.emit('output')
        