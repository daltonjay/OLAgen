import os
from Bio import SeqIO
from PyQt5.QtCore import pyqtSignal, QStringListModel
from PyQt5.QtWidgets import QDialog, QFileDialog, QVBoxLayout, QLabel
from PyQt5 import uic
import sys

# This will get you the directory where fastaWindow.py is located
current_dir = os.path.dirname(os.path.abspath(__file__))

# This will get you the parent directory of current_dir, which is your_project/
project_dir = os.path.dirname(current_dir)

# This checks if project_dir is already in sys.path; if not, it adds it
if project_dir not in sys.path:
    sys.path.insert(0, project_dir)

# Now you can import olagenProcess.py as it is in your_project/
from olagenProcess import *

class fastaWindow(QDialog):
    
    return_global_dict = pyqtSignal(dict)
    return_global_list = pyqtSignal(list)
    
    def __init__(self):
        super(fastaWindow, self).__init__()
        uic.loadUi('GUI/fasta_window.ui', self)
        
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
        main_window = olaGUI()
        widget.addWidget(main_window)
        widget.setCurrentIndex(widget.currentIndex()+1)
    
    def openOutput(self):
        if self.user_fasta_file:
            if (self.clustalCheck.isChecked() or self.mafftCheck.isChecked()):
                output_window = outputWindow()
                widget.addWidget(output_window)
                widget.setCurrentIndex(widget.currentIndex()+1)
        