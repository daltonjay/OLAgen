
# -*- coding: utf-8 -*-
"""
Copyright: Dalton J. Nelson, Vanderbilt University
First Created: April 2023
Last Modified: January 29, 2024
"""
import sys
import os
import subprocess
import csv
from PyQt5 import QtWidgets
from PyQt5.QtCore import (Qt, pyqtSignal, QStringListModel)
from PyQt5.QtWidgets import *
from PyQt5 import uic
from Bio import SeqIO
from PyQt5.QtWidgets import QWidget
import matplotlib.pyplot as plt
from olagenProcess import *

global_muts = None
target_names = None
global_AA_seqs = None
global_SeqIO_seqs = None

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
            print("mafft box is checked")
            self.clustalCheck.setChecked(False)
        else:
            print("mafft is UNchecked")
    
    def clustalBoxClicked(self, state):
        if state == 2: # Checked State
            print("clust box is checked")
            self.mafftCheck.setChecked(False)
        else:
            print("clust is UNchecked")
            
    def alignByCheck(self):
        global global_muts
        global global_AA_seqs
        global global_SeqIO_seqs
        
        if self.mafftCheck.isChecked():
            self.statusLabel.setText('')
            print('I will run mafft alignment')
            
            alignedAAs, alignedNTs, seqIO_data = runMafftAlignment(self.user_fasta_file)
            global_AA_seqs = alignedAAs
            mut_out_AAs = self.getMutations(alignedAAs)
            global_muts = mut_out_AAs
            global_SeqIO_seqs = seqIO_data
            #mut_out_NTs = self.getMutations(alignedNTs)
            self.return_global_dict.emit(global_muts)
            
        elif self.clustalCheck.isChecked():
            self.statusLabel.setText('')
            print('I will run clustal alignment')
            
            alignedAAs, alignedNTs, seqIO_data = runClustAlignment(self.user_fasta_file)
            global_AA_seqs = alignedAAs
            mut_out_AAs = self.getMutations(alignedAAs)
            global_muts = mut_out_AAs
            global_SeqIO_seqs = seqIO_data
            #mut_out_NTs = self.getMutations(alignedNTs)
            self.return_global_dict.emit(global_muts)
            
        else:
            self.statusLabel.setText('Please select an alignment method.')
    
    def getMutations(self, pre_aligned_input):
        mutations = {}
        for item in pre_aligned_input:
            mutations[item] = determine_mutations(pre_aligned_input['Wuhan_strain'].seq, pre_aligned_input[item])
        return mutations
    
    def returnHome(self):
        main_window = olaGUI()
        widget.addWidget(main_window)
        widget.setCurrentIndex(widget.currentIndex()+1)
    
    def openOutput(self):
        output_window = outputWindow()
        widget.addWidget(output_window)
        widget.setCurrentIndex(widget.currentIndex()+1)
        
class outputWindow(QDialog):
    
    def __init__(self):
        super(outputWindow, self).__init__()
        uic.loadUi('GUI/output_window.ui', self)
        
        self.csvExportBtn.clicked.connect(self.exportOutputToCSV)
        self.homeButton.clicked.connect(self.promptMainWindow)
        self.targetLst.addItems(target_names)
        self.targetLst.itemSelectionChanged.connect(self.update_table)
        
        self.probeTable.setColumnWidth(0, 61)  # Adjust the width as needed for the first column
        self.probeTable.setColumnWidth(1, 193)  # Adjust the width as needed for the second column
        self.probeTable.setColumnWidth(2, 193)  # Adjust the width as needed for the third column
        self.probeTable.setColumnWidth(3, 193)  # Adjust the width as needed for the fourth column
    
    def update_table(self):
        selected_target = self.targetLst.selectedItems()
        global target_specific_array
        
        if selected_target:
            selected_item = selected_target[0].text()
            values = viability_test(global_muts, selected_item) # get OLA testable SNPs
            full_region , VP_region, CP_region, WT_region = viableSNP_sequences(global_SeqIO_seqs, selected_item, values)

            # Clear existing table
            self.probeTable.clearContents()
            
            # Set the column count
            self.probeTable.setColumnCount(4)
            
            self.probeTable.setColumnWidth(0, 61)  # Adjust the width as needed for the first column
            self.probeTable.setColumnWidth(1, 193)  # Adjust the width as needed for the second column
            self.probeTable.setColumnWidth(2, 193)  # Adjust the width as needed for the third column
            self.probeTable.setColumnWidth(3, 193)  # Adjust the width as needed for the fourth column
            
            # Set the row count to the length of the list
            self.probeTable.setRowCount(len(values))
            
            snpList = []
            vpList = []
            wtList = []
            cpList = []
            # Set the values of each row in the first column
            for i, value in enumerate(values):
                itemSNP = QTableWidgetItem(value)
                itemVP = QTableWidgetItem(str(VP_region[value]))
                itemCP = QTableWidgetItem(str(CP_region[value]))
                itemWT = QTableWidgetItem(str(WT_region[value]))
                self.probeTable.setItem(i, 0, itemSNP)
                self.probeTable.setItem(i, 3, itemCP)
                self.probeTable.setItem(i, 1, itemVP)
                self.probeTable.setItem(i, 2, itemWT)
                snpList.append(value)
                vpList.append(VP_region[value])
                cpList.append(CP_region[value])
                wtList.append(WT_region[value])
            
            target_specific_array = [snpList, vpList, wtList, cpList] 
            
    
    def exportOutputToCSV(self):
        selected_target = self.targetLst.selectedItems()
        
        if len(selected_target):
            transposed_array = [[row[i] for row in target_specific_array] for i in range(len(target_specific_array[0]))]
            
            headers = ['SNP', 'Variable Probe (5\' to 3\')', 'Common Probe (5\' to 3\')']
            
            csv_file_path = selected_target[0].text() + '_OLAgenOutput.csv'
            print(csv_file_path)
            
            with open(csv_file_path, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(headers)
                writer.writerows(transposed_array)
        else:
            self.plsUploadLbl.setText("<i>Please select a target.</i>")
     
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
dialog = fastaWindow()
#dialog.return_global_var.connect(lambda var: print("Returned global variable:", var))
widget=QtWidgets.QStackedWidget()
main_window = olaGUI()
widget.addWidget(main_window)
widget.setFixedHeight(600)
widget.setFixedWidth(800)
widget.show()
app.exec_()