#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 13:41:12 2024

@author: daltonjaynelson
"""
import sys
import os
import subprocess
from PyQt5.QtWidgets import *
from PyQt5 import uic
from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import SeqIO

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
        

class olaGUI(QMainWindow):
    
    def __init__(self):
        super(olaGUI, self).__init__()
        uic.loadUi("2024 Jan - OLAgen GUI v1.ui", self)
        self.show()
        
        self.fastaInputButton.clicked.connect(self.fastaInit)
        self.helpButton.clicked.connect(self.helpPrompt)
        
    def fastaInit(self):
        options = QFileDialog.Options()
        options |= QFileDialog.ReadOnly 
        
        file_dialog = QFileDialog()
        files, _ = file_dialog.getOpenFileNames(self, 'Select File(s)', '', 'FASTA sequence file(s) (*.fasta or *.fas)', options = options)
        
        if files:
            fastafiles = []

            # Parse each file individually
            for file in files:
                fastafiles.extend(list(SeqIO.parse(file, format='fasta')))
            
            # Confirm the files loaded properly
            for entry in fastafiles:
                print(entry.id)
                
            # Make a dictionary for easy access moving forward.
            sequences = {}
            for entry in fastafiles:
                sequences[entry.id] = entry
                
            self.alignFastaFiles(fastafiles)
            
           
            
    def alignFastaFiles(self, fasta_files):
        
        # Use mafft for simple alignment - we first must make sure mafft is installed
        # use this: conda install -c biocore mafft
        
        cmd = ["mafft", "--auto"] + fasta_files + [">", "ola_align.fasta"]
        subprocess.run(cmd)
        
        print(f"Alignment completed. Alignment saved to ola_align.fasta")
            
    def helpPrompt(self):
        help_window = HelpWindow()
        help_window.exec_()
        

def main():
    app = QApplication([])
    window = olaGUI()
    app.exec_()

if __name__ == '__main__':
    main()
        
    
    