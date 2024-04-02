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
from PyQt5 import QtWidgets
from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import SeqIO
from qtwidgets import AnimatedToggle # pip install qtpy; pip install qtwidgets


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
        print("Loading UI...")
        self.load_ui()
        print("Initializing widgets...")
        self.init_additional_widgets()
        print("Showing...")
        self.show()
        
        
    def load_ui(self):
        # Load the UI from the .ui file
        ui_file_path = "/Users/daltonjaynelson/Documents/Research/OLAgen/GUI/2024 Jan - OLAgen GUI v1.ui"
        uic.loadUi(ui_file_path, self)
        
    def init_additional_widgets(self):
        # Create the toggle_2 object
        toggle_2 = AnimatedToggle(
            checked_color="#FFB000",
            pulse_checked_color="#44FFB000"
        )
        print('number 1')
        # Create a layout to hold the toggle_2 widget and other UI elements
        layout = QVBoxLayout()
        layout.addWidget(toggle_2)
        print('number 2')
        # Find the central widget container in the main window
        central_widget = self.centralWidget()
        print('number 3')
        # If there's no central widget, create one
        if not central_widget:
            central_widget = QWidget()
            self.setCentralWidget(central_widget)
        print('number 4')
        # Add the loaded UI elements to the layout
        layout.addWidget(central_widget)
        print('number 5')
        # Set the layout as the layout for the central widget
        central_widget.setLayout(layout)
        print('number 6')
        
        
    def fastaInit(self):
        options = QFileDialog.Options()
        options |= QFileDialog.ReadOnly 
        
        file_dialog = QFileDialog()
        files, _ = file_dialog.getOpenFileNames(self, 'Select File(s)', '', 'FASTA sequence file(s) (*.fasta or *.fas)', options = options)
        
        if files:
            fastafiles = list(SeqIO.parse(files, format = 'fasta'))

            # Confirm the file loaded properly
            for entry in fastafiles:
                print(entry.id)
                
            # Make a dictionary for easy access moving forward.
            sequences = {}
            for entry in fastafiles:
                sequences[entry.id] = entry
                
            #self.alignFastaFiles(files)
            
            print(f"Alignment completed. Alignment saved to ola_align.fasta")
            
    def alignFastaFiles(self, fasta_files):
        
        # Use mafft for simple alignment - we first must make sure mafft is installed
        # use this: conda install -c biocore mafft
        
        cmd = f"mafft --auto {fasta_files} > ola_align.fasta"
        subprocess.run(cmd, shell=True)
            
    def helpPrompt(self):
        help_window = HelpWindow()
        help_window.exec_()
        

def main():
    app = QApplication([])
    window = olaGUI()
    app.exec_()
    

if __name__ == '__main__':
    main()
        
    
    