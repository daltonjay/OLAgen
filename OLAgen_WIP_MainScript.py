
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
from olagenProcess import runOlagen

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
        self.mafftCheck.stateChanged.connect(self.mafftBoxClicked)
        self.clustalCheck.stateChanged.connect(self.clustalBoxClicked)
        self.alignBtnBox.accepted.connect(self.alignByCheck)
        self.alignBtnBox.rejected.connect(self.returnHome)
        
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
        if self.mafftCheck.isChecked():
            self.statusLabel.setText('')
            print('I will run mafft alignment')
            print(self.user_fasta_file)
            runOlagen(self.user_fasta_file)
        elif self.clustalCheck.isChecked():
            self.statusLabel.setText('')
            print('I will run clustal alignment')
        else:
            self.statusLabel.setText('Please select an alignment method.')
    
    def returnHome(self):
        main_window = olaGUI()
        widget.addWidget(main_window)
        widget.setCurrentIndex(widget.currentIndex()+1)
        
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

        
    
    