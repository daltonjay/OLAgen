
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
import pandas as pd
from PyQt5 import QtWidgets
from PyQt5.QtCore import (Qt, pyqtSignal, QStringListModel)
from PyQt5.QtWidgets import *
from PyQt5 import uic
from Bio import SeqIO
from PyQt5.QtWidgets import QWidget
from olagenProcess import *
from dimerScreening import *

global_muts = None
target_names = None
global_AA_seqs = None
global_SeqIO_seqs = None
global_storage_df = pd.DataFrame(columns = ['Target', 'SNP', 'VP_Probe', 'WT_Probe', 'CP_Probe'])
primer_row_counter = -1

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
        
class outputWindow(QDialog):
    
    def __init__(self):
        super(outputWindow, self).__init__()
        uic.loadUi('GUI/output_windowTAB.ui', self)
        
        self.csvExportBtn.clicked.connect(self.exportToCSV)
        self.homeButton.clicked.connect(self.promptMainWindow)
        self.targetLst.addItems(target_names)
        self.targetLst.itemSelectionChanged.connect(self.update_table)
        self.addSOIBtn.clicked.connect(self.addLigationSet)
        self.genPrimerBtn.clicked.connect(self.generate_primers)
        self.genFullReagent.clicked.connect(self.update_all_reagents)
        self.genAltOrient.clicked.connect(self.alternate_orientation)
        
        self.tabWidget.setCurrentIndex(0)
        
        # Ligation Probe Table Settings (Tab Index 0)
        self.probeTable.setColumnWidth(0, 61)  # Adjust the width as needed for the first column
        self.probeTable.setColumnWidth(1, 193)  # Adjust the width as needed for the second column
        self.probeTable.setColumnWidth(2, 193)  # Adjust the width as needed for the third column
        self.probeTable.setColumnWidth(3, 193)  # Adjust the width as needed for the fourth column
        
        # Set of Interest Table Settings
        self.soiTable.setColumnWidth(0, 80)
        self.soiTable.setColumnWidth(1, 193)
        self.soiTable.setColumnWidth(2, 193)
        self.soiTable.setColumnWidth(3, 193)
        
        # Primer Table Settings
        self.primerTable.setColumnWidth(0, 70)
        self.primerTable.setColumnWidth(1, 100)
        self.primerTable.setColumnWidth(2, 70)
        self.primerTable.setColumnWidth(3, 193)
        self.primerTable.setColumnWidth(4, 193)
        self.primerTable.setColumnWidth(5, 193)
        self.primerTable.setColumnWidth(6, 100)
        self.primerTable.setColumnWidth(7, 300)
        
        # All Reagent Table Settings
        self.reagentTable.setColumnWidth(0, 60)
        self.reagentTable.setColumnWidth(1, 100)
        self.reagentTable.setColumnWidth(2, 60)
        self.reagentTable.setColumnWidth(3, 200)
        self.reagentTable.setColumnWidth(4, 200)
        self.reagentTable.setColumnWidth(5, 200)
        self.reagentTable.setColumnWidth(6, 200)
        self.reagentTable.setColumnWidth(7, 200)
        self.reagentTable.setColumnWidth(8, 200)
        self.reagentTable.setColumnWidth(9, 100)
        self.reagentTable.setColumnWidth(10, 300)
        
    
    def update_table(self):
        selected_target = self.targetLst.selectedItems()
        global target_specific_array
        global global_storage_df
        
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
            
            targetList = []
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
                targetList.append(selected_item)
                snpList.append(value)
                vpList.append(str(VP_region[value]))
                cpList.append(str(CP_region[value]))
                wtList.append(str(WT_region[value]))
            
            target_specific_array = [targetList, snpList, vpList, wtList, cpList] 
            self.store_possibilities(target_specific_array)
             
    def alternate_orientation(self):
        selectedSOI = self.probeTable.selectedItems()
        selected_target = self.targetLst.selectedItems()
        
        if selectedSOI:
            if selected_target:
                selected_item = selected_target[0].text()
                reorient_row_data = [item.text() for item in selectedSOI]
                snp_value = [reorient_row_data[0]]
                print(snp_value)

                current_row_count = self.probeTable.rowCount()
                self.probeTable.setRowCount(current_row_count + 1)
                full_region , VP_region, CP_region, WT_region = viableSNP_sequences(global_SeqIO_seqs, selected_item, snp_value, 'reorient')
                print(VP_region, CP_region, WT_region)
                itemSNP = QTableWidgetItem(snp_value[0]+'(-)')
                itemVP = QTableWidgetItem(str(VP_region[snp_value[0]]))
                itemCP = QTableWidgetItem(str(CP_region[snp_value[0]]))
                itemWT = QTableWidgetItem(str(WT_region[snp_value[0]]))
                self.probeTable.setItem(current_row_count, 0, itemSNP)
                self.probeTable.setItem(current_row_count, 3, itemCP)
                self.probeTable.setItem(current_row_count, 1, itemVP)
                self.probeTable.setItem(current_row_count, 2, itemWT)
                
                targetList = [selected_item]
                snpList = [snp_value[0]+'(-)']
                vpList = [str(VP_region[snp_value[0]])]
                wtList = [str(WT_region[snp_value[0]])]
                cpList = [str(CP_region[snp_value[0]])]
                
                target_specific_array = [targetList, snpList, vpList, wtList, cpList]
                self.store_possibilities(target_specific_array)
    
    def store_possibilities(self, target_array):
        global global_storage_df
        
        transposed_array = [[row[i] for row in target_array] for i in range(len(target_array[0]))]
        new_data = pd.DataFrame(transposed_array, columns = ['Target', 'SNP', 'VP_Probe', 'WT_Probe', 'CP_Probe'])
        
        global_storage_df = pd.concat([global_storage_df, new_data])
            
    def addLigationSet(self):
        selectedSOI = self.probeTable.selectedItems()
        if selectedSOI:
            selected_row_data = [item.text() for item in selectedSOI]
            self.add_row_to_SOI(selected_row_data)
            
    def add_row_to_SOI(self, row_data):
        
        selected_target = self.targetLst.selectedItems()
        
        current_row_count = self.soiTable.rowCount()
        self.soiTable.setRowCount(current_row_count + 1)
        
        itemSNP = QTableWidgetItem(str(row_data[0]))
        itemID = QTableWidgetItem('OLAset_' + str(current_row_count + 1))
        itemTarget = QTableWidgetItem(selected_target[0].text())
        itemGen = QTableWidgetItem('No')
        
        self.soiTable.setItem(current_row_count, 0, itemID)
        self.soiTable.setItem(current_row_count, 1, itemTarget)
        self.soiTable.setItem(current_row_count, 2, itemSNP)
        self.soiTable.setItem(current_row_count, 3, itemGen)
        
    def generate_primers(self):
        global primer_row_counter
        global global_storage_df
        
        soi_for_primers = self.soiTable.selectedItems()
        
        if soi_for_primers:
            
            selected_row_data = [item.text() for item in soi_for_primers]
            
            selected_row = soi_for_primers[0].row()
            itemGen = QTableWidgetItem('Yes')
            
            self.soiTable.setItem(selected_row, 3, itemGen) 
        
        choice_target = str(selected_row_data[1])
        choice_snp = str(selected_row_data[2])
        
        desired_row = global_storage_df[(global_storage_df['Target'] == choice_target) & (global_storage_df['SNP'] == choice_snp)]
        
        success_primers_df = primer_library_test(primer_df, str(desired_row['VP_Probe']), str(desired_row['CP_Probe']))
        print(success_primers_df)
        
        
        for row_indx, row in success_primers_df.iterrows():
            
            current_row_count = self.primerTable.rowCount()
            self.primerTable.setRowCount(current_row_count + 1)
            primer_row_counter += 1
            
            itemID = QTableWidgetItem(str(selected_row_data[0]))
            targItem = QTableWidgetItem(str(selected_row_data[1]))
            snpItem = QTableWidgetItem(str(selected_row_data[2]))
            fwdItem = QTableWidgetItem(str(row.iloc[0]))
            revItem = QTableWidgetItem(str(row.iloc[1]))
            hydrItem = QTableWidgetItem(str(row.iloc[2]))
            typeItem = QTableWidgetItem(str(row.iloc[3]))
            citationItem = QTableWidgetItem(str(row.iloc[4]))
            
            self.primerTable.setItem(primer_row_counter, 0, itemID)
            self.primerTable.setItem(primer_row_counter, 1, targItem)
            self.primerTable.setItem(primer_row_counter, 2, snpItem)
            self.primerTable.setItem(primer_row_counter, 3, fwdItem)
            self.primerTable.setItem(primer_row_counter, 4, revItem)
            self.primerTable.setItem(primer_row_counter, 5, hydrItem)
            self.primerTable.setItem(primer_row_counter, 6, typeItem)
            self.primerTable.setItem(primer_row_counter, 7, citationItem)
    
    def update_all_reagents(self):
        global global_storage_df
        print(global_storage_df.head())
        
        soi_for_reagents = self.primerTable.selectedItems() 
        
        if soi_for_reagents:
            selected_row_data = [item.text() for item in soi_for_reagents]  
            
        choice_ID = str(selected_row_data[0])
        choice_target = str(selected_row_data[1])
        choice_snp = str(selected_row_data[2])
        choice_FWD = str(selected_row_data[3])
        choice_REV = str(selected_row_data[4])
        choice_HYD = str(selected_row_data[5])
        choice_Type = str(selected_row_data[6])
        choice_cite = str(selected_row_data[7])
        choice_revRC = str(Seq(choice_REV).reverse_complement())
        
        desired_row = global_storage_df[(global_storage_df['Target'] == choice_target) & (global_storage_df['SNP'] == choice_snp)]
        print(desired_row['VP_Probe'].values[0])
        
        idItem = QTableWidgetItem(choice_ID)
        targItem = QTableWidgetItem(choice_target)
        snpItem = QTableWidgetItem(choice_snp)
        vpItem = QTableWidgetItem(choice_FWD + 'cgc' + choice_HYD + desired_row['VP_Probe'].values[0])
        vpWtItem = QTableWidgetItem(choice_FWD + 'cgc' + desired_row['WT_Probe'].values[0])
        cpItem = QTableWidgetItem('/5Phos/' + desired_row['CP_Probe'].values[0] + choice_revRC)
        fwdItem = QTableWidgetItem(choice_FWD)
        revItem = QTableWidgetItem(choice_REV)
        hydrItem = QTableWidgetItem(choice_HYD)
        typeItem = QTableWidgetItem(choice_Type)
        citeItem = QTableWidgetItem(choice_cite)
        
        
        current_row_count = self.reagentTable.rowCount()
        self.reagentTable.setRowCount(current_row_count + 1)
        
        self.reagentTable.setItem(current_row_count, 0, idItem)
        self.reagentTable.setItem(current_row_count, 1, targItem)
        self.reagentTable.setItem(current_row_count, 2, snpItem)
        self.reagentTable.setItem(current_row_count, 3, vpItem)
        self.reagentTable.setItem(current_row_count, 4, vpWtItem)
        self.reagentTable.setItem(current_row_count, 5, cpItem)
        self.reagentTable.setItem(current_row_count, 6, fwdItem)
        self.reagentTable.setItem(current_row_count, 7, revItem)
        self.reagentTable.setItem(current_row_count, 8, hydrItem)
        self.reagentTable.setItem(current_row_count, 9, typeItem)
        self.reagentTable.setItem(current_row_count, 10, citeItem)
        
            
    def exportToCSV(self):
        filename, _ = QFileDialog.getSaveFileName(self, "Export Data", "", "CSV Files (*.csv)")

        if filename:
            try:
                # Get the current tab index
                currentTabIndex = self.tabWidget.currentIndex()
                
                # Get the table widget associated with the current tab
                currentTableWidget = self.tabWidget.widget(currentTabIndex).findChild(QTableWidget)

                with open(filename, "w", newline='') as csvfile:
                    writer = csv.writer(csvfile)
                    
                    # Write headers
                    headers = []
                    for column in range(currentTableWidget.columnCount()):
                        headers.append(currentTableWidget.horizontalHeaderItem(column).text())
                    writer.writerow(headers)

                    # Write data
                    for row in range(currentTableWidget.rowCount()):
                        row_data = []
                        for column in range(currentTableWidget.columnCount()):
                            item = currentTableWidget.item(row, column)
                            if item is not None:
                                row_data.append(item.text())
                            else:
                                row_data.append('')
                        writer.writerow(row_data)
                QMessageBox.information(self, "Success", "Data exported successfully!")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Error exporting data: {str(e)}")

     
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
widget = QtWidgets.QStackedWidget()
main_window = olaGUI()
widget.addWidget(main_window)
widget.setFixedHeight(600)
widget.setFixedWidth(800)
widget.setWindowTitle("OLAgen")
widget.show()
sys.exit(app.exec_())