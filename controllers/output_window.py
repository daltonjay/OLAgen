''' Copyright 2024, Dalton J. Nelson
    Contributors: Kunal Chugh'''

import os
import sys
import csv
import pandas as pd

from Bio.Seq import Seq
from PyQt5.QtWidgets import QDialog, QMessageBox, QFileDialog, QTableWidget, QTableWidgetItem
from PyQt5 import uic
from PyQt5.QtCore import pyqtSignal

current_dir = os.path.dirname(os.path.abspath(__file__))
project_dir = os.path.dirname(current_dir)
if project_dir not in sys.path:
    sys.path.insert(0, project_dir)
from utils.olagenProcess import *
from utils.dimerScreening import *

# Output window UI and controller
class OutputWindow(QDialog):

    def __init__(self, global_state):
        super(OutputWindow, self).__init__()
        ui_path = os.path.join(os.path.dirname(__file__), '..', 'views', 'output_windowTAB.ui')
        uic.loadUi(ui_path, self)

        self.global_state = global_state  # Store the GlobalState instance
        
        self.csvExportBtn.clicked.connect(self.exportToCSV)
        self.homeButton.clicked.connect(self.returnHome)
        self.targetLst.addItems(self.global_state.target_names)
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
   
    def display_alignment_result(self, alignment_result):
        # Use this method to display the alignment result in your outputWindow
        # For example, if you have a QTextEdit or similar widget to display the result
        self.alignmentResultTextEdit.setText(alignment_result["alignment_result"])

    def update_table(self):
        selected_target = self.targetLst.selectedItems()
        global target_specific_array
        global global_storage_df
        
        if selected_target:
            selected_item = selected_target[0].text()
            values = viability_test(self.global_state.global_muts, selected_item) # get OLA testable SNPs
            full_region , VP_region, CP_region, WT_region = viableSNP_sequences(self.global_state.global_SeqIO_seqs, selected_item, values)

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
                full_region , VP_region, CP_region, WT_region = viableSNP_sequences(self.global_state.global_SeqIO_seqs, selected_item, snp_value, 'reorient')
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
        transposed_array = [[row[i] for row in target_array] for i in range(len(target_array[0]))]
        new_data = pd.DataFrame(transposed_array, columns=['Target', 'SNP', 'VP_Probe', 'WT_Probe', 'CP_Probe'])

        # Update global_storage_df using the global_state instance
        self.global_state.global_storage_df = pd.concat([self.global_state.global_storage_df, new_data])

            
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
     
    def returnHome(self):
        self.global_state.mainWidget.setCurrentIndex(0)
    
