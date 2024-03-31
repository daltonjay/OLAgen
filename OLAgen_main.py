# -*- coding: utf-8 -*-
"""
Copyright: Dalton J. Nelson, Vanderbilt University
First Created: April 2023
Last Modified: January 29, 2024
"""
import sys
import pandas as pd
from PyQt5.QtWidgets import QApplication
from controllers.main_window import mainWindow

#FIXME: Global Var Relocation
global_muts = None
target_names = None
global_AA_seqs = None
global_SeqIO_seqs = None
global_storage_df = pd.DataFrame(columns = ['Target', 'SNP', 'VP_Probe', 'WT_Probe', 'CP_Probe'])
primer_row_counter = -1        
     
def main():
    app = QApplication(sys.argv)

    mainWindow = mainWindow()
    mainWindow.show()

    sys.exit(app.exec_())

if __name__ == "__main__":
    main()