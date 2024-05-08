''' Copyright 2024, Dalton J. Nelson
    Contributors: Kunal Chugh'''

from PyQt5.QtWidgets import QDialog, QVBoxLayout, QLabel

# Help dialog UI and controller
class HelpWindow(QDialog):
    def __init__(self):
        super().__init__()
        
        self.initUI()
        
    def initUI(self):
        self.setWindowTitle('Help Window')
        self.setFixedSize(500, 550)

        layout = QVBoxLayout()

        label = QLabel(
        '''
        OLAgen is a tool for designing oligonucleotide ligation assay (OLA) 
        and OLA-PCR reagent sequences. 
        
        In depth details will be covered in publication (In Preparation).
        
        A tutorial is provided at the YouTube channel @daltonjaynelson. The
        tutorial covers workflow, file setup and preparation (FASTA), and 
        feature descriptions. 
        
        For initial troubleshooting, please ensure that you have the following 
        dependencies installed:
        
        - Biopython (1.83)
        - numpy (1.26.4)
        - pandas (2.2.2)
        - openpyxl (3.1.2)
        - PyQt5 (5.15.10)
        - Primer3-py (2.0.3)
        
        You must also have operating system specific installations of MAFFT
        and/or Clustal Omega.
        
        Further troubleshooting or feedback, please see github.com/daltonjay 
        or contact dalton.jay.nelson@gmail.com.
        '''
        )
        layout.addWidget(label)

        self.setLayout(layout)
