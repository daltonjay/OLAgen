import os
import requests
from PyQt5.QtWidgets import QDialog, QDialogButtonBox, QMessageBox
from PyQt5 import uic

from .output_window import OutputWindow

# Genbank window UI and controller
class GenbankWindow(QDialog):

    def __init__(self, global_state):
        super(GenbankWindow, self).__init__()
        ui_path = os.path.join(os.path.dirname(__file__), '..', 'views', 'genbank_window.ui')
        uic.loadUi(ui_path, self)
        self.show()

        self.global_state = global_state  # Store the global_state instance

        self.genbuttonBox.setStandardButtons(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        self.genbuttonBox.accepted.connect(self.fetch_genbank_data)
        self.genbuttonBox.rejected.connect(self.returnHome)

    def fetch_genbank_data(self):
        reference_accession = self.lineEdit_ref_genbank.text().strip()
        target_accession = self.lineEdit_target_genbank.text().strip()

        if reference_accession and target_accession:
            reference_sequence = self.get_sequence(reference_accession)
            target_sequence = self.get_sequence(target_accession)

            if reference_sequence and target_sequence:
                # Store sequences in global_state or process them as needed
                self.openOutput()
            else:
                QMessageBox.warning(self, "Error", "Failed to fetch sequences from GenBank.")
        else:
            QMessageBox.warning(self, "Error", "Please enter both reference and target GenBank accessions.")

    def get_sequence(self, accession):
        url = f"https://api.ncbi.nlm.nih.gov/datasets/v1/genome/accession/{accession}/sequence"
        response = requests.get(url)
        if response.ok:
            return response.text
        else:
            QMessageBox.warning(self, "Error", f"Failed to fetch sequence for accession {accession}.")
            return None

    def openOutput(self):
        if not hasattr(self, 'output_window'):
            self.output_window = OutputWindow(global_state=self.global_state)
            self.global_state.mainWidget.addWidget(self.output_window)
        self.global_state.mainWidget.setCurrentIndex(self.global_state.mainWidget.indexOf(self.output_window))

    def returnHome(self):
        self.global_state.mainWidget.setCurrentIndex(0)
