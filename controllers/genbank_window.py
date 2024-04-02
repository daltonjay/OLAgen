import requests
import subprocess
from PyQt5.QtWidgets import QDialog, QVBoxLayout, QLabel, QDialogButtonBox
from PyQt5.QtCore import pyqtSignal
from PyQt5 import uic
import os
from .output_window import OutputWindow

class GenbankWindow(QDialog):
    return_genbank_data = pyqtSignal(dict)
    switch_view = pyqtSignal(str)  # Signal to indicate view switch

    def __init__(self, global_state):
        super(GenbankWindow, self).__init__()
        ui_path = os.path.join(os.path.dirname(__file__), '..', 'views', 'genbank_window.ui')
        uic.loadUi(ui_path, self)

        self.global_state = global_state  # Store the global_state instance

        self.genbuttonBox.accepted.connect(self.on_ok_clicked)
        self.genbuttonBox.rejected.connect(self.on_cancel_clicked)
        self.genbuttonBox.button(QDialogButtonBox.Ok).clicked.connect(self.fetch_genbank_data)

    def fetch_genbank_data(self):
        # Fetch the sequences from GenBank
        reference_data = self.fetch_sequence_data(self.reference_accession)
        target_data = self.fetch_sequence_data(self.target_accession)
        
        # Save the sequences as files
        reference_file = 'reference_sequence.fasta'
        target_file = 'target_sequence.fasta'
        with open(reference_file, 'w') as file:
            file.write(reference_data)
        with open(target_file, 'w') as file:
            file.write(target_data)

        # Run MAFFT for alignment
        alignment_tool = "mafft" if self.mafft_checkbox.isChecked() else "clustalo"
        if alignment_tool == "mafft":
            alignment_command = ["mafft", "--auto", reference_file]
            process = subprocess.run(alignment_command, input=target_data.encode(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
            alignment_result = process.stdout.decode()

        # Assuming alignment_result is in a format that outputWindow can handle
        self.return_genbank_data.emit({"alignment_result": alignment_result})
        self.on_ok_clicked()

    def fetch_sequence_data(self, accession):
        url = f"https://api.ncbi.nlm.nih.gov/datasets/v1/genome/accession/{accession}/download"
        response = requests.get(url, stream=True)
        if response.ok:
            # Assuming the API returns the sequence data directly
            return response.text
        else:
            raise Exception(f"Failed to fetch data for accession {accession}")

    def returnHome(self):
        self.global_state.mainWidget.setCurrentIndex(0)

    def openOutput(self):
        if not hasattr(self, 'output_window'):
            self.output_window = OutputWindow(global_state=self.global_state)
            self.global_state.mainWidget.addWidget(self.output_window)
        self.global_state.mainWidget.setCurrentIndex(self.global_state.mainWidget.indexOf(self.output_window))