import os
from PyQt5.QtWidgets import QDialog
from PyQt5 import uic
from .fasta_window import FastaWindow
from .help_window import HelpWindow
from .genbank_window import GenbankWindow
from .output_window import OutputWindow

class MainWindow(QDialog):
    def __init__(self, global_state, parent=None):
        super(MainWindow, self).__init__(parent)
        ui_path = os.path.join(os.path.dirname(__file__), '..', 'views', 'main_window.ui')
        uic.loadUi(ui_path, self)
        self.setWindowTitle('OLAgen')

        self.global_state = global_state  # Store the GlobalState instance

        self.fastaInputButton.clicked.connect(self.openFastaWindow)
        self.helpButton.clicked.connect(self.openHelpWindow)
        self.genbankInputButton.clicked.connect(self.openGenbankWindow)

    def openFastaWindow(self):
        if not hasattr(self, 'fasta_window'):
            self.fasta_window = FastaWindow(global_state=self.global_state)
            self.global_state.mainWidget.addWidget(self.fasta_window)
        self.global_state.mainWidget.setCurrentIndex(self.global_state.mainWidget.indexOf(self.fasta_window))

    def openHelpWindow(self):
        if not hasattr(self, 'help_window'):
            self.help_window = HelpWindow()
        self.help_window.exec_()

    def openGenbankWindow(self):
        if not hasattr(self, 'genbank_window'):
            self.genbank_window = GenbankWindow(global_state=self.global_state)
            self.global_state.mainWidget.addWidget(self.genbank_window)
        self.global_state.mainWidget.setCurrentIndex(self.global_state.mainWidget.indexOf(self.genbank_window))

    def openOutputWindow(self):
        if not hasattr(self, 'output_window'):
            self.output_window = OutputWindow(global_state=self.global_state)
            self.global_state.mainWidget.addWidget(self.output_window)
        self.global_state.mainWidget.setCurrentIndex(self.global_state.mainWidget.indexOf(self.output_window))