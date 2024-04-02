import os
from PyQt5.QtWidgets import QDialog, QStackedWidget, QVBoxLayout
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

        self.stack = QStackedWidget(self)
        self.layout = QVBoxLayout(self)  # Create a layout
        self.layout.addWidget(self.stack)  # Add the stack to the layout
        self.setLayout(self.layout)  # Set the layout for the main window
        self.setFixedHeight(600)
        self.setFixedWidth(800)

        self.fastaInputButton.clicked.connect(self.openFastaWindow)
        self.helpButton.clicked.connect(self.openHelpWindow)
        self.genbankInputButton.clicked.connect(self.openGenbankWindow)

        self.fasta_window = FastaWindow(global_state=self.global_state)
        self.fasta_window.switch_view.connect(self.handleSwitchView)
        self.stack.addWidget(self.fasta_window)

    def openFastaWindow(self):
        self.stack.setCurrentIndex(self.stack.indexOf(self.fasta_window))

    def openHelpWindow(self):
        if not hasattr(self, 'help_window'):
            self.help_window = HelpWindow()
        self.help_window.exec_()

    def openGenbankWindow(self):
        if not hasattr(self, 'genbank_window'):
            self.genbank_window = GenbankWindow(global_state=self.global_state)
            self.stack.addWidget(self.genbank_window)
        self.stack.setCurrentIndex(self.stack.indexOf(self.genbank_window))

    def openOutputWindow(self):
        if not hasattr(self, 'output_window'):
            self.output_window = OutputWindow(global_state=self.global_state)
            self.stack.addWidget(self.output_window)
        self.stack.setCurrentIndex(self.stack.indexOf(self.output_window))
        
    def handleSwitchView(self, view_name):
        if view_name == 'home':
            self.stack.setCurrentIndex(self.stack.indexOf(self))
        elif view_name == 'output':
            self.openOutputWindow()
