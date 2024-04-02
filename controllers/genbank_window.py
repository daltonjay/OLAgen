from PyQt5.QtWidgets import QDialog, QDialogButtonBox
from PyQt5 import uic
import os
from .output_window import OutputWindow

class GenbankWindow(QDialog):

    def __init__(self, global_state):
        super(GenbankWindow, self).__init__()
        ui_path = os.path.join(os.path.dirname(__file__), '..', 'views', 'genbank_window.ui')
        uic.loadUi(ui_path, self)
        self.show()

        self.global_state = global_state  # Store the global_state instance

        self.genbuttonBox.setStandardButtons(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        self.genbuttonBox.accepted.connect(self.openOutput)
        self.genbuttonBox.rejected.connect(self.returnHome)

    def openOutput(self):
        if not hasattr(self, 'output_window'):
            self.output_window = OutputWindow(global_state=self.global_state)
            self.global_state.mainWidget.addWidget(self.output_window)
        self.global_state.mainWidget.setCurrentIndex(self.global_state.mainWidget.indexOf(self.output_window))

    def returnHome(self):
        self.global_state.mainWidget.setCurrentIndex(0)
