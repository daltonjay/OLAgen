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
        self.setFixedSize(300, 200)

        layout = QVBoxLayout()

        label = QLabel('This is a help window. Yet, there is no help. Good luck.')
        layout.addWidget(label)

        self.setLayout(layout)
