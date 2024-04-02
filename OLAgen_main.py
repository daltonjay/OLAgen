# -*- coding: utf-8 -*-
"""
Copyright: Dalton J. Nelson, Vanderbilt University
First Created: April 2023
Last Modified: January 29, 2024
"""
import sys
from PyQt5.QtWidgets import QApplication
from controllers.main_window import MainWindow 
from models.global_state import GlobalState

def main():
    app = QApplication(sys.argv)

    # Initialize the global state
    global_state = GlobalState()

    # Pass the global state to the main window
    main_window = MainWindow(global_state = global_state)
    main_window.show()

    sys.exit(app.exec_())

if __name__ == "__main__":
    main()