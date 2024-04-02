# -*- coding: utf-8 -*-
"""
Copyright: Dalton J. Nelson, Vanderbilt University
First Created: April 2023
Last Modified: April 2024
"""
import sys
from PyQt5.QtWidgets import QApplication
from controllers.main_window import MainWindow 
from models.global_state import GlobalState

# Main entry point of the application
def main():
    app = QApplication([])

    global_state = GlobalState() #Global State manages all data (global variables)

    main_widget = global_state.mainWidget #access the mainWidget

    main_window = MainWindow(global_state=global_state)
    main_widget.addWidget(main_window) #add the mainwindow to the widget
    main_widget.setFixedHeight(600)
    main_widget.setFixedWidth(800)
    main_widget.show()

    sys.exit(app.exec_())
    
if __name__ == "__main__":
    main()
