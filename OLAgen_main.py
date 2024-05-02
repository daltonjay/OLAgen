"""
Copyright: Dalton J. Nelson, Vanderbilt University 2024
©️ 2024 Dalton J. Nelson and Vanderbilt University. All Rights Reserved. 
Copyrighted material is available free of charge to those individuals who 
have a full-time academic appointment at an academic or other non-profit 
research institution for non-commercial purposes. You may use and make 
copies of the copyrighted material for your research. You may not use the 
copyrighted material for commercial purposes without the written consent 
from the provider. Commercial purposes include sale, lease, license, or 
other transfer of the copyrighted material to a for-profit organization. 
Commercial purposes shall also include uses of the copyrighted material 
or modifications by any organization, to perform contract research, 
to produce or manufacture products for general sale, or to conduct 
research activities that result in any sale, lease, license, or transfer 
of the copyrighted material or modifications to a for-profit organization. 
If you are interested in using the tool for commercial purposes, please 
contact Dalton J. Nelson at dalton.jay.nelson@gmail.com, Rick Haselton 
at rick.haselton@vanderbilt.edu, or Vanderbilt Center for Technology 
Transfer and Commercialization at cttc@vanderbilt.edu.

Contributors: Kunal Chugh
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
