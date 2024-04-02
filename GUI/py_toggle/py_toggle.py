from PySide6.QtCore import *
from PySide6.QtGui import * 
from PySide6.QtWidgets import *

class PyToggle(QCheckBox):
    def __init__(
        self,
        width = 60,
        bg_color = "#777",
        circle_color = "#DDD",
        active_color = "#000BCff"
    ):
        QCheckBox.__init__(self)
        
        # Set Default Paramaters
        self.setFixedSize(width, 28)
        self.setCursor(Qt.PointingHandCursor)