# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'main_windowV2.ui'
##
## Created by: Qt User Interface Compiler version 6.6.1
##
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from PySide6.QtCore import (QCoreApplication, QDate, QDateTime, QLocale,
    QMetaObject, QObject, QPoint, QRect,
    QSize, QTime, QUrl, Qt)
from PySide6.QtGui import (QBrush, QColor, QConicalGradient, QCursor,
    QFont, QFontDatabase, QGradient, QIcon,
    QImage, QKeySequence, QLinearGradient, QPainter,
    QPalette, QPixmap, QRadialGradient, QTransform)
from PySide6.QtWidgets import (QApplication, QDialog, QHBoxLayout, QLabel,
    QPushButton, QSizePolicy, QWidget)

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        if not Dialog.objectName():
            Dialog.setObjectName(u"Dialog")
        Dialog.resize(800, 600)
        self.olaLogo = QLabel(Dialog)
        self.olaLogo.setObjectName(u"olaLogo")
        self.olaLogo.setGeometry(QRect(130, 70, 541, 231))
        self.olaLogo.setPixmap(QPixmap(u"OLAgenLogo.png"))
        self.olaLogo.setScaledContents(True)
        self.horizontalLayoutWidget = QWidget(Dialog)
        self.horizontalLayoutWidget.setObjectName(u"horizontalLayoutWidget")
        self.horizontalLayoutWidget.setGeometry(QRect(170, 300, 461, 41))
        self.horizontalLayout = QHBoxLayout(self.horizontalLayoutWidget)
        self.horizontalLayout.setObjectName(u"horizontalLayout")
        self.horizontalLayout.setContentsMargins(0, 0, 0, 0)
        self.fastaInputButton = QPushButton(self.horizontalLayoutWidget)
        self.fastaInputButton.setObjectName(u"fastaInputButton")
        font = QFont()
        font.setBold(True)
        self.fastaInputButton.setFont(font)
        self.fastaInputButton.setCursor(QCursor(Qt.PointingHandCursor))

        self.horizontalLayout.addWidget(self.fastaInputButton)

        self.helpButton = QPushButton(self.horizontalLayoutWidget)
        self.helpButton.setObjectName(u"helpButton")
        self.helpButton.setFont(font)
        self.helpButton.setCursor(QCursor(Qt.PointingHandCursor))

        self.horizontalLayout.addWidget(self.helpButton)


        self.retranslateUi(Dialog)

        QMetaObject.connectSlotsByName(Dialog)
    # setupUi

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(QCoreApplication.translate("Dialog", u"Dialog", None))
        self.olaLogo.setText("")
        self.fastaInputButton.setText(QCoreApplication.translate("Dialog", u"Upload .fasta", None))
        self.helpButton.setText(QCoreApplication.translate("Dialog", u"Help", None))
    # retranslateUi

