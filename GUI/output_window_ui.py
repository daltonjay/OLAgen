# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'output_window.ui'
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
from PySide6.QtWidgets import (QAbstractItemView, QApplication, QDialog, QHeaderView,
    QLabel, QListWidget, QListWidgetItem, QPushButton,
    QSizePolicy, QTableWidget, QTableWidgetItem, QVBoxLayout,
    QWidget)

class Ui_outputDialog(object):
    def setupUi(self, outputDialog):
        if not outputDialog.objectName():
            outputDialog.setObjectName(u"outputDialog")
        outputDialog.resize(800, 600)
        self.logoLabel = QLabel(outputDialog)
        self.logoLabel.setObjectName(u"logoLabel")
        self.logoLabel.setGeometry(QRect(338, 30, 151, 61))
        self.logoLabel.setPixmap(QPixmap(u"OLAgenLogo.png"))
        self.logoLabel.setScaledContents(True)
        self.probeTable = QTableWidget(outputDialog)
        if (self.probeTable.columnCount() < 4):
            self.probeTable.setColumnCount(4)
        font = QFont()
        font.setBold(True)
        __qtablewidgetitem = QTableWidgetItem()
        __qtablewidgetitem.setFont(font);
        self.probeTable.setHorizontalHeaderItem(0, __qtablewidgetitem)
        __qtablewidgetitem1 = QTableWidgetItem()
        __qtablewidgetitem1.setFont(font);
        self.probeTable.setHorizontalHeaderItem(1, __qtablewidgetitem1)
        __qtablewidgetitem2 = QTableWidgetItem()
        __qtablewidgetitem2.setFont(font);
        self.probeTable.setHorizontalHeaderItem(2, __qtablewidgetitem2)
        __qtablewidgetitem3 = QTableWidgetItem()
        __qtablewidgetitem3.setFont(font);
        self.probeTable.setHorizontalHeaderItem(3, __qtablewidgetitem3)
        self.probeTable.setObjectName(u"probeTable")
        self.probeTable.setGeometry(QRect(70, 270, 661, 261))
        self.verticalLayoutWidget = QWidget(outputDialog)
        self.verticalLayoutWidget.setObjectName(u"verticalLayoutWidget")
        self.verticalLayoutWidget.setGeometry(QRect(300, 110, 211, 151))
        self.verticalLayout_3 = QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout_3.setObjectName(u"verticalLayout_3")
        self.verticalLayout_3.setContentsMargins(0, 0, 0, 0)
        self.label = QLabel(self.verticalLayoutWidget)
        self.label.setObjectName(u"label")
        self.label.setFont(font)
        self.label.setAlignment(Qt.AlignCenter)

        self.verticalLayout_3.addWidget(self.label)

        self.targetLst = QListWidget(self.verticalLayoutWidget)
        self.targetLst.setObjectName(u"targetLst")
        self.targetLst.viewport().setProperty("cursor", QCursor(Qt.PointingHandCursor))
        self.targetLst.setEditTriggers(QAbstractItemView.EditKeyPressed|QAbstractItemView.SelectedClicked)
        self.targetLst.setProperty("showDropIndicator", False)

        self.verticalLayout_3.addWidget(self.targetLst)

        self.csvExportBtn = QPushButton(outputDialog)
        self.csvExportBtn.setObjectName(u"csvExportBtn")
        self.csvExportBtn.setGeometry(QRect(600, 550, 131, 32))
        self.csvExportBtn.setFont(font)
        self.csvExportBtn.setCursor(QCursor(Qt.PointingHandCursor))
        self.homeButton = QPushButton(outputDialog)
        self.homeButton.setObjectName(u"homeButton")
        self.homeButton.setGeometry(QRect(60, 550, 113, 32))
        self.plsUploadLbl = QLabel(outputDialog)
        self.plsUploadLbl.setObjectName(u"plsUploadLbl")
        self.plsUploadLbl.setGeometry(QRect(580, 580, 171, 20))
        self.plsUploadLbl.setAlignment(Qt.AlignCenter)

        self.retranslateUi(outputDialog)

        QMetaObject.connectSlotsByName(outputDialog)
    # setupUi

    def retranslateUi(self, outputDialog):
        outputDialog.setWindowTitle(QCoreApplication.translate("outputDialog", u"Dialog", None))
        self.logoLabel.setText("")
        ___qtablewidgetitem = self.probeTable.horizontalHeaderItem(0)
        ___qtablewidgetitem.setText(QCoreApplication.translate("outputDialog", u"SNP", None));
        ___qtablewidgetitem1 = self.probeTable.horizontalHeaderItem(1)
        ___qtablewidgetitem1.setText(QCoreApplication.translate("outputDialog", u"Variable Probe (5' to 3')", None));
        ___qtablewidgetitem2 = self.probeTable.horizontalHeaderItem(2)
        ___qtablewidgetitem2.setText(QCoreApplication.translate("outputDialog", u"Variable WT Probe (5' to 3')", None));
        ___qtablewidgetitem3 = self.probeTable.horizontalHeaderItem(3)
        ___qtablewidgetitem3.setText(QCoreApplication.translate("outputDialog", u"Common Probe (5' to 3')", None));
        self.label.setText(QCoreApplication.translate("outputDialog", u"Select Target of Interest", None))
        self.csvExportBtn.setText(QCoreApplication.translate("outputDialog", u"Export to .csv", None))
        self.homeButton.setText(QCoreApplication.translate("outputDialog", u"Return Home", None))
        self.plsUploadLbl.setText("")
    # retranslateUi

