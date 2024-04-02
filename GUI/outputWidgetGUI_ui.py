# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'outputWidgetGUI.ui'
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
from PySide6.QtWidgets import (QAbstractItemView, QApplication, QHeaderView, QLabel,
    QListWidget, QListWidgetItem, QPushButton, QSizePolicy,
    QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget)

class Ui_outputWidget(object):
    def setupUi(self, outputWidget):
        if not outputWidget.objectName():
            outputWidget.setObjectName(u"outputWidget")
        outputWidget.resize(800, 600)
        self.verticalLayoutWidget = QWidget(outputWidget)
        self.verticalLayoutWidget.setObjectName(u"verticalLayoutWidget")
        self.verticalLayoutWidget.setGeometry(QRect(80, 130, 211, 261))
        self.verticalLayout_3 = QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout_3.setObjectName(u"verticalLayout_3")
        self.verticalLayout_3.setContentsMargins(0, 0, 0, 0)
        self.label = QLabel(self.verticalLayoutWidget)
        self.label.setObjectName(u"label")
        font = QFont()
        font.setBold(True)
        self.label.setFont(font)
        self.label.setAlignment(Qt.AlignCenter)

        self.verticalLayout_3.addWidget(self.label)

        self.targetLst = QListWidget(self.verticalLayoutWidget)
        self.targetLst.setObjectName(u"targetLst")
        self.targetLst.viewport().setProperty("cursor", QCursor(Qt.PointingHandCursor))
        self.targetLst.setEditTriggers(QAbstractItemView.EditKeyPressed|QAbstractItemView.SelectedClicked)
        self.targetLst.setProperty("showDropIndicator", False)

        self.verticalLayout_3.addWidget(self.targetLst)

        self.probeTable = QTableWidget(outputWidget)
        if (self.probeTable.columnCount() < 3):
            self.probeTable.setColumnCount(3)
        __qtablewidgetitem = QTableWidgetItem()
        __qtablewidgetitem.setFont(font);
        self.probeTable.setHorizontalHeaderItem(0, __qtablewidgetitem)
        __qtablewidgetitem1 = QTableWidgetItem()
        __qtablewidgetitem1.setFont(font);
        self.probeTable.setHorizontalHeaderItem(1, __qtablewidgetitem1)
        __qtablewidgetitem2 = QTableWidgetItem()
        __qtablewidgetitem2.setFont(font);
        self.probeTable.setHorizontalHeaderItem(2, __qtablewidgetitem2)
        self.probeTable.setObjectName(u"probeTable")
        self.probeTable.setGeometry(QRect(320, 130, 401, 261))
        self.csvExportBtn = QPushButton(outputWidget)
        self.csvExportBtn.setObjectName(u"csvExportBtn")
        self.csvExportBtn.setGeometry(QRect(582, 400, 131, 32))
        self.csvExportBtn.setFont(font)
        self.csvExportBtn.setCursor(QCursor(Qt.PointingHandCursor))
        self.logoLabel = QLabel(outputWidget)
        self.logoLabel.setObjectName(u"logoLabel")
        self.logoLabel.setGeometry(QRect(340, 20, 151, 61))
        self.logoLabel.setPixmap(QPixmap(u"OLAgenLogo.png"))
        self.logoLabel.setScaledContents(True)

        self.retranslateUi(outputWidget)

        QMetaObject.connectSlotsByName(outputWidget)
    # setupUi

    def retranslateUi(self, outputWidget):
        outputWidget.setWindowTitle(QCoreApplication.translate("outputWidget", u"Form", None))
        self.label.setText(QCoreApplication.translate("outputWidget", u"Select Target of Interest", None))
        ___qtablewidgetitem = self.probeTable.horizontalHeaderItem(0)
        ___qtablewidgetitem.setText(QCoreApplication.translate("outputWidget", u"SNP", None));
        ___qtablewidgetitem1 = self.probeTable.horizontalHeaderItem(1)
        ___qtablewidgetitem1.setText(QCoreApplication.translate("outputWidget", u"Common Probe", None));
        ___qtablewidgetitem2 = self.probeTable.horizontalHeaderItem(2)
        ___qtablewidgetitem2.setText(QCoreApplication.translate("outputWidget", u"Variable Probe", None));
        self.csvExportBtn.setText(QCoreApplication.translate("outputWidget", u"Export to .csv", None))
        self.logoLabel.setText("")
    # retranslateUi

