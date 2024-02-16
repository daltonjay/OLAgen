# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'output_windowTAB.ui'
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
    QSizePolicy, QTabWidget, QTableWidget, QTableWidgetItem,
    QVBoxLayout, QWidget)

class Ui_outputDialog(object):
    def setupUi(self, outputDialog):
        if not outputDialog.objectName():
            outputDialog.setObjectName(u"outputDialog")
        outputDialog.resize(800, 600)
        self.logoLabel = QLabel(outputDialog)
        self.logoLabel.setObjectName(u"logoLabel")
        self.logoLabel.setGeometry(QRect(330, 10, 151, 61))
        self.logoLabel.setPixmap(QPixmap(u"OLAgenLogo.png"))
        self.logoLabel.setScaledContents(True)
        self.verticalLayoutWidget = QWidget(outputDialog)
        self.verticalLayoutWidget.setObjectName(u"verticalLayoutWidget")
        self.verticalLayoutWidget.setGeometry(QRect(290, 90, 211, 121))
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
        self.tabWidget = QTabWidget(outputDialog)
        self.tabWidget.setObjectName(u"tabWidget")
        self.tabWidget.setGeometry(QRect(40, 220, 721, 301))
        self.ligTab = QWidget()
        self.ligTab.setObjectName(u"ligTab")
        self.probeTable = QTableWidget(self.ligTab)
        if (self.probeTable.columnCount() < 4):
            self.probeTable.setColumnCount(4)
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
        self.probeTable.setGeometry(QRect(20, 0, 671, 231))
        self.addSOIBtn = QPushButton(self.ligTab)
        self.addSOIBtn.setObjectName(u"addSOIBtn")
        self.addSOIBtn.setGeometry(QRect(502, 240, 151, 32))
        self.tabWidget.addTab(self.ligTab, "")
        self.soiTab = QWidget()
        self.soiTab.setObjectName(u"soiTab")
        self.soiTable = QTableWidget(self.soiTab)
        if (self.soiTable.columnCount() < 4):
            self.soiTable.setColumnCount(4)
        __qtablewidgetitem4 = QTableWidgetItem()
        __qtablewidgetitem4.setFont(font);
        self.soiTable.setHorizontalHeaderItem(0, __qtablewidgetitem4)
        __qtablewidgetitem5 = QTableWidgetItem()
        __qtablewidgetitem5.setFont(font);
        self.soiTable.setHorizontalHeaderItem(1, __qtablewidgetitem5)
        __qtablewidgetitem6 = QTableWidgetItem()
        __qtablewidgetitem6.setFont(font);
        self.soiTable.setHorizontalHeaderItem(2, __qtablewidgetitem6)
        __qtablewidgetitem7 = QTableWidgetItem()
        __qtablewidgetitem7.setFont(font);
        self.soiTable.setHorizontalHeaderItem(3, __qtablewidgetitem7)
        self.soiTable.setObjectName(u"soiTable")
        self.soiTable.setGeometry(QRect(20, 0, 671, 231))
        self.genPrimerBtn = QPushButton(self.soiTab)
        self.genPrimerBtn.setObjectName(u"genPrimerBtn")
        self.genPrimerBtn.setGeometry(QRect(500, 240, 151, 32))
        self.removeSetBtn = QPushButton(self.soiTab)
        self.removeSetBtn.setObjectName(u"removeSetBtn")
        self.removeSetBtn.setGeometry(QRect(60, 240, 151, 32))
        self.tabWidget.addTab(self.soiTab, "")
        self.primerTab = QWidget()
        self.primerTab.setObjectName(u"primerTab")
        self.primerTable = QTableWidget(self.primerTab)
        if (self.primerTable.columnCount() < 6):
            self.primerTable.setColumnCount(6)
        __qtablewidgetitem8 = QTableWidgetItem()
        __qtablewidgetitem8.setFont(font);
        self.primerTable.setHorizontalHeaderItem(0, __qtablewidgetitem8)
        __qtablewidgetitem9 = QTableWidgetItem()
        __qtablewidgetitem9.setFont(font);
        self.primerTable.setHorizontalHeaderItem(1, __qtablewidgetitem9)
        __qtablewidgetitem10 = QTableWidgetItem()
        __qtablewidgetitem10.setFont(font);
        self.primerTable.setHorizontalHeaderItem(2, __qtablewidgetitem10)
        __qtablewidgetitem11 = QTableWidgetItem()
        __qtablewidgetitem11.setFont(font);
        self.primerTable.setHorizontalHeaderItem(3, __qtablewidgetitem11)
        __qtablewidgetitem12 = QTableWidgetItem()
        __qtablewidgetitem12.setFont(font);
        self.primerTable.setHorizontalHeaderItem(4, __qtablewidgetitem12)
        __qtablewidgetitem13 = QTableWidgetItem()
        __qtablewidgetitem13.setFont(font);
        self.primerTable.setHorizontalHeaderItem(5, __qtablewidgetitem13)
        self.primerTable.setObjectName(u"primerTable")
        self.primerTable.setGeometry(QRect(20, 0, 671, 231))
        self.removePrimerBtn = QPushButton(self.primerTab)
        self.removePrimerBtn.setObjectName(u"removePrimerBtn")
        self.removePrimerBtn.setGeometry(QRect(60, 240, 151, 32))
        self.genFullReagent = QPushButton(self.primerTab)
        self.genFullReagent.setObjectName(u"genFullReagent")
        self.genFullReagent.setGeometry(QRect(490, 240, 171, 32))
        self.tabWidget.addTab(self.primerTab, "")
        self.allTab = QWidget()
        self.allTab.setObjectName(u"allTab")
        self.reagentTable = QTableWidget(self.allTab)
        if (self.reagentTable.columnCount() < 11):
            self.reagentTable.setColumnCount(11)
        __qtablewidgetitem14 = QTableWidgetItem()
        __qtablewidgetitem14.setFont(font);
        self.reagentTable.setHorizontalHeaderItem(0, __qtablewidgetitem14)
        __qtablewidgetitem15 = QTableWidgetItem()
        __qtablewidgetitem15.setFont(font);
        self.reagentTable.setHorizontalHeaderItem(1, __qtablewidgetitem15)
        __qtablewidgetitem16 = QTableWidgetItem()
        __qtablewidgetitem16.setFont(font);
        self.reagentTable.setHorizontalHeaderItem(2, __qtablewidgetitem16)
        __qtablewidgetitem17 = QTableWidgetItem()
        __qtablewidgetitem17.setFont(font);
        self.reagentTable.setHorizontalHeaderItem(3, __qtablewidgetitem17)
        __qtablewidgetitem18 = QTableWidgetItem()
        __qtablewidgetitem18.setFont(font);
        self.reagentTable.setHorizontalHeaderItem(4, __qtablewidgetitem18)
        __qtablewidgetitem19 = QTableWidgetItem()
        __qtablewidgetitem19.setFont(font);
        self.reagentTable.setHorizontalHeaderItem(5, __qtablewidgetitem19)
        __qtablewidgetitem20 = QTableWidgetItem()
        __qtablewidgetitem20.setFont(font);
        self.reagentTable.setHorizontalHeaderItem(6, __qtablewidgetitem20)
        __qtablewidgetitem21 = QTableWidgetItem()
        __qtablewidgetitem21.setFont(font);
        self.reagentTable.setHorizontalHeaderItem(7, __qtablewidgetitem21)
        __qtablewidgetitem22 = QTableWidgetItem()
        __qtablewidgetitem22.setFont(font);
        self.reagentTable.setHorizontalHeaderItem(8, __qtablewidgetitem22)
        __qtablewidgetitem23 = QTableWidgetItem()
        __qtablewidgetitem23.setFont(font);
        self.reagentTable.setHorizontalHeaderItem(9, __qtablewidgetitem23)
        __qtablewidgetitem24 = QTableWidgetItem()
        __qtablewidgetitem24.setFont(font);
        self.reagentTable.setHorizontalHeaderItem(10, __qtablewidgetitem24)
        self.reagentTable.setObjectName(u"reagentTable")
        self.reagentTable.setGeometry(QRect(20, 10, 671, 231))
        self.tabWidget.addTab(self.allTab, "")

        self.retranslateUi(outputDialog)

        self.tabWidget.setCurrentIndex(2)


        QMetaObject.connectSlotsByName(outputDialog)
    # setupUi

    def retranslateUi(self, outputDialog):
        outputDialog.setWindowTitle(QCoreApplication.translate("outputDialog", u"Dialog", None))
        self.logoLabel.setText("")
        self.label.setText(QCoreApplication.translate("outputDialog", u"Select Target of Interest", None))
        self.csvExportBtn.setText(QCoreApplication.translate("outputDialog", u"Export to .csv", None))
        self.homeButton.setText(QCoreApplication.translate("outputDialog", u"Return Home", None))
        self.plsUploadLbl.setText("")
        ___qtablewidgetitem = self.probeTable.horizontalHeaderItem(0)
        ___qtablewidgetitem.setText(QCoreApplication.translate("outputDialog", u"SNP", None));
        ___qtablewidgetitem1 = self.probeTable.horizontalHeaderItem(1)
        ___qtablewidgetitem1.setText(QCoreApplication.translate("outputDialog", u"Variable Probe (5' to 3')", None));
        ___qtablewidgetitem2 = self.probeTable.horizontalHeaderItem(2)
        ___qtablewidgetitem2.setText(QCoreApplication.translate("outputDialog", u"Variable WT Probe (5' to 3')", None));
        ___qtablewidgetitem3 = self.probeTable.horizontalHeaderItem(3)
        ___qtablewidgetitem3.setText(QCoreApplication.translate("outputDialog", u"Common Probe (5' to 3')", None));
        self.addSOIBtn.setText(QCoreApplication.translate("outputDialog", u"Add Set of Interest", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.ligTab), QCoreApplication.translate("outputDialog", u"Ligation Probes", None))
        ___qtablewidgetitem4 = self.soiTable.horizontalHeaderItem(0)
        ___qtablewidgetitem4.setText(QCoreApplication.translate("outputDialog", u"Set ID", None));
        ___qtablewidgetitem5 = self.soiTable.horizontalHeaderItem(1)
        ___qtablewidgetitem5.setText(QCoreApplication.translate("outputDialog", u"Target", None));
        ___qtablewidgetitem6 = self.soiTable.horizontalHeaderItem(2)
        ___qtablewidgetitem6.setText(QCoreApplication.translate("outputDialog", u"SNP", None));
        ___qtablewidgetitem7 = self.soiTable.horizontalHeaderItem(3)
        ___qtablewidgetitem7.setText(QCoreApplication.translate("outputDialog", u"Primers Generated?", None));
        self.genPrimerBtn.setText(QCoreApplication.translate("outputDialog", u"Generate Primers", None))
        self.removeSetBtn.setText(QCoreApplication.translate("outputDialog", u"Remove Set", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.soiTab), QCoreApplication.translate("outputDialog", u"Sets of Interest", None))
        ___qtablewidgetitem8 = self.primerTable.horizontalHeaderItem(0)
        ___qtablewidgetitem8.setText(QCoreApplication.translate("outputDialog", u"Set ID", None));
        ___qtablewidgetitem9 = self.primerTable.horizontalHeaderItem(1)
        ___qtablewidgetitem9.setText(QCoreApplication.translate("outputDialog", u"Forward Primer", None));
        ___qtablewidgetitem10 = self.primerTable.horizontalHeaderItem(2)
        ___qtablewidgetitem10.setText(QCoreApplication.translate("outputDialog", u"Reverse Primer", None));
        ___qtablewidgetitem11 = self.primerTable.horizontalHeaderItem(3)
        ___qtablewidgetitem11.setText(QCoreApplication.translate("outputDialog", u"Hydrolysis Probe", None));
        ___qtablewidgetitem12 = self.primerTable.horizontalHeaderItem(4)
        ___qtablewidgetitem12.setText(QCoreApplication.translate("outputDialog", u"Primer Type", None));
        ___qtablewidgetitem13 = self.primerTable.horizontalHeaderItem(5)
        ___qtablewidgetitem13.setText(QCoreApplication.translate("outputDialog", u"Citation", None));
        self.removePrimerBtn.setText(QCoreApplication.translate("outputDialog", u"Remove Candidate", None))
        self.genFullReagent.setText(QCoreApplication.translate("outputDialog", u"Generate Full Reagents", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.primerTab), QCoreApplication.translate("outputDialog", u"Primer Candidates", None))
        ___qtablewidgetitem14 = self.reagentTable.horizontalHeaderItem(0)
        ___qtablewidgetitem14.setText(QCoreApplication.translate("outputDialog", u"Set ID", None));
        ___qtablewidgetitem15 = self.reagentTable.horizontalHeaderItem(1)
        ___qtablewidgetitem15.setText(QCoreApplication.translate("outputDialog", u"Target", None));
        ___qtablewidgetitem16 = self.reagentTable.horizontalHeaderItem(2)
        ___qtablewidgetitem16.setText(QCoreApplication.translate("outputDialog", u"SNP", None));
        ___qtablewidgetitem17 = self.reagentTable.horizontalHeaderItem(3)
        ___qtablewidgetitem17.setText(QCoreApplication.translate("outputDialog", u"VP+Fwd+Hydr", None));
        ___qtablewidgetitem18 = self.reagentTable.horizontalHeaderItem(4)
        ___qtablewidgetitem18.setText(QCoreApplication.translate("outputDialog", u"VP WT+Fwd", None));
        ___qtablewidgetitem19 = self.reagentTable.horizontalHeaderItem(5)
        ___qtablewidgetitem19.setText(QCoreApplication.translate("outputDialog", u"CP+Rev", None));
        ___qtablewidgetitem20 = self.reagentTable.horizontalHeaderItem(6)
        ___qtablewidgetitem20.setText(QCoreApplication.translate("outputDialog", u"Forward Primer", None));
        ___qtablewidgetitem21 = self.reagentTable.horizontalHeaderItem(7)
        ___qtablewidgetitem21.setText(QCoreApplication.translate("outputDialog", u"Reverse Primer", None));
        ___qtablewidgetitem22 = self.reagentTable.horizontalHeaderItem(8)
        ___qtablewidgetitem22.setText(QCoreApplication.translate("outputDialog", u"Hydrolysis Probe", None));
        ___qtablewidgetitem23 = self.reagentTable.horizontalHeaderItem(9)
        ___qtablewidgetitem23.setText(QCoreApplication.translate("outputDialog", u"Primer Type", None));
        ___qtablewidgetitem24 = self.reagentTable.horizontalHeaderItem(10)
        ___qtablewidgetitem24.setText(QCoreApplication.translate("outputDialog", u"Citation", None));
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.allTab), QCoreApplication.translate("outputDialog", u"All Reagents", None))
    # retranslateUi

