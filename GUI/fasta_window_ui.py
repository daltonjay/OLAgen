# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'fasta_window.ui'
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
from PySide6.QtWidgets import (QAbstractButton, QApplication, QCheckBox, QDialog,
    QDialogButtonBox, QLabel, QListView, QPushButton,
    QSizePolicy, QVBoxLayout, QWidget)

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        if not Dialog.objectName():
            Dialog.setObjectName(u"Dialog")
        Dialog.resize(800, 600)
        font = QFont()
        font.setItalic(False)
        Dialog.setFont(font)
        self.logoLabel = QLabel(Dialog)
        self.logoLabel.setObjectName(u"logoLabel")
        self.logoLabel.setGeometry(QRect(310, 20, 191, 81))
        self.logoLabel.setPixmap(QPixmap(u"OLAgenLogoC.png"))
        self.logoLabel.setScaledContents(True)
        self.fileUplBtn = QPushButton(Dialog)
        self.fileUplBtn.setObjectName(u"fileUplBtn")
        self.fileUplBtn.setGeometry(QRect(340, 100, 113, 32))
        self.layoutWidget = QWidget(Dialog)
        self.layoutWidget.setObjectName(u"layoutWidget")
        self.layoutWidget.setGeometry(QRect(250, 140, 301, 421))
        self.verticalLayout_3 = QVBoxLayout(self.layoutWidget)
        self.verticalLayout_3.setObjectName(u"verticalLayout_3")
        self.verticalLayout_3.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout = QVBoxLayout()
        self.verticalLayout.setObjectName(u"verticalLayout")
        self.label = QLabel(self.layoutWidget)
        self.label.setObjectName(u"label")
        font1 = QFont()
        font1.setPointSize(18)
        font1.setBold(True)
        font1.setItalic(False)
        self.label.setFont(font1)
        self.label.setAlignment(Qt.AlignCenter)

        self.verticalLayout.addWidget(self.label)

        self.userFileLabel = QLabel(self.layoutWidget)
        self.userFileLabel.setObjectName(u"userFileLabel")
        font2 = QFont()
        font2.setFamilies([u".AppleSystemUIFont"])
        font2.setPointSize(15)
        font2.setBold(False)
        font2.setItalic(False)
        self.userFileLabel.setFont(font2)
        self.userFileLabel.setStyleSheet(u"font: 15pt \".AppleSystemUIFont\";")
        self.userFileLabel.setAlignment(Qt.AlignCenter)

        self.verticalLayout.addWidget(self.userFileLabel)


        self.verticalLayout_3.addLayout(self.verticalLayout)

        self.label_2 = QLabel(self.layoutWidget)
        self.label_2.setObjectName(u"label_2")
        font3 = QFont()
        font3.setPointSize(16)
        font3.setBold(True)
        font3.setItalic(False)
        self.label_2.setFont(font3)
        self.label_2.setAlignment(Qt.AlignCenter)

        self.verticalLayout_3.addWidget(self.label_2)

        self.label_4 = QLabel(self.layoutWidget)
        self.label_4.setObjectName(u"label_4")
        self.label_4.setAlignment(Qt.AlignCenter)

        self.verticalLayout_3.addWidget(self.label_4)

        self.entryIDList = QListView(self.layoutWidget)
        self.entryIDList.setObjectName(u"entryIDList")
        font4 = QFont()
        font4.setPointSize(14)
        font4.setItalic(False)
        self.entryIDList.setFont(font4)

        self.verticalLayout_3.addWidget(self.entryIDList)

        self.verticalLayout_2 = QVBoxLayout()
        self.verticalLayout_2.setObjectName(u"verticalLayout_2")
        self.label_3 = QLabel(self.layoutWidget)
        self.label_3.setObjectName(u"label_3")
        self.label_3.setFont(font1)
        self.label_3.setAlignment(Qt.AlignCenter)

        self.verticalLayout_2.addWidget(self.label_3)

        self.mafftCheck = QCheckBox(self.layoutWidget)
        self.mafftCheck.setObjectName(u"mafftCheck")
        font5 = QFont()
        font5.setPointSize(15)
        font5.setItalic(False)
        self.mafftCheck.setFont(font5)
        self.mafftCheck.setAutoFillBackground(False)
        self.mafftCheck.setTristate(False)

        self.verticalLayout_2.addWidget(self.mafftCheck, 0, Qt.AlignHCenter)

        self.clustalCheck = QCheckBox(self.layoutWidget)
        self.clustalCheck.setObjectName(u"clustalCheck")
        self.clustalCheck.setFont(font5)
        self.clustalCheck.setTristate(False)

        self.verticalLayout_2.addWidget(self.clustalCheck, 0, Qt.AlignHCenter)

        self.alignBtnBox = QDialogButtonBox(self.layoutWidget)
        self.alignBtnBox.setObjectName(u"alignBtnBox")
        self.alignBtnBox.setOrientation(Qt.Vertical)
        self.alignBtnBox.setStandardButtons(QDialogButtonBox.Cancel|QDialogButtonBox.Ok)

        self.verticalLayout_2.addWidget(self.alignBtnBox, 0, Qt.AlignHCenter)


        self.verticalLayout_3.addLayout(self.verticalLayout_2)

        self.statusLabel = QLabel(self.layoutWidget)
        self.statusLabel.setObjectName(u"statusLabel")
        self.statusLabel.setFont(font2)
        self.statusLabel.setStyleSheet(u"font: 15pt \".AppleSystemUIFont\";")
        self.statusLabel.setAlignment(Qt.AlignCenter)

        self.verticalLayout_3.addWidget(self.statusLabel)


        self.retranslateUi(Dialog)

        QMetaObject.connectSlotsByName(Dialog)
    # setupUi

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(QCoreApplication.translate("Dialog", u"Dialog", None))
        self.logoLabel.setText("")
        self.fileUplBtn.setText(QCoreApplication.translate("Dialog", u"Upload .fasta", None))
        self.label.setText(QCoreApplication.translate("Dialog", u".fasta file:", None))
        self.userFileLabel.setText(QCoreApplication.translate("Dialog", u"...", None))
        self.label_2.setText(QCoreApplication.translate("Dialog", u"User Input Entry IDs", None))
        self.label_4.setText(QCoreApplication.translate("Dialog", u"Please Select the Reference Sequence", None))
        self.label_3.setText(QCoreApplication.translate("Dialog", u"Alignment Options:", None))
        self.mafftCheck.setText(QCoreApplication.translate("Dialog", u"mafft Alignment", None))
        self.clustalCheck.setText(QCoreApplication.translate("Dialog", u"Clustal Omega Alignment", None))
        self.statusLabel.setText("")
    # retranslateUi

