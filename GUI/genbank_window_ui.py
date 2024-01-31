# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'genbank_window.ui'
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
    QDialogButtonBox, QLabel, QLineEdit, QSizePolicy,
    QWidget)

class Ui_genbankWindow(object):
    def setupUi(self, genbankWindow):
        if not genbankWindow.objectName():
            genbankWindow.setObjectName(u"genbankWindow")
        genbankWindow.resize(800, 600)
        self.lineEdit_ref_genbank = QLineEdit(genbankWindow)
        self.lineEdit_ref_genbank.setObjectName(u"lineEdit_ref_genbank")
        self.lineEdit_ref_genbank.setGeometry(QRect(480, 170, 113, 21))
        self.genbuttonBox = QDialogButtonBox(genbankWindow)
        self.genbuttonBox.setObjectName(u"genbuttonBox")
        self.genbuttonBox.setGeometry(QRect(320, 320, 164, 32))
        self.genbuttonBox.setStandardButtons(QDialogButtonBox.Cancel|QDialogButtonBox.Ok)
        self.label = QLabel(genbankWindow)
        self.label.setObjectName(u"label")
        self.label.setGeometry(QRect(220, 170, 261, 20))
        self.checkBox = QCheckBox(genbankWindow)
        self.checkBox.setObjectName(u"checkBox")
        self.checkBox.setGeometry(QRect(410, 260, 87, 20))
        self.checkBox_2 = QCheckBox(genbankWindow)
        self.checkBox_2.setObjectName(u"checkBox_2")
        self.checkBox_2.setGeometry(QRect(410, 280, 121, 20))
        self.label_2 = QLabel(genbankWindow)
        self.label_2.setObjectName(u"label_2")
        self.label_2.setGeometry(QRect(230, 220, 241, 20))
        self.label_3 = QLabel(genbankWindow)
        self.label_3.setObjectName(u"label_3")
        self.label_3.setGeometry(QRect(300, 270, 101, 20))
        self.lineEdit_target_genbank = QLineEdit(genbankWindow)
        self.lineEdit_target_genbank.setObjectName(u"lineEdit_target_genbank")
        self.lineEdit_target_genbank.setGeometry(QRect(470, 220, 113, 21))

        self.retranslateUi(genbankWindow)

        QMetaObject.connectSlotsByName(genbankWindow)
    # setupUi

    def retranslateUi(self, genbankWindow):
        genbankWindow.setWindowTitle(QCoreApplication.translate("genbankWindow", u"Dialog", None))
        self.label.setText(QCoreApplication.translate("genbankWindow", u"Reference Sequence GenBank Accession:", None))
        self.checkBox.setText(QCoreApplication.translate("genbankWindow", u"mafft", None))
        self.checkBox_2.setText(QCoreApplication.translate("genbankWindow", u"Clustal Omega", None))
        self.label_2.setText(QCoreApplication.translate("genbankWindow", u"Target Sequence GenBank Accession:", None))
        self.label_3.setText(QCoreApplication.translate("genbankWindow", u"Alignment Tool:", None))
    # retranslateUi

