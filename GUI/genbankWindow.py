# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file '2024 Jan - OLAgen GUI GenBank v3.ui'
##
## Created by: Qt User Interface Compiler version 6.6.1
##
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from PySide6.QtCore import (QCoreApplication, QDate, QDateTime, QLocale,
    QMetaObject, QObject, QPoint, QRect,
    QSize, QTime, QUrl, Qt)
from PySide6.QtGui import (QAction, QBrush, QColor, QConicalGradient,
    QCursor, QFont, QFontDatabase, QGradient,
    QIcon, QImage, QKeySequence, QLinearGradient,
    QPainter, QPalette, QPixmap, QRadialGradient,
    QTransform)
from PySide6.QtWidgets import (QAbstractButton, QApplication, QCheckBox, QDialogButtonBox,
    QLabel, QLineEdit, QMainWindow, QMenu,
    QMenuBar, QSizePolicy, QStatusBar, QWidget)

class Ui_gbMainWindow(object):
    def setupUi(self, gbMainWindow):
        if not gbMainWindow.objectName():
            gbMainWindow.setObjectName(u"gbMainWindow")
        gbMainWindow.resize(800, 600)
        self.actionClose = QAction(gbMainWindow)
        self.actionClose.setObjectName(u"actionClose")
        self.centralwidget = QWidget(gbMainWindow)
        self.centralwidget.setObjectName(u"centralwidget")
        self.buttonBox = QDialogButtonBox(self.centralwidget)
        self.buttonBox.setObjectName(u"buttonBox")
        self.buttonBox.setGeometry(QRect(310, 290, 164, 32))
        self.buttonBox.setStandardButtons(QDialogButtonBox.Cancel|QDialogButtonBox.Ok)
        self.lineEdit_ref_genbank = QLineEdit(self.centralwidget)
        self.lineEdit_ref_genbank.setObjectName(u"lineEdit_ref_genbank")
        self.lineEdit_ref_genbank.setGeometry(QRect(470, 140, 113, 21))
        self.lineEdit_target_genbank = QLineEdit(self.centralwidget)
        self.lineEdit_target_genbank.setObjectName(u"lineEdit_target_genbank")
        self.lineEdit_target_genbank.setGeometry(QRect(460, 190, 113, 21))
        self.label = QLabel(self.centralwidget)
        self.label.setObjectName(u"label")
        self.label.setGeometry(QRect(210, 140, 261, 20))
        self.label_2 = QLabel(self.centralwidget)
        self.label_2.setObjectName(u"label_2")
        self.label_2.setGeometry(QRect(220, 190, 241, 20))
        self.checkBox = QCheckBox(self.centralwidget)
        self.checkBox.setObjectName(u"checkBox")
        self.checkBox.setGeometry(QRect(400, 230, 87, 20))
        self.checkBox_2 = QCheckBox(self.centralwidget)
        self.checkBox_2.setObjectName(u"checkBox_2")
        self.checkBox_2.setGeometry(QRect(400, 250, 121, 20))
        self.label_3 = QLabel(self.centralwidget)
        self.label_3.setObjectName(u"label_3")
        self.label_3.setGeometry(QRect(290, 240, 101, 20))
        gbMainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QMenuBar(gbMainWindow)
        self.menubar.setObjectName(u"menubar")
        self.menubar.setGeometry(QRect(0, 0, 800, 36))
        self.menuFile = QMenu(self.menubar)
        self.menuFile.setObjectName(u"menuFile")
        gbMainWindow.setMenuBar(self.menubar)
        self.statusbar = QStatusBar(gbMainWindow)
        self.statusbar.setObjectName(u"statusbar")
        gbMainWindow.setStatusBar(self.statusbar)

        self.menubar.addAction(self.menuFile.menuAction())
        self.menuFile.addAction(self.actionClose)

        self.retranslateUi(gbMainWindow)

        QMetaObject.connectSlotsByName(gbMainWindow)
    # setupUi

    def retranslateUi(self, gbMainWindow):
        gbMainWindow.setWindowTitle(QCoreApplication.translate("gbMainWindow", u"MainWindow", None))
        self.actionClose.setText(QCoreApplication.translate("gbMainWindow", u"Close", None))
        self.label.setText(QCoreApplication.translate("gbMainWindow", u"Reference Sequence GenBank Accession:", None))
        self.label_2.setText(QCoreApplication.translate("gbMainWindow", u"Target Sequence GenBank Accession:", None))
        self.checkBox.setText(QCoreApplication.translate("gbMainWindow", u"mafft", None))
        self.checkBox_2.setText(QCoreApplication.translate("gbMainWindow", u"Clustal Omega", None))
        self.label_3.setText(QCoreApplication.translate("gbMainWindow", u"Alignment Tool:", None))
        self.menuFile.setTitle(QCoreApplication.translate("gbMainWindow", u"File", None))
    # retranslateUi

