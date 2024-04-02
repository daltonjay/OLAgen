# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file '2024 Jan - OLAgen GUI v2.ui'
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
from PySide6.QtWidgets import (QApplication, QLabel, QMainWindow, QMenu,
    QMenuBar, QPushButton, QSizePolicy, QStatusBar,
    QWidget)

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        if not MainWindow.objectName():
            MainWindow.setObjectName(u"MainWindow")
        MainWindow.resize(800, 600)
        self.actionClose = QAction(MainWindow)
        self.actionClose.setObjectName(u"actionClose")
        self.centralwidget = QWidget(MainWindow)
        self.centralwidget.setObjectName(u"centralwidget")
        self.olagenLogo = QLabel(self.centralwidget)
        self.olagenLogo.setObjectName(u"olagenLogo")
        self.olagenLogo.setGeometry(QRect(190, 40, 461, 201))
        self.olagenLogo.setPixmap(QPixmap(u"OLAgenLogo.png"))
        self.olagenLogo.setScaledContents(True)
        self.fastaInputButton = QPushButton(self.centralwidget)
        self.fastaInputButton.setObjectName(u"fastaInputButton")
        self.fastaInputButton.setGeometry(QRect(330, 230, 141, 26))
        self.helpButton = QPushButton(self.centralwidget)
        self.helpButton.setObjectName(u"helpButton")
        self.helpButton.setGeometry(QRect(360, 260, 81, 26))
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QMenuBar(MainWindow)
        self.menubar.setObjectName(u"menubar")
        self.menubar.setGeometry(QRect(0, 0, 800, 36))
        self.menuFile = QMenu(self.menubar)
        self.menuFile.setObjectName(u"menuFile")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QStatusBar(MainWindow)
        self.statusbar.setObjectName(u"statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.menubar.addAction(self.menuFile.menuAction())
        self.menuFile.addAction(self.actionClose)

        self.retranslateUi(MainWindow)

        QMetaObject.connectSlotsByName(MainWindow)
    # setupUi

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QCoreApplication.translate("MainWindow", u"MainWindow", None))
        self.actionClose.setText(QCoreApplication.translate("MainWindow", u"Close", None))
        self.olagenLogo.setText("")
        self.fastaInputButton.setText(QCoreApplication.translate("MainWindow", u"Select .fasta Files", None))
        self.helpButton.setText(QCoreApplication.translate("MainWindow", u"Help", None))
        self.menuFile.setTitle(QCoreApplication.translate("MainWindow", u"File", None))
    # retranslateUi

