# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file '2024 Jan - OLAgen GUI Alignment Progress v3.ui'
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
    QMenuBar, QProgressBar, QSizePolicy, QStatusBar,
    QWidget)

class Ui_progressMainWindow(object):
    def setupUi(self, progressMainWindow):
        if not progressMainWindow.objectName():
            progressMainWindow.setObjectName(u"progressMainWindow")
        progressMainWindow.resize(800, 600)
        self.actionClose = QAction(progressMainWindow)
        self.actionClose.setObjectName(u"actionClose")
        self.centralwidget = QWidget(progressMainWindow)
        self.centralwidget.setObjectName(u"centralwidget")
        self.progressBar = QProgressBar(self.centralwidget)
        self.progressBar.setObjectName(u"progressBar")
        self.progressBar.setGeometry(QRect(140, 170, 521, 151))
        self.progressBar.setValue(24)
        self.label = QLabel(self.centralwidget)
        self.label.setObjectName(u"label")
        self.label.setGeometry(QRect(340, 190, 121, 20))
        progressMainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QMenuBar(progressMainWindow)
        self.menubar.setObjectName(u"menubar")
        self.menubar.setGeometry(QRect(0, 0, 800, 36))
        self.menuFile = QMenu(self.menubar)
        self.menuFile.setObjectName(u"menuFile")
        progressMainWindow.setMenuBar(self.menubar)
        self.statusbar = QStatusBar(progressMainWindow)
        self.statusbar.setObjectName(u"statusbar")
        progressMainWindow.setStatusBar(self.statusbar)

        self.menubar.addAction(self.menuFile.menuAction())
        self.menuFile.addAction(self.actionClose)

        self.retranslateUi(progressMainWindow)

        QMetaObject.connectSlotsByName(progressMainWindow)
    # setupUi

    def retranslateUi(self, progressMainWindow):
        progressMainWindow.setWindowTitle(QCoreApplication.translate("progressMainWindow", u"MainWindow", None))
        self.actionClose.setText(QCoreApplication.translate("progressMainWindow", u"Close", None))
        self.label.setText(QCoreApplication.translate("progressMainWindow", u"Alignment Progress", None))
        self.menuFile.setTitle(QCoreApplication.translate("progressMainWindow", u"File", None))
    # retranslateUi

