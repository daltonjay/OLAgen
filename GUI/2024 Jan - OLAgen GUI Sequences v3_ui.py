# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file '2024 Jan - OLAgen GUI Sequences v3.ui'
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
from PySide6.QtWidgets import (QApplication, QMainWindow, QMenu, QMenuBar,
    QSizePolicy, QStatusBar, QWidget)

class Ui_genbankWindow(object):
    def setupUi(self, genbankWindow):
        if not genbankWindow.objectName():
            genbankWindow.setObjectName(u"genbankWindow")
        genbankWindow.resize(800, 600)
        self.actionClose = QAction(genbankWindow)
        self.actionClose.setObjectName(u"actionClose")
        self.centralwidget = QWidget(genbankWindow)
        self.centralwidget.setObjectName(u"centralwidget")
        genbankWindow.setCentralWidget(self.centralwidget)
        self.menubar = QMenuBar(genbankWindow)
        self.menubar.setObjectName(u"menubar")
        self.menubar.setGeometry(QRect(0, 0, 800, 36))
        self.menuFile = QMenu(self.menubar)
        self.menuFile.setObjectName(u"menuFile")
        genbankWindow.setMenuBar(self.menubar)
        self.statusbar = QStatusBar(genbankWindow)
        self.statusbar.setObjectName(u"statusbar")
        genbankWindow.setStatusBar(self.statusbar)

        self.menubar.addAction(self.menuFile.menuAction())
        self.menuFile.addAction(self.actionClose)

        self.retranslateUi(genbankWindow)

        QMetaObject.connectSlotsByName(genbankWindow)
    # setupUi

    def retranslateUi(self, genbankWindow):
        genbankWindow.setWindowTitle(QCoreApplication.translate("genbankWindow", u"MainWindow", None))
        self.actionClose.setText(QCoreApplication.translate("genbankWindow", u"Close", None))
        self.menuFile.setTitle(QCoreApplication.translate("genbankWindow", u"File", None))
    # retranslateUi

