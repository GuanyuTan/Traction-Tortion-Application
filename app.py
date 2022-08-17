import sys
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQt5.QtWebEngineWidgets import *

class MainWindow(QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        self.browser = QWebEngineView()
        self.browser.setUrl(QUrl('http://192.168.100.16:8501'))
        self.setCentralWidget(self.browser)
        self.showMaximized()
app = QApplication(sys.argv)
QApplication.setApplicationName('Application')
window = MainWindow()
app.exec_()