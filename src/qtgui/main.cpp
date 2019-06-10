#include <QtWidgets/QApplication>
#include <QMessageBox>

#include "mainwindow.h"

#include <ostream>
#include <sstream>

#include <streambuf>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);

	QIcon appicon(":/tarquinicon.png");
	a.setWindowIcon(appicon);

	// the main application window
	MainWindow w;
	w.showNormal();
	//w.showMaximized();

	// and pass control the message processing loop
	return a.exec();
}

