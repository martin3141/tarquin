#include "guicommon.h"
#include "CBoswell.hpp"

#include <QtGui>

/*!
 * \brief Obtain the application log for the GUI.
 */
tarquin::CBoswell& GetLog()
{
	static tarquin::CBoswell log(tarquin::LOG_STDOUT);
	return log;
}

/*!
 * \brief The dialog to show when an error occurs, usually as the result of an exception.
 */
void ErrorDialog(QWidget* parent, const QString& title, const QString& message)
{
	QMessageBox dlg(parent);

	dlg.setText(message);
	dlg.setWindowTitle(title);

	// TODO: have logo on display here as well

	dlg.exec();
}

/*!
 * \brief A dialog to display when the user needs to see something useful.
 */
void InfoDialog(QWidget* parent, const QString& title, const QString& message)
{
	QMessageBox dlg(parent);

	dlg.setText(message);
	dlg.setWindowTitle(title);

	// TODO: have logo on display here as well

	dlg.exec();
}
