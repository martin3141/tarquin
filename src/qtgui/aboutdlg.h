#ifndef TARQUIN_ABOUT_DIALOG_INCLUDED
#define TARQUIN_ABOUT_DIALOG_INCLUDED

#include <QDialog>

#include "ui_aboutdlg.h"

class AboutDlg : public QDialog
{
	Q_OBJECT 

	public:

		AboutDlg(QWidget* parent);

	private:

		//! The user interface for the dialog.
		Ui::AboutDlg m_ui;

};

#endif // TARQUIN_ABOUT_DIALOG_INCLUDED
