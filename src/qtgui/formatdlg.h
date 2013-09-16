#ifndef FORMATDLG_H
#define FORMATDLG_H

#include <QDialog>

#include "ui_formatdlg.h"

class Session;

class formatdlg : public QDialog
{
    Q_OBJECT

	public:
		
		formatdlg(QWidget* parent, Session* session);
	
	protected:

		virtual void accept();

	private:

		//! The user interface for this dialog.
		Ui::formatdlg m_ui;

		//! The session we are loading data for.
		Session* m_session;

};

#endif // FORMATDLG_H
