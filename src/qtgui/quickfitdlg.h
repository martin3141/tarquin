#ifndef QUICKFITDLG_H
#define QUICKFITDLG_H

#include <QDialog>

#include "ui_quickfitdlg.h"

class Session;

class quickfitdlg : public QDialog
{
    Q_OBJECT

	public:
		
		quickfitdlg(QWidget* parent, Session* session);
	
	protected:

		virtual void accept();

		enum fid_type_e
		{
			WATER_SUPPRESSED_FID,
			WATER_UNSUPPRESSED_FID
		};


		void LoadFID(QString filename, fid_type_e fid_type);

	private slots:

		void OnBtnOpenWS();
		
	    void OnBtnOpenWU();

	private:

		QString GetPathAndSave(QString path, QString settings_name);
		
		//! The user interface for this dialog.
		Ui::quickfitdlg m_ui;

		//! The session we are loading data for.
		Session* m_session;

		//! Has the water suppressed data been loaded.
		bool m_loaded_ws;

};

#endif // QUICKFITDLG_H
