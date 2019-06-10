#ifndef TARQUIN_INPUT_DLG_INCLUDED
#define TARQUIN_INPUT_DLG_INCLUDED

#include <QDialog>

#include "ui_inputdlg.h"

class Session;


class InputDlg : public QDialog
{
	Q_OBJECT

	public:

		InputDlg(QWidget* parent, Session* session);

		bool startFit() const;

	protected:

		//virtual void accept();

		enum fid_type_e
		{
			WATER_SUPPRESSED_FID,
			WATER_UNSUPPRESSED_FID
		};


		void LoadFID(QString filename, fid_type_e fid_type);

    private:

        bool CheckDlg();
        
        void UpdateDlg();

	private slots:

		void OnBtnOpenWS();
		
	    void OnBtnOpenWU();

        void OnBtnOpenXML();

		void OnBtnOpenCSV();

        void OnBtnSaveLCM();
		
		void OnBtnOpenParaFile();

        void OnBtnLoad();
        
        void OnBtnFit();

	private:

		QString GetPathAndSave(QString path, QString settings_name);

	private:

		//! The session we are loading data for.
		Session* m_session;

		//! The user interface for this dialog.
		Ui::InputDlg m_ui;

		//! Has the water suppressed data been loaded.
		bool m_loaded_ws;

		//! Has a precompiled basis been loaded.
		bool m_loaded_basis;

		//! Should fitting be started once the dialog is finished/
		bool m_start_fit;
};

#endif // TARQUIN_INPUT_DLG_INCLUDED
