#ifndef RESULTS_DLG_INCLUDED
#define RESULTS_DLG_INCLUDED

#include "ui_resultsdlg.h"

class QStandardItemModel;

namespace tarquin
{
	class Workspace;

} // namespace tarquin

class Session;

class ResultsDlg : public QDialog
{
	Q_OBJECT

	public:

		ResultsDlg(QWidget* parent, Session* session);

		void InitFromWorkspace(const tarquin::Workspace& workspace);

	private:

		//! The user interface for this dialog.
		Ui::ResultsDlg m_ui;

        //! The session we are loading data for.
		Session* m_session;

		//! The model we are using to display the results.
		QStandardItemModel* m_model;
        
};


#endif // RESULTS_DLG_INCLUDED
