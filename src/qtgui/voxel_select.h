#ifndef VOXELSELECT_H
#define VOXELSELECT_H

#include <QDialog>

#include "ui_voxel_select.h"

class Session;

class voxel_select : public QDialog
{
    Q_OBJECT

	public:
		
		voxel_select(QWidget* parent, Session* session);
	
	protected:

		virtual void accept();

	private:

		//! The user interface for this dialog.
		Ui::voxel_select m_ui;

		//! The session we are loading data for.
		Session* m_session;

};

#endif // VOXELSELECT_H
