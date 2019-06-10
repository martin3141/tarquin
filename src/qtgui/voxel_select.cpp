#include "session.h"
#include "exception.hpp"
#include "voxel_select.h"
#include "guicommon.h"
#include "CFID.hpp"

#include <boost/filesystem.hpp>
#include <QtWidgets>


voxel_select::voxel_select(QWidget *parent, Session* session) :
    QDialog(parent),
	m_session(session)
{
	assert( session );
    m_ui.setupUi(this);

	tarquin::CFID& fid = m_session->GetWorkspace().GetFID();
	tarquin::Options& options = m_session->GetWorkspace().GetOptions();

	m_ui.txtRowsAvail->setText(QString::number(fid.GetRows()));
	m_ui.txtColsAvail->setText(QString::number(fid.GetCols()));
	m_ui.txtSlicesAvail->setText(QString::number(fid.GetSlices()));

	m_ui.txtRows->setText(QString::number(options.GetFitRows()));
	m_ui.txtCols->setText(QString::number(options.GetFitCols()));
	m_ui.txtSlices->setText(QString::number(options.GetFitSlices()));
}

void voxel_select::accept()
{
	bool ok = false;
    int rows = m_ui.txtRows->text().toInt(&ok);
    int rows_avail = m_ui.txtRowsAvail->text().toInt();
	if( !ok || rows > rows_avail || rows < 0 )
    {
        QMessageBox::information(this, tr("Error"), tr("Please correct the number of rows selected for fitting."));
        return;
    }

    int cols = m_ui.txtCols->text().toInt(&ok);
    int cols_avail = m_ui.txtColsAvail->text().toInt();
	if( !ok || cols > cols_avail || cols < 0 )
    {
        QMessageBox::information(this, tr("Error"), tr("Please correct the number of columns selected for fitting."));
        return;
    }

    int slices = m_ui.txtSlices->text().toInt(&ok);
    int slices_avail = m_ui.txtSlicesAvail->text().toInt();
	if( !ok || slices > slices_avail || slices < 0 )
    {
        QMessageBox::information(this, tr("Error"), tr("Please correct the number of slices selected for fitting."));
        return;
    }
    
    tarquin::coord_vec fit_list;
    fit_list.clear();
        
    int row_start = floor(rows_avail/2.0 - rows/2.0 + 1.0);
    int row_end = floor(rows_avail/2.0 + rows/2.0);
    int col_start = floor(cols_avail/2.0 - cols/2.0 + 1.0);
    int col_end = floor(cols_avail/2.0 + cols/2.0);
    int slice_start = floor(slices_avail/2.0 - slices/2.0 + 1.0);
    int slice_end = floor(slices_avail/2.0 + slices/2.0);

    for ( int slice = slice_start; slice < slice_end + 1; slice++ )
    {
        for ( int col = col_start; col < col_end + 1; col++ )
        {
            for ( int row = row_start; row < row_end + 1; row++ )
            {
                tarquin::coord fit_spec(row, col, slice); 
                fit_list.push_back(fit_spec);
            }
        }
    }

	tarquin::Options& options = m_session->GetWorkspace().GetOptions();
	options.SetFitList(fit_list);

	QDialog::accept();
}
