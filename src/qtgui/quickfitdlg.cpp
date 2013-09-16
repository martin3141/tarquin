#include "proto_buf.hpp"
#include "session.h"
#include "exception.hpp"
#include "quickfitdlg.h"
#include "voxel_select.h"
#include "guicommon.h"
#include "formatdlg.h"

#include <boost/filesystem.hpp>
#include <QtGui>


quickfitdlg::quickfitdlg(QWidget *parent, Session* session) :
    QDialog(parent),
	m_session(session),
	m_loaded_ws(false)
{
	assert( session );
    m_ui.setupUi(this);
	connect(m_ui.btnOpenWS, SIGNAL(clicked()), this, SLOT(OnBtnOpenWS()));
	connect(m_ui.btnOpenWU, SIGNAL(clicked()), this, SLOT(OnBtnOpenWU()));
}

void quickfitdlg::OnBtnOpenWS()
{
    // reset the data format
	tarquin::Options& opts = m_session->GetWorkspace().GetOptions();
    opts.SetFormat(tarquin::NOTSET);
	// retrieve the path that was last used from the settings
	QSettings settings;

	QString path = QFileDialog::getOpenFileName(this, 
			tr("Open Water Suppressed Data"), 
			settings.value("ws_input_dir", QString()).toString(),
			tr("Files (*)"));

	// nothing selected, they must have pushed cancel
	if( !path.size() )
		return;

	GetPathAndSave(path, "ws_input_dir");

	// set the file on the form
	m_ui.txtWS->setText(path);

	// and load the FID
	LoadFID(path, WATER_SUPPRESSED_FID);
}

void quickfitdlg::OnBtnOpenWU()
{
	if( !m_loaded_ws )
	{
		InfoDialog(this, tr("Water Suppressed Missing"), tr("The water suppressed data must be loaded first."));
		return;
	}

	// retrieve the path that was last used from the settings
	QSettings settings;

	QString path = QFileDialog::getOpenFileName(this, 
			tr("Open Water Unsuppressed Data"), 
			settings.value("ws_input_dir", QString()).toString(),
			tr("Files (*)"));

	// nothing selected, they must have pushed cancel
	if( !path.size() )
		return;

	// set the file on the form
	m_ui.txtWU->setText(path);

	// and load the FID
	LoadFID(path, WATER_UNSUPPRESSED_FID);
}

QString quickfitdlg::GetPathAndSave(
		QString path,          //!< the full filename
		QString settings_name  //!< the name of the settings we are going to save it as 
		)
{
	// turn it into a path
	boost::filesystem::path dir = path.toStdString();

	// get the parent directory as a QString
	QString parent_dir = QString::fromStdString(dir.parent_path().string());

	// write it to the registry
	QSettings settings;
	settings.setValue(settings_name, parent_dir);

	// and return the directory in case they want it
	return QString::fromStdString(dir.string());
}

void quickfitdlg::LoadFID(QString filename, fid_type_e fid_type)
{
	// set the options on the session
	tarquin::Options& opts = m_session->GetWorkspace().GetOptions();
	std::string fn = filename.toStdString();
	tarquin::CFID& fid = m_session->GetWorkspace().GetFID();

	// try and guess the format based on the filename
    bool dcm = true;
	if ( opts.GetFormat() == tarquin::NOTSET )
	{
        // try to load as DICOM
        try
        {
            opts.SetFormat(tarquin::DCM);
            fid.Load(filename.toStdString(), opts, m_session->GetWorkspace(), GetLog());
        }
        catch( const tarquin::Exception& e )
        {
            // remove any bad data fields from the previous load attempt
            tarquin::CFID newfid;
            fid = newfid;

            dcm = false;
            if ( (fn.substr(fn.size()-5, filename.size()) == ".sdat") 
                    || (fn.substr(fn.size()-5, filename.size()) == ".SDAT")
                    || (fn.substr(fn.size()-5, filename.size()) == ".spar")
                    || (fn.substr(fn.size()-5, filename.size()) == ".SPAR" ) )
                opts.SetFormat(tarquin::PHILIPS);
            else if ( fn.substr(fn.size()-4, filename.size()) == ".shf" )
                opts.SetFormat(tarquin::SHF);
            else if ( fn.substr(fn.size()-2, filename.size()) == ".7" )
                opts.SetFormat(tarquin::GE);
            else if ( fn.substr(fn.size()-4, filename.size()) == ".rda"
                    || fn.substr(fn.size()-4, filename.size()) == ".RDA" )
                opts.SetFormat(tarquin::RDA);
            else
                opts.SetFormat(tarquin::SIEMENS);
        }
	}

	try
	{
		if( WATER_SUPPRESSED_FID == fid_type )
		{
            if ( !dcm )
                fid.Load(filename.toStdString(), opts, m_session->GetWorkspace(), GetLog());

			// did the FID contain the water data as well?
			if( fid.GetCWF() )
			{
				m_ui.txtWU->setText(filename);

                // copy some parameters
			    m_session->GetWorkspace().GetFIDWater().SetParametersFromFID(fid);
				m_session->GetWorkspace().GetFIDWater().LoadW(filename.toStdString(), opts, GetLog());

				// disable the button for loading the separate water unsuprressed FID
				m_ui.btnOpenWU->setEnabled(false);
			}

			// flag as loaded
			m_loaded_ws = true;
		}
		// the water reference file is separate
		else
		{
			m_session->GetWorkspace().GetFIDWater().Load(filename.toStdString(), opts, m_session->GetWorkspace(), GetLog());
		}

	}
	catch( const tarquin::Exception& e )
	{
		QMessageBox::information(this, tr("Data format"), tr("The data format was not detected correctly. Please select manually."));
        opts.SetFormat(tarquin::NOTSET);
		formatdlg fmtdlg(this, m_session);
		if( QDialog::Accepted != fmtdlg.exec() )
			return;
		
		// MW - not sure if we need this?
		/*if( WATER_SUPPRESSED_FID == fid_type )
			m_ui.txtWS->setText(QString());
		else
			m_ui.txtWU->setText(QString());
			*/
		
		LoadFID(filename, fid_type);
	}
	m_ui.btnOpenWS->setEnabled(false);
}

void quickfitdlg::accept()
{
	QString water_suppressed   = m_ui.txtWS->text();
	QString water_unsuppressed = m_ui.txtWU->text();

	if( !m_loaded_ws )
	{
		QMessageBox::information(this, tr("Missing FID"), tr("You must provide a path for the water suppressed FID data."));
		return;
	}
	
	tarquin::Options& opts = m_session->GetWorkspace().GetOptions();
	opts.SetFilenameWater(water_unsuppressed.toStdString());
	opts.SetFilename(water_suppressed.toStdString());

	tarquin::CFID& fid = m_session->GetWorkspace().GetFID();

    if ( fid.GetVoxelCount() > 1 )
    {
        voxel_select voxel_select_dlg(this, m_session);
        if( QDialog::Accepted != voxel_select_dlg.exec() )
            return;
    }

    m_session->data_loaded = true;
	QDialog::accept();
}
