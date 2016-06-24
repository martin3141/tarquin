#include "proto_buf.hpp"
#include "session.h"
#include "mainwindow.h"
#include "version/version.h"
#include "aboutdlg.h"
#include "inputdlg.h"
#include "quickfitdlg.h"
#include "tarquinplotwidget.h"
#include "guicommon.h"
#include "fidio/CDICOMFile.hpp"
#include <boost/filesystem.hpp>

#include "export_data.hpp"
#include "preprocess.hpp"
#include "tarquin.hpp"
#include "common.hpp"

#include "MRI.hpp"

#include <QtGui>
#include <qwt_plot.h>
#include <qwt_legend.h>
#include <qwt_plot_picker.h>

MainWindow::MainWindow(QWidget* parent, Qt::WFlags flags) :
	QMainWindow(parent, flags),
	m_session(NULL)
{
	m_ui.setupUi(this);

    QList<int> ratio;
    ratio.append(1);
    ratio.append(0);
    m_ui.splitter_main->setSizes(ratio);
    m_ui.splitter_map->setSizes(ratio);

	// construct the title bar
	std::string title = tarquin::version::display_name() + " " + tarquin::version::version_string();

	//setWindowTitle(  QString::fromStdString(title) ); 

	// connect up actions on the file menu
	connect(m_ui.actionQuick_Fit, SIGNAL(triggered()), this, SLOT(OnFileQF()));
	connect(m_ui.actionOpen_FID,        SIGNAL(triggered()), this, SLOT(OnFileOpenFID()));
	connect(m_ui.actionOpen_fit,        SIGNAL(triggered()), this, SLOT(OnFileOpenFit()));
	connect(m_ui.actionSave_fit,        SIGNAL(triggered()), this, SLOT(OnFileSaveFit()));
	connect(m_ui.actionAbout,           SIGNAL(triggered()), this, SLOT(OnFileAbout()));
	connect(m_ui.actionPrint_Plot,      SIGNAL(triggered()), this, SLOT(OnFilePrint()));
	connect(m_ui.actionPrint_Localisation,      SIGNAL(triggered()), this, SLOT(OnFilePrintLocalisation()));
	connect(m_ui.actionPDF,      SIGNAL(triggered()), this, SLOT(OnFileExportPDF()));
	connect(m_ui.actionPNG,      SIGNAL(triggered()), this, SLOT(OnFileExportPNG()));
	connect(m_ui.actionJPEG,      SIGNAL(triggered()), this, SLOT(OnFileExportJPEG()));
	connect(m_ui.action_localisation_PDF,      SIGNAL(triggered()), this, SLOT(OnFileExportLocalisationPDF()));
	connect(m_ui.action_localisation_PNG,      SIGNAL(triggered()), this, SLOT(OnFileExportLocalisationPNG()));
	connect(m_ui.action_localisation_JPEG,      SIGNAL(triggered()), this, SLOT(OnFileExportLocalisationJPEG()));
	connect(m_ui.actionExport_window,      SIGNAL(triggered()), this, SLOT(OnFileExportWindowPNG()));
	connect(m_ui.actionExit,            SIGNAL(triggered()), this, SLOT(OnFileExit()));

	// connect up actions on the analysis menu
	connect(m_ui.actionRun_TARQUIN,     SIGNAL(triggered()), this, SLOT(OnAnalysisRunTARQUIN()));
	connect(m_ui.actionPreprocess,     SIGNAL(triggered()), this, SLOT(OnAnalysisPreprocess()));

	// connect up actions on the plot menu (formerly called 'view')
	connect(m_ui.actionFID_seconds,			    SIGNAL(triggered()), this, SLOT(OnViewTimeDomainSec()));
	connect(m_ui.actionFID_points,			    SIGNAL(triggered()), this, SLOT(OnViewTimeDomainPts()));
	connect(m_ui.actionSpectrum_PPM,			SIGNAL(triggered()), this, SLOT(OnViewFreqDomainPPM()));
	connect(m_ui.actionSpectrum_Hz,			SIGNAL(triggered()), this, SLOT(OnViewFreqDomainHz()));
	connect(m_ui.actionSpectrum_points,			SIGNAL(triggered()), this, SLOT(OnViewFreqDomainPts()));
	connect(m_ui.actionReal,			SIGNAL(triggered()), this, SLOT(OnViewReal()));
	connect(m_ui.actionImaginary,		SIGNAL(triggered()), this, SLOT(OnViewImag()));
	connect(m_ui.actionAbsolute,		SIGNAL(triggered()), this, SLOT(OnViewAbs()));
	connect(m_ui.actionNext_fit,		SIGNAL(triggered()), this, SLOT(OnNextFit()));
	connect(m_ui.actionPrevious_fit,	SIGNAL(triggered()), this, SLOT(OnPreviousFit()));
	connect(m_ui.actionGotoCenterVoxel,	SIGNAL(triggered()), this, SLOT(OnCenterVoxel()));
	connect(m_ui.actionSubtract_baseline,	    SIGNAL(triggered()), this, SLOT(OnSubtractBaseline()));
	connect(m_ui.actionShow_raw_data,	    SIGNAL(triggered()), this, SLOT(OnShowRawData()));
	connect(m_ui.actionAAToggle,	    SIGNAL(triggered()), this, SLOT(OnAAToggle()));
	connect(m_ui.actionPlot_WUS_data,	    SIGNAL(triggered()), this, SLOT(OnWUSToggle()));
	connect(m_ui.actionSet_line_broadening,	    SIGNAL(triggered()), this, SLOT(OnSetLB()));
	connect(m_ui.actionSet_zero_filling_factor,	    SIGNAL(triggered()), this, SLOT(OnSetZF()));
	connect(m_ui.actionSet_left_ppm_limit,	    SIGNAL(triggered()), this, SLOT(OnSetLeftppm()));
	connect(m_ui.actionSet_right_ppm_limit,	    SIGNAL(triggered()), this, SLOT(OnSetRightppm()));
	connect(m_ui.actionSet_baseline_smoothing,	SIGNAL(triggered()), this, SLOT(OnSetBS()));

	// connect up the actions on the results menu
	connect(m_ui.actViewAmplitudes, SIGNAL(triggered()), this, SLOT(OnViewAmplitudes()));
	connect(m_ui.actExport_PDF,     SIGNAL(triggered()), this, SLOT(OnExportPDF()));
	connect(m_ui.actExport_TXT,     SIGNAL(triggered()), this, SLOT(OnExportTXT()));
	connect(m_ui.actExport_CSV,     SIGNAL(triggered()), this, SLOT(OnExportCSV()));
	connect(m_ui.actExport_CSV_SV,     SIGNAL(triggered()), this, SLOT(OnExportCSV_SV()));
	connect(m_ui.actExport_CSV_fit_SV,     SIGNAL(triggered()), this, SLOT(OnExportCSVFit_SV()));
	connect(m_ui.actExport_CSV_fit,     SIGNAL(triggered()), this, SLOT(OnExportCSVFit()));
	connect(m_ui.actExport_CSV_fit_mag,     SIGNAL(triggered()), this, SLOT(OnExportCSVFitMag()));
	connect(m_ui.actExport_CSV_spectra,     SIGNAL(triggered()), this, SLOT(OnExportCSVSpectra()));
	connect(m_ui.actExport_CSV_spectra_mag,     SIGNAL(triggered()), this, SLOT(OnExportCSVSpectraMag()));
	connect(m_ui.actExport_DPT_raw,     SIGNAL(triggered()), this, SLOT(OnExportDPTRaw()));
	connect(m_ui.actExport_DPT_proc,     SIGNAL(triggered()), this, SLOT(OnExportDPTProc()));
	connect(m_ui.actExport_DPT_water,     SIGNAL(triggered()), this, SLOT(OnExportDPTWater()));
    
    connect(m_ui.actionAuto_scale_x_axis, SIGNAL(triggered()), this, SLOT(OnAutoX()));
    connect(m_ui.actionAuto_scale_y_axis, SIGNAL(triggered()), this, SLOT(OnAutoY()));
    
    connect(m_ui.spinSlices, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));
    connect(m_ui.spinRows, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));
    connect(m_ui.spinCols, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));

    connect(m_ui.spinCols, SIGNAL(editingFinished()), this, SLOT(OnEditFin()));
    connect(m_ui.spinRows, SIGNAL(editingFinished()), this, SLOT(OnEditFin()));
    connect(m_ui.spinSlices, SIGNAL(editingFinished()), this, SLOT(OnEditFin()));

    connect(m_ui.check_box_fit, SIGNAL(clicked(bool)), this, SLOT(OnFitChange()));

    connect(m_ui.sel_none, SIGNAL(pressed()), this, SLOT(OnSelNone()));
    connect(m_ui.sel_all, SIGNAL(pressed()), this, SLOT(OnSelAll()));
    
    connect(m_ui.close_mri, SIGNAL(pressed()), this, SLOT(OnCloseMRI()));
    connect(m_ui.load_mri, SIGNAL(pressed()), this, SLOT(OnLoadMRI()));
    connect(m_ui.ww_wl, SIGNAL(pressed()), this, SLOT(OnWWLL()));
    connect(m_ui.hide_grid, SIGNAL(clicked()), this, SLOT(OnHideGrid()));
    connect(m_ui.hide_lines, SIGNAL(clicked()), this, SLOT(OnHideLines()));
    connect(m_ui.hide_mri, SIGNAL(clicked()), this, SLOT(OnHideMRI()));
    connect(m_ui.set_trans, SIGNAL(clicked()), this, SLOT(OnSetTrans()));
    
    connect(m_ui.phi0, SIGNAL(pressed()), this, SLOT(OnPhi0()));
    connect(m_ui.phi1, SIGNAL(pressed()), this, SLOT(OnPhi1()));
    connect(m_ui.cal_ppm, SIGNAL(pressed()), this, SLOT(OnCalPPM()));
    connect(m_ui.set_pivot, SIGNAL(pressed()), this, SLOT(OnSetPivot()));
    
    // disable spinners and fit check box
    m_ui.spinRows->setEnabled(false);
    m_ui.spinRows->setKeyboardTracking(false);
    m_ui.spinCols->setEnabled(false);
    m_ui.spinCols->setKeyboardTracking(false);
    m_ui.spinSlices->setEnabled(false);
    m_ui.spinSlices->setKeyboardTracking(false);
    m_ui.check_box_fit->setEnabled(false);
    m_ui.sel_all->setEnabled(false);
    m_ui.sel_none->setEnabled(false);

    // disable MRI stuff
    m_ui.ww_wl->setEnabled(false);
    m_ui.load_mri->setEnabled(false);
    m_ui.close_mri->setEnabled(false);
    m_ui.slice_slider->setEnabled(false);
    m_ui.slice_spin->setEnabled(false);
    m_ui.hide_grid->setEnabled(false);
    m_ui.hide_lines->setEnabled(false);
    m_ui.hide_mri->setEnabled(false);
    m_ui.set_trans->setEnabled(false);


	// make the plotting widget
	m_plot = new TarquinPlotWidget(m_ui.plot_frame);
	QHBoxLayout* layout = new QHBoxLayout(m_ui.plot_frame);
	layout->addWidget(m_plot);

    // MRS Geom display 
    m_view = new MyGraphicsView(m_ui.grid_frame);
    m_view->setRenderHints(QPainter::Antialiasing | QPainter::HighQualityAntialiasing | QPainter::TextAntialiasing);
	QHBoxLayout* layout_view = new QHBoxLayout(m_ui.grid_frame);
    layout_view->addWidget(m_view);
    m_view->showMaximized();
    m_scene = new MyGraphicsScene();
    
    m_scene->setBackgroundBrush(Qt::black);

    m_view->scale(1.75,1.75);
    m_view->setScene(m_scene);
    
    m_ui.phi0->setVisible(false);
    m_ui.phi1->setVisible(false);
    m_ui.set_pivot->setVisible(false);
    m_ui.cal_ppm->setVisible(false);

    m_ui.slice_slider->setMinimum(0);
    m_ui.slice_slider->setMaximum(0);
    m_ui.slice_slider->setValue(0);

    m_ui.slice_spin->setMinimum(0);
    m_ui.slice_spin->setMaximum(0);
    m_ui.slice_spin->setValue(0);

    m_phi0 = false;
    m_phi1 = false;
    m_ww_wl = false;
}

void MainWindow::OnExportTXT()
{
    if( !m_session )
    {
    	QMessageBox::information(this, tr("No Data Ready"), 
				tr("You need to generate some results before you can export."));
		return;
    }

	if( !m_session->FitAvailable() )
	{
		QMessageBox::information(this, tr("No Data Ready"), 
				tr("You need to generate some results before you can export."));
		return;
	}

	QString file = QFileDialog::getSaveFileName(this, tr("Export TXT"), "", tr("TXT Files (*.txt)"));

	if( !file.size() )
		return;

    
    // export the data
	try
	{
		ExportTxtResults(file.toStdString(), m_session->GetWorkspace());

		InfoDialog(this, tr("Success"), tr("File was exported to: ") + file);
	}
	catch( const std::exception& e )
	{
		ErrorDialog(this, tr("Error Exporting File"), tr("There was an error: ") + e.what());
	}
}

void MainWindow::OnExportPDF()
{
    if( !m_session )
    {
    	QMessageBox::information(this, tr("No Data Ready"), 
				tr("You need to generate some results before you can export."));
		return;
    }

	if( !m_session->FitAvailable() )
	{
		QMessageBox::information(this, tr("No Data Ready"), 
				tr("You need to generate some results before you can export."));
		return;
	}

    if ( !m_session->m_in_fit_list )
	{
		QMessageBox::information(this, tr("No Data Ready"), 
				tr("This voxel has not been analysed."));
		return;
	}

	QString file = QFileDialog::getSaveFileName(this, tr("Export PDF"), "", tr("PDF Files (*.pdf)"));

	if( !file.size() )
		return;

    m_session->GetWorkspace().GetOptions().SetPdfStack(true);
    
    // export the data
	try
	{
	    int fit_no = m_session->m_fit_number;
		ExportPdfResults(file.toStdString(), m_session->GetWorkspace(),fit_no);

		InfoDialog(this, tr("Success"), tr("File was exported to: ") + file);
	}
	catch( const std::exception& e )
	{
		ErrorDialog(this, tr("Error Exporting File"), tr("There was an error: ") + e.what());
	}
}

void MainWindow::OnExportPDFStack()
{
    if( !m_session )
    {
    	QMessageBox::information(this, tr("No Data Ready"), 
				tr("You need to generate some results before you can export."));
		return;
    }

	if( !m_session->FitAvailable() )
	{
		QMessageBox::information(this, tr("No Data Ready"), 
				tr("You need to generate some results before you can export."));
		return;
	}

    if ( !m_session->m_in_fit_list )
	{
		QMessageBox::information(this, tr("No Data Ready"), 
				tr("This voxel has not been analysed."));
		return;
	}

	QString file = QFileDialog::getSaveFileName(this, tr("Export PDF"), "", tr("PDF Files (*.pdf)"));

	if( !file.size() )
		return;

    
    m_session->GetWorkspace().GetOptions().SetPdfStack(true);
    
    // export the data
	try
	{
	    int fit_no = m_session->m_fit_number;
		ExportPdfResults(file.toStdString(), m_session->GetWorkspace(),fit_no);

		InfoDialog(this, tr("Success"), tr("File was exported to: ") + file);
	}
	catch( const std::exception& e )
	{
		ErrorDialog(this, tr("Error Exporting File"), tr("There was an error: ") + e.what());
	}
}

void MainWindow::OnExportDPTWater()
{
    if( !m_session )
    {
    	QMessageBox::information(this, tr("No Data Ready"), 
				tr("You need to load some data before you can export."));
		return;
    }
    

    if ( m_session->GetWorkspace().GetOptions().GetFilenameWater() == "" && !m_session->GetWorkspace().GetFIDRaw().GetCWF() )
    {
		QMessageBox::information(this, tr("Can't export"), 
				tr("Water data is not avaiable for this voxel."));
		return;
    }

    QString file = QFileDialog::getSaveFileName(this, tr("Export DPT"), "", tr("DPT Files (*.dpt)"));

	if( !file.size() )
		return;

    tarquin::CFID fidwater = m_session->GetWorkspace().GetFIDWater();
    int num = fidwater.vox2ind(m_session->m_voxel);
    fidwater.SaveToFile(file.toStdString(),num);
}

void MainWindow::OnExportDPTRaw()
{
    OnExportDPT(false);
}

void MainWindow::OnExportDPTProc()
{
    OnExportDPT(true);
}

void MainWindow::OnExportDPT(bool processed)
{
    if( !m_session )
    {
    	QMessageBox::information(this, tr("No Data Ready"), 
				tr("You need to load some data before you can export."));
		return;
    }

    if ( !m_session->data_preprocessed && processed )
	{
		QMessageBox::information(this, tr("No Data Ready"), 
				tr("You need to pre-process spectra before you can export to dpt."));
		return;
	}

    // check this voxel is in the fit list
    if ( !m_session->m_in_fit_list && processed )
	{
		QMessageBox::information(this, tr("Can't export"), 
				tr("This voxel has not been pre-processed."));
		return;
	}

	QString file = QFileDialog::getSaveFileName(this, tr("Export DPT"), "", tr("DPT Files (*.dpt)"));

	if( !file.size() )
		return;

    // export the data
	try
	{
        tarquin::Workspace& workspace = m_session->GetWorkspace();

        tarquin::CFID fid = workspace.GetFIDRaw();
        tarquin::CFID fidproc = workspace.GetFIDProc();

        if ( processed )
        {
            // preprocessing and interactive phasing parameters are
            // "hard" applied to fidproc
            fidproc.SetPhi0(m_session->m_voxel,0.0);
            fidproc.SetPhi1(m_session->m_voxel,0.0);
            int num = fidproc.vox2ind(m_session->m_voxel);
            fidproc.SaveToFile(file.toStdString(),num);
        }
        else
        {
            fid.SetPhi0(m_session->m_voxel,fidproc.GetPhi0(m_session->m_voxel));
            fid.SetPhi1(m_session->m_voxel,fidproc.GetPhi1(m_session->m_voxel));
            fid.SetPPMRef(m_session->m_voxel,fidproc.GetPPMRef(m_session->m_voxel));
            int num = fid.vox2ind(m_session->m_voxel);
            fid.SaveToFile(file.toStdString(),num);
        }

		InfoDialog(this, tr("Success"), tr("File was exported to: ") + file);
	}
	catch( const std::exception& e )
	{
		ErrorDialog(this, tr("Error Exporting File"), tr("There was an error: ") + e.what());
	}
}

void MainWindow::OnExportCSV_SV()
{
    if( !m_session )
    {
    	QMessageBox::information(this, tr("No Data Ready"), 
				tr("You need to generate some results before you can export."));
		return;
    }

	if( !m_session->FitAvailable() )
	{
		QMessageBox::information(this, tr("No Data Ready"), 
				tr("You need to generate some results before you can export."));
		return;
	}

    if ( !m_session->m_in_fit_list )
	{
		QMessageBox::information(this, tr("No Data Ready"), 
				tr("This voxel has not been analysed."));
		return;
	}

    
	QString file = QFileDialog::getSaveFileName(this, tr("Export CSV"), "", tr("CSV Files (*.csv)"));

	if( !file.size() )
		return;

    // export the data
	try
	{

	    int num = m_session->m_fit_number;
		ExportCsvResults(file.toStdString(), m_session->GetWorkspace(), num);

		InfoDialog(this, tr("Success"), tr("File was exported to: ") + file);
	}
	catch( const std::exception& e )
	{
		ErrorDialog(this, tr("Error Exporting File"), tr("There was an error: ") + e.what());
	}

}

void MainWindow::OnExportCSV()
{
    if( !m_session )
    {
    	QMessageBox::information(this, tr("No Data Ready"), 
				tr("You need to generate some results before you can export."));
		return;
    }

	if( !m_session->FitAvailable() )
	{
		QMessageBox::information(this, tr("No Data Ready"), 
				tr("You need to generate some results before you can export."));
		return;
	}
    
	QString file = QFileDialog::getSaveFileName(this, tr("Export CSV"), "", tr("CSV Files (*.csv)"));

	if( !file.size() )
		return;

    // export the data
	try
	{
		ExportCsvResults(file.toStdString(), m_session->GetWorkspace());

		InfoDialog(this, tr("Success"), tr("File was exported to: ") + file);
	}
	catch( const std::exception& e )
	{
		ErrorDialog(this, tr("Error Exporting File"), tr("There was an error: ") + e.what());
	}

}

void MainWindow::OnExportCSVSpectraMag()
{
    if( !m_session )
    {
    	QMessageBox::information(this, tr("No Data Ready"), 
				tr("You need to generate some results before you can export."));
		return;
    }

    // check there is at least one voxel to plot
    tarquin::Workspace& workspace = m_session->GetWorkspace();
	tarquin::Options& options = workspace.GetOptions();
	std::vector<tarquin::coord>& fit_list = options.GetFitList();
    if ( fit_list.size() == 0 )
    {
        QMessageBox::information(this, tr("CSV export stopped"), 
			tr("Fit list empty. Please choose some voxels to be fit first."));
        return;
    }

    QString file = QFileDialog::getSaveFileName(this, tr("Export CSV"), "", tr("CSV Files (*.csv)"));

	if( !file.size() )
		return;
    
    // export the data
	try
	{
        ExportCsvSpectraAligned(file.toStdString(), m_session->GetWorkspace(), true);

		InfoDialog(this, tr("Success"), tr("File was exported to: ") + file);
	}
	catch( const std::exception& e )
	{
		ErrorDialog(this, tr("Error Exporting File"), tr("There was an error: ") + e.what());
	}

}

void MainWindow::OnExportCSVSpectra()
{
    if( !m_session )
    {
    	QMessageBox::information(this, tr("No Data Ready"), 
				tr("You need to generate some results before you can export."));
		return;
    }

    // check there is at least one voxel to plot
    tarquin::Workspace& workspace = m_session->GetWorkspace();
	tarquin::Options& options = workspace.GetOptions();
	std::vector<tarquin::coord>& fit_list = options.GetFitList();
    if ( fit_list.size() == 0 )
    {
        QMessageBox::information(this, tr("CSV export stopped"), 
			tr("Fit list empty. Please choose some voxels to be fit first."));
        return;
    }

    QString file = QFileDialog::getSaveFileName(this, tr("Export CSV"), "", tr("CSV Files (*.csv)"));

	if( !file.size() )
		return;
    
    // export the data
	try
	{
        ExportCsvSpectraAligned(file.toStdString(), m_session->GetWorkspace());

		InfoDialog(this, tr("Success"), tr("File was exported to: ") + file);
	}
	catch( const std::exception& e )
	{
		ErrorDialog(this, tr("Error Exporting File"), tr("There was an error: ") + e.what());
	}

}

void MainWindow::OnExportCSVFit_SV()
{
    if( !m_session )
    {
    	QMessageBox::information(this, tr("No Data Ready"), 
				tr("You need to generate some results before you can export."));
		return;
    }
    
    // have we done a fit yet? 
    if( !m_session->data_fitted )
	{
		QMessageBox::information(this, tr("No Data Ready"), 
				tr("You need to generate some results before you can export."));
		return;
	}

    if ( !m_session->m_in_fit_list )
	{
		QMessageBox::information(this, tr("No Data Ready"), 
				tr("This voxel has not been analysed."));
		return;
	}
    
    // check there is at least one voxel to plot
    tarquin::Workspace& workspace = m_session->GetWorkspace();
	tarquin::Options& options = workspace.GetOptions();
	std::vector<tarquin::coord>& fit_list = options.GetFitList();
    if ( fit_list.size() == 0 )
    {
        QMessageBox::information(this, tr("CSV export stopped"), 
			tr("Fit list empty. Please choose some voxels to be fit first."));
        return;
    }

	QString file = QFileDialog::getSaveFileName(this, tr("Export CSV"), "", tr("CSV Files (*.csv)"));

	if( !file.size() )
		return;
    
    // export the data
	try
	{
	    int num = m_session->m_fit_number;
        ExportCsvFit(file.toStdString(), m_session->GetWorkspace(), num);

		InfoDialog(this, tr("Success"), tr("File was exported to: ") + file);
	}
	catch( const std::exception& e )
	{
		ErrorDialog(this, tr("Error Exporting File"), tr("There was an error: ") + e.what());
	}

}


void MainWindow::OnExportCSVFit()
{
    if( !m_session )
    {
    	QMessageBox::information(this, tr("No Data Ready"), 
				tr("You need to generate some results before you can export."));
		return;
    }
    
    // have we done a fit yet? 
    if( !m_session->data_fitted )
	{
		QMessageBox::information(this, tr("No Data Ready"), 
				tr("You need to generate some results before you can export."));
		return;
	}
    
    // check there is at least one voxel to plot
    tarquin::Workspace& workspace = m_session->GetWorkspace();
	tarquin::Options& options = workspace.GetOptions();
	std::vector<tarquin::coord>& fit_list = options.GetFitList();
    if ( fit_list.size() == 0 )
    {
        QMessageBox::information(this, tr("CSV export stopped"), 
			tr("Fit list empty. Please choose some voxels to be fit first."));
        return;
    }

	QString file = QFileDialog::getSaveFileName(this, tr("Export CSV"), "", tr("CSV Files (*.csv)"));

	if( !file.size() )
		return;
    
    // export the data
	try
	{
        ExportCsvFit(file.toStdString(), m_session->GetWorkspace());

		InfoDialog(this, tr("Success"), tr("File was exported to: ") + file);
	}
	catch( const std::exception& e )
	{
		ErrorDialog(this, tr("Error Exporting File"), tr("There was an error: ") + e.what());
	}

}

void MainWindow::OnExportCSVFitMag()
{
    if( !m_session )
    {
    	QMessageBox::information(this, tr("No Data Ready"), 
				tr("You need to generate some results before you can export."));
		return;
    }
    
    // have we done a fit yet? 
    if( !m_session->data_fitted )
	{
		QMessageBox::information(this, tr("No Data Ready"), 
				tr("You need to generate some results before you can export."));
		return;
	}
    
    // check there is at least one voxel to plot
    tarquin::Workspace& workspace = m_session->GetWorkspace();
	tarquin::Options& options = workspace.GetOptions();
	std::vector<tarquin::coord>& fit_list = options.GetFitList();
    if ( fit_list.size() == 0 )
    {
        QMessageBox::information(this, tr("CSV export stopped"), 
			tr("Fit list empty. Please choose some voxels to be fit first."));
        return;
    }

	QString file = QFileDialog::getSaveFileName(this, tr("Export CSV"), "", tr("CSV Files (*.csv)"));

	if( !file.size() )
		return;
    
    // export the data
	try
	{
        ExportCsvFit(file.toStdString(), m_session->GetWorkspace(), -1, true);

		InfoDialog(this, tr("Success"), tr("File was exported to: ") + file);
	}
	catch( const std::exception& e )
	{
		ErrorDialog(this, tr("Error Exporting File"), tr("There was an error: ") + e.what());
	}

}

void MainWindow::OnFilePrint()
{
	if( !m_session ) 
		return;
	
	// the print object
    QPrinter printer;
    printer.setOrientation(QPrinter::Landscape);

	// bring up the setup dialog
	QPrintDialog dlg(&printer);

	if( !dlg.exec() )
		return;

	// and print
	m_plot->print(printer);
}

void MainWindow::OnFilePrintLocalisation()
{
	if( !m_session ) 
		return;
	
	// the print object
    QPrinter printer;
    printer.setOrientation(QPrinter::Landscape);

	// bring up the setup dialog
	QPrintDialog dlg(&printer);

	if( !dlg.exec() )
		return;
    
    QPainter painter(&printer);
    painter.setRenderHints(QPainter::Antialiasing | QPainter::HighQualityAntialiasing | QPainter::TextAntialiasing);

	// and print
    m_view->render(&painter);
}


void MainWindow::OnFileExportPDF()
{
	if( !m_session ) 
		return;

	QString file = QFileDialog::getSaveFileName(this, tr("Export PDF"), "", tr("PDF Files (*.pdf)"));

	if( !file.size() )
		return;
			
	try
	{
		// print to a PDF file
		// reduce the font size first
		QwtLegend* legend = m_plot->legend();
		legend->setFont(QFont("Courier", 6));
		QPrinter printer(QPrinter::HighResolution);
		printer.setOrientation( QPrinter::Landscape );
		printer.setColorMode( QPrinter::Color );
		printer.setOutputFileName(file);
		m_plot->print(printer);
		legend->setFont(QFont("Courier", 8));

		InfoDialog(this, tr("Success"), tr("Saved successfully to: ") + file);
	}
	catch( const std::exception& /*e*/ )
	{
		ErrorDialog(this, tr("Failure"), tr("Failed to export file."));
	}
}

void MainWindow::OnFileExportJPEG()
{
	if( !m_session ) 
		return;

	QString file = QFileDialog::getSaveFileName(this, tr("Export JPEG"), "", tr("JPEG Files (*.jpg)"));

	if( !file.size() )
		return;
			
	try
	{
        QPixmap qPix = QPixmap::grabWidget(m_plot);
        qPix.save(file,"JPEG");
        InfoDialog(this, tr("Success"), tr("Saved successfully to: ") + file);
	}
	catch( const std::exception& /*e*/ )
	{
		ErrorDialog(this, tr("Failure"), tr("Failed to export file."));
	}
}

void MainWindow::OnFileExportPNG()
{
	if( !m_session ) 
		return;

	//QwtLegend* legend = m_plot->legend();
    //legend->contentsWidget()->setVisible(false);

	QString file = QFileDialog::getSaveFileName(this, tr("Export PNG"), "", tr("PNG Files (*.png)"));

	if( !file.size() )
		return;
			
	try
	{
        QPixmap qPix = QPixmap::grabWidget(m_plot);
        qPix.save(file,"PNG");
        InfoDialog(this, tr("Success"), tr("Saved successfully to: ") + file);
	}
	catch( const std::exception& /*e*/ )
	{
		ErrorDialog(this, tr("Failure"), tr("Failed to export file."));
	}
}

void MainWindow::OnFileExportLocalisationPDF()
{
	if( !m_session ) 
		return;

	QString file = QFileDialog::getSaveFileName(this, tr("Export PDF"), "", tr("PDF Files (*.pdf)"));

	if( !file.size() )
		return;
			
	try
	{
		// print to a PDF file
		// reduce the font size first
		QPrinter printer(QPrinter::HighResolution);
		printer.setOrientation( QPrinter::Portrait );
		printer.setColorMode( QPrinter::Color );
		printer.setOutputFileName(file);
        QPainter painter(&printer);
        painter.setRenderHints(QPainter::Antialiasing | QPainter::HighQualityAntialiasing | QPainter::TextAntialiasing);

        m_view->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        m_view->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);

        // and print
        m_view->render(&painter);
        m_view->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
        m_view->setVerticalScrollBarPolicy(Qt::ScrollBarAsNeeded);

		InfoDialog(this, tr("Success"), tr("Saved successfully to: ") + file);
	}
	catch( const std::exception& /*e*/ )
	{
		ErrorDialog(this, tr("Failure"), tr("Failed to export file."));
	}
}


void MainWindow::OnFileExportLocalisationPNG()
{
	if( !m_session ) 
		return;

	//QwtLegend* legend = m_plot->legend();
    //legend->contentsWidget()->setVisible(false);

	QString file = QFileDialog::getSaveFileName(this, tr("Export PNG"), "", tr("PNG Files (*.png)"));

	if( !file.size() )
		return;
			
	try
	{
        m_view->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        m_view->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        QPixmap qPix = QPixmap::grabWidget(m_view);
        m_view->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
        m_view->setVerticalScrollBarPolicy(Qt::ScrollBarAsNeeded);
        qPix.save(file,"PNG");
        InfoDialog(this, tr("Success"), tr("Saved successfully to: ") + file);
	}
	catch( const std::exception& /*e*/ )
	{
		ErrorDialog(this, tr("Failure"), tr("Failed to export file."));
	}
}

void MainWindow::OnFileExportWindowPNG()
{
	if( !m_session ) 
		return;

	//QwtLegend* legend = m_plot->legend();
    //legend->contentsWidget()->setVisible(false);

	QString file = QFileDialog::getSaveFileName(this, tr("Export PNG"), "", tr("PNG Files (*.png)"));

	if( !file.size() )
		return;
			
	try
	{
        m_view->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        m_view->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        QPixmap qPix = QPixmap::grabWidget(this);
        m_view->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
        m_view->setVerticalScrollBarPolicy(Qt::ScrollBarAsNeeded);
        qPix.save(file,"PNG");
        InfoDialog(this, tr("Success"), tr("Saved successfully to: ") + file);
	}
	catch( const std::exception& /*e*/ )
	{
		ErrorDialog(this, tr("Failure"), tr("Failed to export file."));
	}
}



void MainWindow::OnFileExportLocalisationJPEG()
{
	if( !m_session ) 
		return;

	//QwtLegend* legend = m_plot->legend();
    //legend->contentsWidget()->setVisible(false);

	QString file = QFileDialog::getSaveFileName(this, tr("Export JPEG"), "", tr("JPEG Files (*.jpg)"));

	if( !file.size() )
		return;
			
	try
	{
        m_view->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        m_view->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        QPixmap qPix = QPixmap::grabWidget(m_view);
        m_view->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
        m_view->setVerticalScrollBarPolicy(Qt::ScrollBarAsNeeded);
        qPix.save(file,"JPEG");
        InfoDialog(this, tr("Success"), tr("Saved successfully to: ") + file);
	}
	catch( const std::exception& /*e*/ )
	{
		ErrorDialog(this, tr("Failure"), tr("Failed to export file."));
	}
}



/*!
 * Make a new session and prompt the user if session already open.
 */
Session* MainWindow::MakeNewSession()
{
    // if an existing session exists, ask the user if they are
	// sure they want to start again
	if( m_session ) 
	{
		int ret = QMessageBox::question(this, tr("Are you sure?"),
				tr("Do you want to abandon your existing session?"),
				QMessageBox::Yes|QMessageBox::No);

		if( ret != QMessageBox::Yes )
			return NULL;
	
		// they are sure
		delete m_session;
		m_session = NULL;
        
        // clear these
        disconnect(m_ui.listWidget_numer, SIGNAL(itemChanged(QListWidgetItem *)), this, SLOT(OnListChange()));
        disconnect(m_ui.listWidget_denom, SIGNAL(itemChanged(QListWidgetItem *)), this, SLOT(OnListChange()));
        m_ui.listWidget_numer->clear();
        m_ui.listWidget_denom->clear();

	    m_plot->clear();
        qDeleteAll(m_scene->items());

        // disable spinners
        m_ui.spinRows->setEnabled(false);
        m_ui.spinCols->setEnabled(false);
        m_ui.spinSlices->setEnabled(false);
        m_ui.check_box_fit->setEnabled(false);
        m_ui.sel_all->setEnabled(false);
        m_ui.sel_none->setEnabled(false);

        // disable MRI stuff
        m_ui.ww_wl->setEnabled(false);
        m_ui.load_mri->setEnabled(false);
        m_ui.close_mri->setEnabled(false);
        m_ui.slice_slider->setEnabled(false);
        m_ui.slice_spin->setEnabled(false);
        m_ui.hide_grid->setEnabled(false);
        m_ui.hide_lines->setEnabled(false);
        m_ui.hide_mri->setEnabled(false);
        m_ui.set_trans->setEnabled(false);
	}
	m_session = new Session(this);

	m_session->m_fit_number = 0;
	
	// set the default domain to be frequency with units of PPM
	m_session->m_show_flags.domain       = Session::FREQUENCY_DOMAIN;
	m_session->m_show_flags.units_fd     = Session::PPM;

    m_ui.phi0->setVisible(false);
    m_ui.phi1->setVisible(false);
    m_ui.set_pivot->setVisible(false);
    m_ui.cal_ppm->setVisible(false);

	return m_session;
}

Session* MainWindow::KillSession()
{
    if( m_session ) 
	{
		// they are sure
		delete m_session;
		m_session = NULL;
    }

	return NULL;
}

void MainWindow::OnFileExit()
{
	int ret = QMessageBox::information(this, tr("Are you sure?"), tr("Are you sure you want to exit?"),
			QMessageBox::Yes|QMessageBox::Cancel);

	if( QMessageBox::Yes == ret )
		QCoreApplication::exit(0);
}

void MainWindow::OnFileOpenFID()
{
	// make a new session, returning it
	if( !MakeNewSession() )
		return;

	InputDlg dlg(this, m_session);

	if( QDialog::Accepted != dlg.exec() )
    {
        KillSession();
		return;
    }

    // looks like the data loaded ok
    // set current voxel and spinners to be in the middle of the grid

    int start_row = 1;
    int start_col = 1;
    int start_slice = 1;

    disconnect(m_ui.spinSlices, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));
    disconnect(m_ui.spinRows, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));
    disconnect(m_ui.spinCols, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));

    m_ui.spinRows->setValue(start_row);
    m_ui.spinCols->setValue(start_col);
    m_ui.spinSlices->setValue(start_slice);

    connect(m_ui.spinSlices, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));
    connect(m_ui.spinRows, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));
    connect(m_ui.spinCols, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));

    coord current(start_row, start_col, start_slice); 
	m_session->SetVoxel(current);

    // set correct max vals
    m_ui.spinRows->setMaximum(m_session->GetWorkspace().GetFID().GetRows());
    m_ui.spinCols->setMaximum(m_session->GetWorkspace().GetFID().GetCols());
    m_ui.spinSlices->setMaximum(m_session->GetWorkspace().GetFID().GetSlices());
       
    // enable spinners
    m_ui.spinRows->setEnabled(true);
    m_ui.spinCols->setEnabled(true);
    m_ui.spinSlices->setEnabled(true);
    m_ui.check_box_fit->setEnabled(true);
    m_ui.sel_all->setEnabled(true);
    m_ui.sel_none->setEnabled(true);

    // enable MRI stuff
    m_ui.ww_wl->setEnabled(true);
    m_ui.load_mri->setEnabled(true);
    m_ui.close_mri->setEnabled(true);
    m_ui.slice_slider->setEnabled(true);
    m_ui.slice_spin->setEnabled(true);
    m_ui.hide_grid->setEnabled(true);
    m_ui.hide_lines->setEnabled(true);
    m_ui.hide_mri->setEnabled(true);
    m_ui.set_trans->setEnabled(true);

	// set the window title
	SetWindowTitle(QString::fromStdString(m_session->GetWorkspace().GetFID().GetFilename()));

    UpdateGeom();
    
    int rows = m_session->GetWorkspace().GetFID().GetRows();
    int cols = m_session->GetWorkspace().GetFID().GetCols();
    int slices = m_session->GetWorkspace().GetFID().GetSlices();
        
    if ( ( rows == 1 ) && ( cols == 1 ) && ( slices == 1 ) )
    {
        QList<int> ratio;
        ratio.append(1);
        ratio.append(0);
        m_ui.splitter_main->setSizes(ratio);
    }
    else
    {
        QList<int> ratio;
        ratio.append(1);
        ratio.append(1);
        m_ui.splitter_main->setSizes(ratio);
    }

    /*QList<int> ratio;
    ratio.append(1);
    ratio.append(1);
    m_ui.splitter_main->setSizes(ratio);*/

    OnVoxelChange();
    //OnEditFin();

	// do the want to start fitting as well?
	if( dlg.startFit() )
	{
		OnAnalysisRunTARQUIN();
	}
	// no, just load the data
	else
	{
		// initial display of the data
		//m_session->m_show_flags.show_input = true;
        m_session->Update();
    }
    OnAutoX();
    OnAutoY();
}

void MainWindow::SetWindowTitle(QString filename)
{
	QString title = tr("TARQUIN - ");
	title += filename;

	this->setWindowTitle(title);
}

/*void MainWindow::OnFileExportFID()
{
}*/

void MainWindow::OnFileOpenFit()
{
	if( !MakeNewSession() )
		return;

	// bring up a dialog
	QString file = QFileDialog::getOpenFileName(this, tr("Load Workspace"), "", tr("XML Files (*.xml)"));

	if( !file.size() )
		return;
			
	try
	{
		m_session->Load(file);
    }
	catch( const std::exception& /*e*/ )
	{
		ErrorDialog(this, tr("Failure"), tr("Failed to open file."));
	}
        
    
    // looks like the data loaded ok
    // set current voxel and spinners to be in the middle of the grid
    int start_row = 1;
    int start_col = 1;
    int start_slice = 1;

    disconnect(m_ui.spinSlices, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));
    disconnect(m_ui.spinRows, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));
    disconnect(m_ui.spinCols, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));

    m_ui.spinRows->setValue(start_row);
    m_ui.spinCols->setValue(start_col);
    m_ui.spinSlices->setValue(start_slice);

    connect(m_ui.spinSlices, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));
    connect(m_ui.spinRows, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));
    connect(m_ui.spinCols, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));

    coord current(start_row, start_col, start_slice); 
	m_session->SetVoxel(current);


    // set correct max vals
    m_ui.spinRows->setMaximum(m_session->GetWorkspace().GetFID().GetRows());
    m_ui.spinCols->setMaximum(m_session->GetWorkspace().GetFID().GetCols());
    m_ui.spinSlices->setMaximum(m_session->GetWorkspace().GetFID().GetSlices());

    OnVoxelChange();
    //OnEditFin();

    // enable spinners
    m_ui.spinRows->setEnabled(true);
    m_ui.spinCols->setEnabled(true);
    m_ui.spinSlices->setEnabled(true);

    SetWindowTitle(QString::fromStdString(m_session->GetWorkspace().GetFID().GetFilename()));

}

void MainWindow::OnFileSaveFit()
{
	if( !m_session )
		return;

	// bring up a dialog
	QString file = QFileDialog::getSaveFileName(this, tr("Save Workspace"), "", tr("XML Files (*.xml)"));

	if( !file.size() )
		return;
			
	try
	{
		m_session->Save(file);
		InfoDialog(this, tr("Success"), tr("Saved successfully to: ") + file);
	}
	catch( const std::exception& /*e*/ )
	{
		ErrorDialog(this, tr("Failure"), tr("Failed to save file."));
	}
}

void MainWindow::OnFileAbout()
{
	AboutDlg dlg(this);
	dlg.exec();
}

void MainWindow::OnFileQF()
{
	// make a new session, returning it
	if( !MakeNewSession() )
		return;

	quickfitdlg dlg(this, m_session);
	
	if( QDialog::Accepted != dlg.exec() )
    {
        KillSession();
		return;
    }
    
    // looks like the data loaded ok
    // set current voxel and spinners to be in the middle of the grid
    int start_row = 1;
    int start_col = 1;
    int start_slice = 1;

    disconnect(m_ui.spinSlices, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));
    disconnect(m_ui.spinRows, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));
    disconnect(m_ui.spinCols, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));

    m_ui.spinRows->setValue(start_row);
    m_ui.spinCols->setValue(start_col);
    m_ui.spinSlices->setValue(start_slice);

    connect(m_ui.spinSlices, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));
    connect(m_ui.spinRows, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));
    connect(m_ui.spinCols, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));

    coord current(start_row, start_col, start_slice); 
	m_session->SetVoxel(current);

    // set correct max vals
    m_ui.spinRows->setMaximum(m_session->GetWorkspace().GetFID().GetRows());
    m_ui.spinCols->setMaximum(m_session->GetWorkspace().GetFID().GetCols());
    m_ui.spinSlices->setMaximum(m_session->GetWorkspace().GetFID().GetSlices());
    
    OnVoxelChange();
    //OnEditFin();

    // enable spinners
    m_ui.spinRows->setEnabled(true);
    m_ui.spinCols->setEnabled(true);
    m_ui.spinSlices->setEnabled(true);

	// initial display of the data
	//m_session->m_show_flags.show_input = true;
	m_session->Update();

	// set the window title
	SetWindowTitle(QString::fromStdString(m_session->GetWorkspace().GetFID().GetFilename()));
	
	// run tarquin
	OnAnalysisRunTARQUIN();
}

/*!
 * Run the preprocessor. This function exists so that the preprocessor can be run in a separate thread
 * so that we can show a progress dialog while it is running.
 */
bool MainWindow::run_preprocessor()
{
	try
	{
		tarquin::Workspace& workspace = m_session->GetWorkspace();

		tarquin::Preprocessor preprocess(workspace, GetLog());

		preprocess();
        
        
		return true;
	}
	catch( const std::exception& e )
	{
		GetLog().LogMessage(tarquin::LOG_ERROR, "Error: '%s'", e.what());
	}

	// if we got here then error
	return false;
}

bool MainWindow::simulate_basis()
{
	tarquin::Workspace& workspace = m_session->GetWorkspace();

	// this is the FID after it has been preprocessed
	tarquin::CFID& fidproc = workspace.GetFIDProc();

	// get the basis
	tarquin::CBasis& basis = workspace.GetBasis();
	
	// convenience reference to algorithm options
	tarquin::Options& options = workspace.GetOptions();

	
	// simulated a basis from a csv dir 
	if ( options.GetBasisPath() != "" )
	{
		if( false == basis.Simulate(options.GetBasisPath(), fidproc, options, GetLog()) ) 
			return false;
	}
	else
	{ 
		// simulate an internal basis
		if( false == basis.Simulate(fidproc, options, GetLog()) ) 
			return false;
	}

	return true;
}

bool MainWindow::run_tarquin()
{
	// convenience reference to the workspace
	tarquin::Workspace& workspace = m_session->GetWorkspace();

	// go ahead and run TARQUIN 
	tarquin::RunTARQUIN(workspace, GetLog());

	return true;
}

void MainWindow::OnAnalysisPreprocess()
{
    // check a session exists
	if( !m_session )
    {
        ErrorDialog(this, tr("Fitting failed"), 
				tr("Please load data first."));
		return;
    }

    if ( !m_session->data_loaded )
    {
        ErrorDialog(this, tr("Fitting failed"), 
				tr("Please load data first."));
		return;
    }
    
    // check we're not repeating the processing 
    if ( m_session->data_preprocessed )
    {
        ErrorDialog(this, tr("Preprocessing stopped"), 
			tr("Data has already been preprocessed."));
        return;
    }
    
    // check there is at least one voxel to be fit
    tarquin::Workspace& workspace = m_session->GetWorkspace();
	tarquin::Options& options = workspace.GetOptions();
	std::vector<tarquin::coord>& fit_list = options.GetFitList();

    if ( fit_list.size() == 0 )
    {
        ErrorDialog(this, tr("Preprocessing stopped"), 
			tr("Fit list empty. Please choose some voxels to be fit first."));
        return;
    }

	// setup the progress dialog that will be displayed whilst we are processing
	QProgressDialog dialog(NULL, Qt::CustomizeWindowHint|Qt::WindowTitleHint);
	dialog.setWindowTitle(QString::fromStdString(tarquin::version::display_name()));
	dialog.setLabelText(QString("Running preprocessor..."));
	dialog.setCancelButton(NULL);
    dialog.setWindowModality(Qt::ApplicationModal);

	// the future watcher that handles interaction while the processing is happening
	QFutureWatcher<bool> watcher;
	QObject::connect(&watcher, SIGNAL(finished()),                    &dialog,   SLOT(reset()));
	QObject::connect(&watcher, SIGNAL(progressRangeChanged(int,int)), &dialog,   SLOT(setRange(int,int)));
	//QObject::connect(&watcher, SIGNAL(progressValueChanged(int)),     &dialog,   SLOT(setValue(int)));

	// start processing
	watcher.setFuture(QtConcurrent::run(this, &MainWindow::run_preprocessor));
	dialog.exec();
	watcher.waitForFinished();

	// make sure that step succeeded
	if( !watcher.result() )
	{
		ErrorDialog(this, tr("Preprocessor Failed"), 
				tr("The preprocessing step failed. Please check the log for details."));
		return;
	}
        
    m_ui.phi0->setVisible(true);
    m_ui.phi1->setVisible(true);
    m_ui.set_pivot->setVisible(true);
    m_ui.cal_ppm->setVisible(true);

    m_session->data_preprocessed = true;

    // disable the fit box
    m_ui.check_box_fit->setEnabled(false);
    m_ui.sel_all->setEnabled(false);
    m_ui.sel_none->setEnabled(false);

	// plot the preprocessed signal
	m_session->m_show_flags.show_preproc = true;
    
    m_ui.phi0->setVisible(true);
    m_ui.phi1->setVisible(true);
    m_ui.set_pivot->setVisible(true);
    m_ui.cal_ppm->setVisible(true);


	// update the plot
	//m_session->Update();
    OnVoxelChange();

}

/*!
 * Run the TARQUIN algorithm, in all its glory, bringing up various progress dialogs to ease the
 * pain of waiting.
 */
void MainWindow::OnAnalysisRunTARQUIN()
{
    // check a session exists
    if( !m_session )
    {
        ErrorDialog(this, tr("Fitting failed"), 
				tr("Please load data first."));
		return;
    }

    // preprocess the data first if needed
    if ( !m_session->data_preprocessed )
        OnAnalysisPreprocess();

    // check we're not repeating the fit
    if ( m_session->data_fitted )
    {
        ErrorDialog(this, tr("Fitting stopped"), 
			tr("Data has already been fitted."));
        return;
    }

	QProgressDialog dialog(NULL, Qt::CustomizeWindowHint|Qt::WindowTitleHint);
    dialog.setWindowTitle(QString::fromStdString(tarquin::version::display_name()));
	dialog.setLabelText(QString("Running preprocessor..."));
	dialog.setCancelButton(NULL);
    dialog.setWindowModality(Qt::ApplicationModal);

	// the future watcher that handles interaction while the processing is happening
	QFutureWatcher<bool> watcher;
	QObject::connect(&watcher, SIGNAL(finished()),                    &dialog,   SLOT(reset()));
	QObject::connect(&watcher, SIGNAL(progressRangeChanged(int,int)), &dialog,   SLOT(setRange(int,int)));
	//QObject::connect(&watcher, SIGNAL(progressValueChanged(int)),     &dialog,   SLOT(setValue(int)));


	// simulate a basis set if required
	tarquin::Options& opts = m_session->GetWorkspace().GetOptions();
	tarquin::CBasis& basis = m_session->GetWorkspace().GetBasis();
    tarquin::CFID& fidproc = m_session->GetWorkspace().GetFIDProc();

	if ( !opts.GetUsePrecompiled() )
	{
		dialog.setLabelText(QString("Simulating basis..."));
		watcher.setFuture(QtConcurrent::run(this, &MainWindow::simulate_basis));
		dialog.exec();
		watcher.waitForFinished();

		// make sure that step succeeded
		if( !watcher.result() )
		{
			ErrorDialog(this, tr("Basis Simulation Failed"), 
					tr("The simulation step failed. Check the csv directory for bad files."));
			return;
		}

		// save the basis if specified
		if ( opts.GetBasisSaveFile() != "")
		{
			basis.Initialise(4.65f);
			if( !save_basis(opts.GetBasisSaveFile().c_str(), basis) )
				throw std::runtime_error("The save_basis function generally failed.");
		}
	}

    // save as .basis format
    if( opts.GetBasisSaveFileLCM() != "" ) 
    {

        if( opts.GetUsePrecompiled() ) 
        {
            basis.SaveLCM(opts.GetBasisSaveFileLCM().c_str(),fidproc);
        }
        else
        {
            CFID lcmfid = fidproc;
            CBasis lcm_basis;
            lcmfid.SetNumberOfPoints(opts.GetLCMBasisPts());

            // resimulate with different N before writing file
            if ( opts.GetBasisPath() != "" )
            {
                if( false == lcm_basis.Simulate(opts.GetBasisPath(), lcmfid, opts, GetLog()) ) 
                {
                    GetLog().Out(LOG_ERROR) << "simulation failed, no basis available so quitting";
                    return;
                }
            }
            else
            {
                // internal basis set
                if( false == lcm_basis.Simulate(lcmfid, opts, GetLog()) ) 
                {
                    GetLog().Out(LOG_ERROR) << "internal simulation failed, no basis available so quitting";
                    return;
                }
            }
            lcm_basis.SaveLCM(opts.GetBasisSaveFileLCM().c_str(),lcmfid);
        }
    }


	// start tarquin
	dialog.setLabelText(QString("Fitting..."));
	watcher.setFuture(QtConcurrent::run(this, &MainWindow::run_tarquin));
	dialog.exec();
	watcher.waitForFinished();

	// hide the input signal 
	//m_session->m_show_flags.show_input = false;
	
	// show the model
	m_session->m_show_flags.show_model = true;

	// add all the vectors to display
	m_session->m_show_flags.basis_indices.clear();
	const cmat_stdvec& S_vec = m_session->GetWorkspace().GetBasisMatrix();
	const cvm::cmatrix& S = S_vec[m_session->m_fit_number];
	for( int i = 0; i < S.nsize(); ++i )
		m_session->m_show_flags.basis_indices.push_back(i);

	// update the plot
	m_session->Update();

    m_session->data_fitted = true;

    m_ui.phi0->setVisible(false);
    m_ui.phi1->setVisible(false);
    m_ui.set_pivot->setVisible(false);
    m_ui.cal_ppm->setVisible(false);

	// bring up the results dialog
	m_session->ShowResults();
    
    //UpdateGeom();
    //OnVoxelChange();
    
    // create tick list for metabolite maps

    // clear previous
    m_ui.listWidget_numer->clear();
    m_ui.listWidget_denom->clear();

    QStringList  itemLabels;
	std::vector < std::string > signal_names = basis.GetSignalNames();	

    for ( size_t n = 0; n < signal_names.size(); n++ )
        itemLabels.append(QString::fromStdString(signal_names[n]));

    QList<QString>::iterator it;
    for (it = itemLabels.begin(); it != itemLabels.end(); ++it)
    {
        QListWidgetItem *listItem_numer = new QListWidgetItem((*it),m_ui.listWidget_numer);
        if ( (*it) == "NAA" || (*it) == "NAAG" || (*it) == "TNAA" )
            listItem_numer->setCheckState(Qt::Checked);
        else
            listItem_numer->setCheckState(Qt::Unchecked);

        m_ui.listWidget_numer->addItem(listItem_numer);

        QListWidgetItem *listItem_denom = new QListWidgetItem((*it),m_ui.listWidget_denom);
        listItem_denom->setCheckState(Qt::Unchecked);
        m_ui.listWidget_denom->addItem(listItem_denom);
    }

    connect(m_ui.listWidget_numer, SIGNAL(itemChanged(QListWidgetItem *)), this, SLOT(OnListChange()));
    connect(m_ui.listWidget_denom, SIGNAL(itemChanged(QListWidgetItem *)), this, SLOT(OnListChange()));

    QList<int> ratio;
    ratio.append(1);
    ratio.append(1);
    m_ui.splitter_map->setSizes(ratio);


    OnListChange();

}

void MainWindow::OnListChange()
{
    std::vector<int> numer_list;
    std::vector<int> denom_list;

    int count = m_ui.listWidget_numer->count();
    for (int index = 0; index < count; index++)
    {
        QListWidgetItem * item_numer = m_ui.listWidget_numer->item(index);
        QListWidgetItem * item_denom = m_ui.listWidget_denom->item(index);

        if ( item_numer->checkState() == Qt::Checked )
            numer_list.push_back(index);

        if ( item_denom->checkState() == Qt::Checked )
            denom_list.push_back(index);
    }

	tarquin::Workspace& workspace = m_session->GetWorkspace();
	const rvec_stdvec ahat  = workspace.GetAmplitudesNormalised();
    
    cvm::rvector new_map(ahat.size());

    bool inf_warning = false;

    for ( size_t n = 0; n < ahat.size(); n++ )
    {
        double numer = 0;
        double denom = 0;
        
        if ( numer_list.size() == 0 )
        {
            numer = 1;
        }
        else
        {
            for ( size_t ind_numer = 0; ind_numer < numer_list.size(); ind_numer++ )
                numer += ahat[n](numer_list[ind_numer]+1);
        }
        
        if ( denom_list.size() == 0 )
        {
            denom = 1;
        }
        else
        {
            for ( size_t ind_denom = 0; ind_denom < denom_list.size(); ind_denom++ )
                denom += ahat[n](denom_list[ind_denom]+1);
        }
        new_map(n+1) = numer/denom;

        if ( denom == 0 )
            inf_warning = true;
    }

    //std::cout << new_map << std::endl;
    
    // fix any inf elements by setting to max finite value
    if ( inf_warning )
    {
        double max_val = -std::numeric_limits<treal>::infinity();
        for ( int n = 1; n < new_map.size()+1; n++ )
            if ( new_map(n) == std::numeric_limits<treal>::infinity() )
                new_map(n) = -std::numeric_limits<treal>::infinity();
            else
                if (new_map(n) > max_val )
                    max_val = new_map(n);

        for ( int n = 1; n < new_map.size()+1; n++ )
            if ( new_map(n) == -std::numeric_limits<treal>::infinity() )
                new_map(n) = max_val;
    }

    //std::cout << new_map << std::endl;

    m_session->map.resize(new_map.size());
    m_session->map = new_map;
    
    if ( inf_warning )
        QMessageBox::information(this, tr("Warning."), tr("Warning : one or more of the denominator values is zero."));


    UpdateGeom();
    OnVoxelChange();
}

void MainWindow::OnViewAmplitudes()
{
	if( !m_session )
    {
    	QMessageBox::information(this, tr("No Data Ready"), 
				tr("You need to generate some results before you can view the amplitudes."));
		return;
    }

	if( !m_session->FitAvailable() )
	{
		QMessageBox::information(this, tr("No Data Ready"), 
				tr("You need to generate some results before you can view the amplitudes."));
		return;
	}

	// got here, we can show them
	m_session->ShowResults();
}

void MainWindow::OnViewReal()
{
	if( !m_session )
		return;

	m_session->m_show_flags.mode = Session::REAL;
	m_session->Update();
}

void MainWindow::OnAutoX()
{
	if( !m_session )
		return;

    m_plot->setAxisAutoScale(QwtPlot::xBottom);
    //m_plot->setAxisScale(QwtPlot::xBottom,4.0,0.5,0);
	m_session->Update();
}

void MainWindow::OnAutoY()
{
	if( !m_session )
		return;

    m_plot->setAxisAutoScale(QwtPlot::yLeft);
	m_session->Update();
}

void MainWindow::OnPhi0()
{
    if( !m_session )
		return;

    m_phi0 = true;
    m_phi1 = false;
    QPoint pos = QCursor::pos();
    m_start_pos_y = pos.y();

    qApp->installEventFilter(this);
}

void MainWindow::OnPhi1()
{
    if( !m_session )
		return;

    m_phi0 = false;
    m_phi1 = true;
    QPoint pos = QCursor::pos();
    m_start_pos_y = pos.y();

    qApp->installEventFilter(this);
}

void MainWindow::OnCalPPM()
{
    // go to frequency domain PPM
	m_session->m_show_flags.domain   = Session::FREQUENCY_DOMAIN;
	m_session->m_show_flags.units_fd = Session::PPM;
    m_session->Update(false);

    // disable the zoomer 
    m_plot->DisableZoomer();

    connect(m_plot->m_picker, SIGNAL(selected(const QwtDoublePoint &)), SLOT(PPMSelected(const QwtDoublePoint &)));
}

void MainWindow::OnSetPivot()
{
    // go to frequency domain PPM
	m_session->m_show_flags.domain   = Session::FREQUENCY_DOMAIN;
	m_session->m_show_flags.units_fd = Session::PPM;
    m_session->Update(false);

    // disable the zoomer 
    m_plot->DisableZoomer();

    connect(m_plot->m_picker, SIGNAL(selected(const QwtDoublePoint &)), SLOT(PivotSelected(const QwtDoublePoint &)));
}


void MainWindow::PPMSelected(const QwtDoublePoint &point)
{
    bool ok;
    double d = QInputDialog::getDouble(this, tr("Set calibration point"), tr("Frequency (PPM):"), point.x(), -std::numeric_limits<treal>::infinity(), std::numeric_limits<treal>::infinity(), 6, &ok);
    
    if (ok)
    {
        tarquin::Workspace& workspace = m_session->GetWorkspace();
        tarquin::CFID& fid = workspace.GetFIDProc();

        fid.SetPPMRef(m_session->m_voxel, fid.GetPPMRef(m_session->m_voxel) - point.x() + d);

        disconnect(m_plot->m_picker,SIGNAL(selected(const QwtDoublePoint &)),0,0);
        m_session->Update(false);
        m_plot->EnableZoomer();
    }
}

void MainWindow::PivotSelected(const QwtDoublePoint &point)
{
    // convert ppm point to hz
    tarquin::Workspace& workspace = m_session->GetWorkspace();
    tarquin::CFID& fid = workspace.GetFIDProc();
    double ppm2hz = fid.GetTransmitterFrequency()/1e6;
    double freq_hz = ( - point.x() + fid.GetPPMRef(m_session->m_voxel) ) * ppm2hz;
	m_session->m_show_flags.pivot   = freq_hz;

    disconnect(m_plot->m_picker,SIGNAL(selected(const QwtDoublePoint &)),0,0);
	m_session->Update(false);
    m_plot->EnableZoomer();
}


bool MainWindow::eventFilter(QObject *obj, QEvent *event)
{
    if (event->type() == QEvent::MouseMove && m_ww_wl && m_ui.ww_wl->isChecked() )
    {
        QPoint pos = QCursor::pos();
        double delta_y = m_start_pos_y - pos.y();
        double delta_x = m_start_pos_x - pos.x();
        m_start_pos_y = pos.y();
        m_start_pos_x = pos.x();
        
        double range = m_session->image.get_max_val() - m_session->image.get_min_val();

        double new_ww = m_session->image.get_ww() - delta_x/range*5000.0;
        double new_wc = m_session->image.get_wc() + delta_y/range*5000.0;

        //std::cout << "New ww : " << new_ww << std::endl;
        //std::cout << "New wc : " << new_wc << std::endl;

        m_session->image.set_ww(new_ww);
        m_session->image.set_wc(new_wc);

        m_session->image.generate_slice(m_ui.slice_slider->value());

        UpdateGeom();
        OnVoxelChange(false);
    }

    if (event->type() == QEvent::MouseButtonRelease && !m_ui.ww_wl->isChecked())
    {
        m_phi0 = false;
        m_phi1 = false;
        qApp->removeEventFilter(this);
    }
    else if (event->type() == QEvent::MouseButtonRelease && m_ui.ww_wl->isChecked()) 
    {
        qApp->removeEventFilter(this);
        m_session->image.generate_slices();
        m_ww_wl = false;
    }
    else if (event->type() == QEvent::MouseMove && m_phi0 && !m_phi1 && !m_ui.ww_wl->isChecked())
    {
        //statusBar()->showMessage(QString("Mouse move (%1,%2)").arg(mouseEvent->pos().x()).arg(mouseEvent->pos().y()));

        tarquin::Workspace& workspace = m_session->GetWorkspace();
        tarquin::CFID& fid = workspace.GetFIDProc();
        cvm::cvector&  y   = fid.GetVectorFID(m_session->m_voxel);

        QPoint pos = QCursor::pos();
        double dphi0 = ( m_start_pos_y - pos.y())/100.0;
        y = y * exp( tcomplex(0, dphi0 ) );

        fid.SetPhi0(m_session->m_voxel, fid.GetPhi0(m_session->m_voxel)+dphi0);

        m_start_pos_y = pos.y();
	    m_session->Update(false);
    }
    else if (event->type() == QEvent::MouseMove && !m_phi0 && m_phi1 && !m_ui.ww_wl->isChecked())
    {
        tarquin::Workspace& workspace = m_session->GetWorkspace();
	    tarquin::Options& options = workspace.GetOptions();
        tarquin::CFID& fid = workspace.GetFIDProc();
        //std::cout << m_session->m_voxel.row << "," << m_session->m_voxel.col << "," << m_session->m_voxel.slice << std::endl;
        cvm::cvector&  y   = fid.GetVectorFID(m_session->m_voxel);
        cvm::cvector Y = fft(y);
        Y = fftshift(Y);
        cvm::rvector freq_scale = fid.GetFreqScale();
        QPoint pos = QCursor::pos();

        double pivot = m_session->m_show_flags.pivot;

        double dphi0 = -( m_start_pos_y - pos.y())*pivot/1000000.0;
        double dphi1 = ( m_start_pos_y - pos.y())/1000000.0;

        for (int n = 1; n < Y.size() + 1; n++ )
        {
            //Y(n) = Y(n) * exp(tcomplex(0, ( m_start_pos - pos.y())*(-pivot+freq_scale(n))/100000.0));
            
            Y(n) = Y(n) * exp(tcomplex(0, ( dphi0 + dphi1*freq_scale(n))));
        }

        Y = fftshift(Y);
        y = ifft(Y);

        fid.SetPhi0(m_session->m_voxel, fid.GetPhi0(m_session->m_voxel)+dphi0);
        fid.SetPhi1(m_session->m_voxel, fid.GetPhi1(m_session->m_voxel)+dphi1/(2.0*M_PI));

        m_start_pos_y = pos.y();

	    m_session->Update(false);
    }
    
    return false;
}

void MainWindow::mouseMoveEvent()
{

    QPoint pos = QCursor::pos();
    std::cout << pos.x() << ", " << pos.y() << std::endl;
    //m_cursor.move( e.pos() );
}


/*
void MainWindow::wheelEvent(QWheelEvent* event)
{
    m_view->setTransformationAnchor(QGraphicsView::AnchorUnderMouse);
 
    // Scale the view / do the zoom
    double scaleFactor = 1.15;
    if(event->delta() > 0) {
        // Zoom in
        m_view->scale(scaleFactor, scaleFactor);
    } else {
        // Zooming out
        m_view->scale(1.0 / scaleFactor, 1.0 / scaleFactor);
    }
}
*/

void MainWindow::OnSelNone()
{
    tarquin::Workspace& workspace = m_session->GetWorkspace();
	tarquin::Options& options = workspace.GetOptions();
	std::vector<tarquin::coord>& fit_list = options.GetFitList();
    fit_list.clear();

    UpdateGeom();
    OnVoxelChange();
}

void MainWindow::OnSelAll()
{
    tarquin::Workspace& workspace = m_session->GetWorkspace();
	tarquin::Options& options = workspace.GetOptions();
	std::vector<tarquin::coord>& fit_list = options.GetFitList();

    // clear old list
    fit_list.clear();

    // make a new one including all voxels
    int rows = m_session->GetWorkspace().GetFID().GetRows();
    int cols = m_session->GetWorkspace().GetFID().GetCols();
    int slices = m_session->GetWorkspace().GetFID().GetSlices();

    
    for ( int slice = 1; slice < slices + 1; slice++ )
        {
            for ( int col = 1; col < cols + 1; col++ )
            {
                for ( int row = 1; row < rows + 1; row++ )
                {
                    coord fit_spec(row, col, slice); 
                    fit_list.push_back(fit_spec);
                }
            }
        }


    m_ui.check_box_fit->setChecked(true);

    UpdateGeom();
    OnVoxelChange();
}

void MainWindow::OnCloseMRI()
{
    m_session->image = MRI();
    m_session->mri_loaded = false;

    disconnect(m_ui.slice_slider, SIGNAL(valueChanged(int)), this, SLOT(OnSliderChange()));
    disconnect(m_ui.slice_spin, SIGNAL(valueChanged(int)), this, SLOT(OnSpinChange()));

    UpdateGeom();
    OnVoxelChange();

    m_ui.slice_slider->setMinimum(0);
    m_ui.slice_slider->setMaximum(0);
    m_ui.slice_slider->setValue(0);

    m_ui.slice_spin->setMinimum(0);
    m_ui.slice_spin->setMaximum(0);
    m_ui.slice_spin->setValue(0);

}

void MainWindow::OnHideGrid()
{
    UpdateGeom();
    OnVoxelChange(false);
}

void MainWindow::OnHideLines()
{
    if ( m_ui.hide_lines->isChecked() )
    {
        m_view->std_vox_pen = Qt::NoPen;
        m_view->in_fit_list_pen = Qt::NoPen;
    }
    else
    {
        m_view->std_vox_pen = QPen(Qt::white,0.2);
        m_view->std_vox_pen.setJoinStyle(Qt::MiterJoin);
        m_view->std_vox_pen.setStyle(Qt::DotLine);
        m_view->in_fit_list_pen = QPen(Qt::white,0.3);
        m_view->in_fit_list_pen.setJoinStyle(Qt::MiterJoin);
    }

    UpdateGeom();
    OnVoxelChange(false);
}

void MainWindow::OnHideMRI()
{
    UpdateGeom();
    OnVoxelChange(false);
}

void MainWindow::OnSetTrans()
{
    bool ok;
    int d = QInputDialog::getInteger(this, tr("Set transparency (%)"),
            tr("Percent:"), m_session->grid_trans, 0, 100, 1, &ok);
    if (ok)
        m_session->grid_trans = d;

    UpdateGeom();
    OnVoxelChange(false);
}

void MainWindow::OnWWLL()
{
    m_ww_wl = false;
}

void MainWindow::OnLoadMRI()
{
    QSettings settings;

	QString path = QFileDialog::getOpenFileName(this, 
			tr("Select MRI file."), 
			settings.value("mri_input_dir", QString()).toString(),
			tr("all files (*)"));
    
    // not sure why this is required but otherwise the button
    // gets stuck after loading on OSX
    m_ui.load_mri->setDown(false); 

	// nothing selected, they must have pushed cancel
	if( !path.size() )
		return;
	
    m_session->mri_loaded = false;

	// they have loaded the file, so keep the directory for next time
	settings.setValue("mri_input_dir", path);

    //std::string strFilename = "/Users/martin/Dropbox/home/python/PYQT_dcm_viewer/eg_case/T2_TSE_AX";

    std::string strFilename = path.toStdString();

    //m_session->image.load_from_file(strFilename);
    
    // clear any previous MRI
    m_session->image = MRI();

    if ( !m_session->image.load_from_file(strFilename) )
    {
        ErrorDialog(this, tr("Problem loading MRI."), tr("Problem loading MRI."));
        return;
    }

    // have a look for other matching MRI files
     // search for other files in the same dir that belong to the same series
    boost::filesystem::path p(strFilename);
    boost::filesystem::path dir = p.parent_path();

    bool image_appended = false;
    // check we found a UID
    if ( m_session->image.get_uid().size() > 0 )
    {
        boost::filesystem::directory_iterator iterator(dir);
        for(; iterator != boost::filesystem::directory_iterator(); ++iterator)
        {
            if ( (iterator->path() != p) && m_session->image.match_uid(iterator->path().string()) ) // don't want the same file twice
            {
                MRI temp_mri;
                if  ( temp_mri.load_from_file(iterator->path().string()) )
                {
                    std::cout << "Loading : " << iterator->path().string() << std::endl;
                    m_session->image.append_slices(temp_mri);
                    image_appended = true;
                }
            }
        }
    }

    if ( image_appended )
    {
        // reorder slices
        m_session->image.reorder_slices();

        // generate qt imaages
        m_session->image.generate_slices();
    }

    std::vector<double> pos = m_session->GetWorkspace().GetFID().GetPos();
    Eigen::Vector3d mrs_pos;
    mrs_pos << pos[0], pos[1], pos[2];

    size_t current_slice = m_session->image.get_closest_slice(mrs_pos);
    
    m_ui.slice_slider->setMinimum(1);
    m_ui.slice_slider->setMaximum(m_session->image.get_num_slices());
    m_ui.slice_slider->setValue(current_slice);
    
    m_ui.slice_spin->setMinimum(1);
    m_ui.slice_spin->setMaximum(m_session->image.get_num_slices());
    m_ui.slice_spin->setValue(current_slice);

    connect(m_ui.slice_slider, SIGNAL(valueChanged(int)), this, SLOT(OnSliderChange()));
    connect(m_ui.slice_spin, SIGNAL(valueChanged(int)), this, SLOT(OnSpinChange()));

    //m_ui.slice_slider->show();
    m_session->mri_loaded = true;

    m_session->image.print_paras();


    UpdateGeom();
    OnVoxelChange();
}


void MainWindow::OnFitChange()
{
    int row = m_ui.spinRows->value();
    int col = m_ui.spinCols->value();
    int slice = m_ui.spinSlices->value();
    
    coord current(row, col, slice); 
    
    // find out if current voxel is in the fit list, and if so what index
    tarquin::Workspace& workspace = m_session->GetWorkspace();
	tarquin::Options& options = workspace.GetOptions();
	std::vector<tarquin::coord>& fit_list = options.GetFitList();

    std::vector<tarquin::coord>::iterator it;
    it = std::find(fit_list.begin(), fit_list.end(), current);

    if ( it != fit_list.end() )
    {
        // voxel is in the fit list so better remove it
	    //int ind = it - fit_list.begin();
        fit_list.erase(it);
        m_ui.check_box_fit->setChecked(false);
    }
    else
    {
        // voxel is not in the fit list so better append it
        fit_list.push_back(current);
        m_ui.check_box_fit->setChecked(true);

        // sort fit list
        std::sort(fit_list.begin(), fit_list.end());
    }
}

void MainWindow::OnEditFin()
{
    m_ui.plot_frame->setFocus();
}

void MainWindow::OnVoxelChange(bool plot_update)
{

    // if voxel tracking is enabled
    UpdateGeom();

    // return the old voxel back to default

    if ( m_ui.check_box_fit->isChecked() )
        m_scene->SetVoxelFormat(m_session->GetVoxel(),m_view->in_fit_list_pen,2);
    else
        m_scene->SetVoxelFormat(m_session->GetVoxel(),m_view->std_vox_pen,1);

    int row = m_ui.spinRows->value();
    int col = m_ui.spinCols->value();
    int slice = m_ui.spinSlices->value();

    //std::cout << "R:" << row << ", C:" << col << ", S:" << slice << std::endl;
    coord current(row, col, slice); 
    
    // find out if current voxel is in the fit list, and if so what index
    tarquin::Workspace& workspace = m_session->GetWorkspace();
	tarquin::Options& options = workspace.GetOptions();
	std::vector<tarquin::coord>& fit_list = options.GetFitList();

    std::vector<tarquin::coord>::iterator it;
    it = std::find(fit_list.begin(), fit_list.end(), current);

    //std::cout << "Fit list size : " << fit_list.size() << std::endl;
    //std::cout << "1 m_fit_number : " << m_session->m_fit_number << std::endl;

    if ( it != fit_list.end() )
    {

        m_session->m_in_fit_list = true;
	    m_session->m_fit_number = it - fit_list.begin();
        m_ui.check_box_fit->setChecked(true);
    }
    else
    {
        m_session->m_in_fit_list = false;
        m_ui.check_box_fit->setChecked(false);
	    m_session->m_fit_number = 0;
    }

    //std::cout << "2 m_fit_number : " << m_session->m_fit_number << std::endl;

	m_session->SetVoxel(current);

    m_scene->SetVoxelFormat(m_session->GetVoxel(),m_view->selected_pen,3);
    
    if ( plot_update )
        m_session->Update();

    OnEditFin();
    
}

void MainWindow::OnNextFit()
{
	if( !m_session )
		return;

    tarquin::Workspace& workspace = m_session->GetWorkspace();
	tarquin::Options& options = workspace.GetOptions();
	std::vector<tarquin::coord>& fit_list = options.GetFitList();

    if ( fit_list.size() == 0 )
    {
        ErrorDialog(this, tr("Can't advance to next fit"), 
			tr("Fit list empty. Please choose some voxels to be fit first."));
        return;
    }

	m_session->m_fit_number++;
	
	if ( size_t(m_session->m_fit_number) >= fit_list.size() )
		m_session->m_fit_number = 0;

    //std::sort(fit_list.begin(), fit_list.end());

    int row = fit_list[m_session->m_fit_number].row;
    int col = fit_list[m_session->m_fit_number].col;
    int slice = fit_list[m_session->m_fit_number].slice;
    
    /*
    for ( size_t n = 0; n < fit_list.size(); n++ )
    {
        std::cout << "I:" << n << std::endl;
        std::cout << "R:" << fit_list[n].row << std::endl;
        std::cout << "C:" << fit_list[n].col << std::endl;
        //std::cout << "S:" << fit_list[n].slice << std::endl;
        std::cout << std::endl;
    }
    */

    disconnect(m_ui.spinSlices, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));
    disconnect(m_ui.spinRows, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));
    disconnect(m_ui.spinCols, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));

    m_ui.spinRows->setValue(row);
    m_ui.spinCols->setValue(col);
    m_ui.spinSlices->setValue(slice);

    connect(m_ui.spinSlices, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));
    connect(m_ui.spinRows, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));
    connect(m_ui.spinCols, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));
    
    OnVoxelChange();

    //coord current(row, col, slice); 
	//m_session->Update();
}

void MainWindow::OnCenterVoxel()
{
	if( !m_session )
		return;
    
    int start_row = ceil(m_session->GetWorkspace().GetFID().GetRows()/2.0);
    int start_col = ceil(m_session->GetWorkspace().GetFID().GetCols()/2.0);
    int start_slice = ceil(m_session->GetWorkspace().GetFID().GetSlices()/2.0);

    disconnect(m_ui.spinSlices, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));
    disconnect(m_ui.spinRows, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));
    disconnect(m_ui.spinCols, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));

    m_ui.spinRows->setValue(start_row);
    m_ui.spinCols->setValue(start_col);
    m_ui.spinSlices->setValue(start_slice);

    connect(m_ui.spinSlices, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));
    connect(m_ui.spinRows, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));
    connect(m_ui.spinCols, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));

    //coord current(start_row, start_col, start_slice); 
	//m_session->SetVoxel(current);

    OnVoxelChange();
}

void MainWindow::OnPreviousFit()
{
	if( !m_session )
		return;

    tarquin::Workspace& workspace = m_session->GetWorkspace();
	const tarquin::Options& options = workspace.GetOptions();
	const std::vector<tarquin::coord>& fit_list = options.GetFitList();

    if ( fit_list.size() == 0 )
    {
        ErrorDialog(this, tr("Can't advance to next fit"), 
			tr("Fit list empty. Please choose some voxels to be fit first."));
        return;
    }

	m_session->m_fit_number--;

	if ( m_session->m_fit_number <= -1 )
	{
		m_session->m_fit_number = fit_list.size()-1;
	}
    
    int row = fit_list[m_session->m_fit_number].row;
    int col = fit_list[m_session->m_fit_number].col;
    int slice = fit_list[m_session->m_fit_number].slice;

    disconnect(m_ui.spinSlices, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));
    disconnect(m_ui.spinRows, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));
    disconnect(m_ui.spinCols, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));

    m_ui.spinRows->setValue(row);
    m_ui.spinCols->setValue(col);
    m_ui.spinSlices->setValue(slice);

    connect(m_ui.spinSlices, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));
    connect(m_ui.spinRows, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));
    connect(m_ui.spinCols, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));

    OnVoxelChange();

    //coord current(row, col, slice); 
	//m_session->Update();
}


void MainWindow::OnViewImag()
{
	if( !m_session )
		return;

	m_session->m_show_flags.mode = Session::IMAG;
	m_session->Update();
}

void MainWindow::OnViewAbs()
{
	if( !m_session )
		return;

	m_session->m_show_flags.mode = Session::ABS;
	m_session->Update();
}


void MainWindow::OnViewTimeDomainSec()
{
	if( !m_session )
		return;

	m_session->m_show_flags.domain   = Session::TIME_DOMAIN;
	m_session->m_show_flags.units_td = Session::SEC;
    m_plot->setAxisAutoScale(QwtPlot::xBottom);
    m_plot->setAxisAutoScale(QwtPlot::yLeft);
	m_session->Update();
}

void MainWindow::OnViewTimeDomainPts()
{
	if( !m_session )
		return;

	m_session->m_show_flags.domain   = Session::TIME_DOMAIN;
	m_session->m_show_flags.units_td = Session::TD_PTS;
    m_plot->setAxisAutoScale(QwtPlot::xBottom);
    m_plot->setAxisAutoScale(QwtPlot::yLeft);
	m_session->Update();
}


void MainWindow::OnViewFreqDomainPPM()
{
	if( !m_session )
		return;

	m_session->m_show_flags.domain   = Session::FREQUENCY_DOMAIN;
	m_session->m_show_flags.units_fd = Session::PPM;
    m_plot->setAxisAutoScale(QwtPlot::xBottom);
    m_plot->setAxisAutoScale(QwtPlot::yLeft);
	m_session->Update();
}

void MainWindow::OnViewFreqDomainHz()
{
	if( !m_session )
		return;

	m_session->m_show_flags.domain   = Session::FREQUENCY_DOMAIN;
	m_session->m_show_flags.units_fd = Session::HERTZ;
    m_plot->setAxisAutoScale(QwtPlot::xBottom);
    m_plot->setAxisAutoScale(QwtPlot::yLeft);
	m_session->Update();
}

void MainWindow::OnViewFreqDomainPts()
{
	if( !m_session )
		return;

	m_session->m_show_flags.domain   = Session::FREQUENCY_DOMAIN;
	m_session->m_show_flags.units_fd = Session::FD_PTS;
    m_plot->setAxisAutoScale(QwtPlot::xBottom);
    m_plot->setAxisAutoScale(QwtPlot::yLeft);
	m_session->Update();
}


void MainWindow::OnSetLB()
{
	if( !m_session )
		return;

    tarquin::Workspace& workspace = m_session->GetWorkspace();
	tarquin::Options& options = workspace.GetOptions();

     bool ok;
     double d = QInputDialog::getDouble(this, tr("Set line broadening"), tr("Line broadining (Hz):"), options.Getlb(), 0, 100000, 1, &ok);

     if (ok)
         options.Setlb(d);
   
	m_session->Update();
}

void MainWindow::OnSetZF()
{
	if( !m_session )
		return;
    
    tarquin::Workspace& workspace = m_session->GetWorkspace();
	tarquin::Options& options = workspace.GetOptions();

     bool ok;
     int d = QInputDialog::getInteger(this, tr("Set zero-filling factor"),
                                        tr("Factor:"), options.GetZF(), 1, 1000000, 1, &ok);
     if (ok)
         options.SetZF(d);
   
	m_session->Update();
}

void MainWindow::OnSetBS()
{
	if( !m_session )
		return;
    
    tarquin::Workspace& workspace = m_session->GetWorkspace();
	tarquin::Options& options = workspace.GetOptions();

     bool ok;
     double d = QInputDialog::getInteger(this, tr("Set number of data points for smoothing"),
                                        tr("Data points:"), options.GetBL(), 1, 100000, 2, &ok);
     if (ok)
         options.SetBL(d);
   
	m_session->Update();
}

void MainWindow::OnSubtractBaseline()
{
	if( !m_session )
		return;

	m_session->m_show_flags.SB   = !m_session->m_show_flags.SB;

	m_session->Update();
}

void MainWindow::OnShowRawData()
{
	if( !m_session )
		return;

	m_session->m_show_flags.show_input   = !m_session->m_show_flags.show_input;

	m_session->Update();
}

void MainWindow::OnAAToggle()
{
	if( !m_session )
		return;

	m_session->m_show_flags.AA   = !m_session->m_show_flags.AA;

	m_session->Update();
}

void MainWindow::OnWUSToggle()
{
	if( !m_session )
		return;

	m_session->m_show_flags.WUS = !m_session->m_show_flags.WUS;

	m_session->Update();
}

void MainWindow::OnSetLeftppm()
{
    if( !m_session )
        return;

    tarquin::Workspace& workspace = m_session->GetWorkspace();
    tarquin::Options& options = workspace.GetOptions();

    bool ok;
    double d = QInputDialog::getDouble(this, tr("Set the left limit"),
            tr("Limit (ppm):"), options.GetPPMend(), options.GetPPMstart(), 100000, 2, &ok);
    if (ok)
        options.SetPPMend(d);

    m_session->Update();
}

void MainWindow::OnSetRightppm()
{
    if( !m_session )
        return;

    tarquin::Workspace& workspace = m_session->GetWorkspace();
    tarquin::Options& options = workspace.GetOptions();

    bool ok;
    double d = QInputDialog::getDouble(this, tr("Set the right limit"),
            tr("Limit (ppm):"), options.GetPPMstart(), -100000, options.GetPPMend(), 2, &ok);
    if (ok)
        options.SetPPMstart(d);

    m_session->Update();
}


MyGraphicsScene::MyGraphicsScene ( QWidget * parent )
: QGraphicsScene(parent) 
{}


MyGraphicsView::MyGraphicsView ( QWidget * parent )
: QGraphicsView(parent) 
{
    selected_pen = QPen(Qt::red,0.6);
    selected_pen.setStyle(Qt::DotLine);
    selected_pen.setJoinStyle(Qt::MiterJoin);
    selected_brush = QBrush(Qt::black);

    std_vox_pen = QPen(Qt::white,0.2);
    //std_vox_pen = Qt::NoPen;
    std_vox_pen.setJoinStyle(Qt::MiterJoin);
    std_vox_pen.setStyle(Qt::DotLine);
    std_vox_brush = QBrush(Qt::black);

    in_fit_list_pen = QPen(Qt::white,0.3);
    //in_fit_list_pen = Qt::NoPen;
    in_fit_list_pen.setJoinStyle(Qt::MiterJoin);
    in_fit_list_brush = QBrush(Qt::yellow);
}

MyGraphicsView::MyGraphicsView ( QGraphicsScene * scene, QWidget * parent)
: QGraphicsView(scene,parent) {}

void MyGraphicsScene::SetVoxelFormat(const tarquin::coord& voxel, const QPen& pen, int Z)
{
    QList <QGraphicsItem *> list;
    list = this->items();

    QList<QGraphicsItem *>::iterator i;

    for (i = list.begin(); i != list.end(); ++i)
    {
        int row = (*i)->data(0).toInt();
        int col = (*i)->data(1).toInt();
        int slice = (*i)->data(2).toInt();
        if ( row == voxel.row && col == voxel.col && slice == voxel.slice )
        {
            QGraphicsPolygonItem *r = dynamic_cast<QGraphicsPolygonItem*>(*i);
            if (r != NULL)
            {
                r->setPen(pen);
                r->setZValue(Z);
            }
            return;
        }
    }
}

void MainWindow::mousePressEvent ( QMouseEvent * event )
{
    if ( !m_session )
        return;

    Qt::MouseButtons mouseButtons = event->buttons();
    Qt::KeyboardModifiers key_mod  = event->modifiers();
    if( mouseButtons == Qt::LeftButton )
    {
        if ( (m_session->mri_loaded) )
        {
            //std::cout << m_view->mapToScene(m_view->mapFrom(this,event->pos())).x() << " " << m_view->mapToScene(m_view->mapFrom(this,event->pos())).y() << std::endl;
            //std::cout << event->pos().x() << " " << event->pos().y() << std::endl;
        }
        
        if ( !m_ui.ww_wl->isChecked() && m_ui.grid_frame->underMouse() ) // also check it is within the MRSI window
        {
            QGraphicsPolygonItem *r = dynamic_cast<QGraphicsPolygonItem*>(m_view->itemAt(m_view->mapFrom(this,event->pos())));
            if (r != NULL)
            {
                /*
                   r->setPen(m_view->selected_pen);
                   r->setZValue(2);
                 */
                int row = r->data(0).toInt();
                int col = r->data(1).toInt();
                int slice = r->data(2).toInt();

                //coord current(row, col, slice); 
                //m_session->SetVoxel(current);

                disconnect(m_ui.spinSlices, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));
                disconnect(m_ui.spinRows, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));
                disconnect(m_ui.spinCols, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));

                m_ui.spinRows->setValue(row);
                m_ui.spinCols->setValue(col);
                m_ui.spinSlices->setValue(slice);

                connect(m_ui.spinSlices, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));
                connect(m_ui.spinRows, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));
                connect(m_ui.spinCols, SIGNAL(valueChanged(int)), this, SLOT(OnVoxelChange()));

                OnVoxelChange();
                //OnEditFin();
            }
        }
        else if (m_session->mri_loaded) 
        {
                //std::cout << "Installing filter" << std::endl;
                qApp->installEventFilter(this);
                // record where we are in pixels
                QPoint pos = QCursor::pos();
                m_start_pos_y = pos.y();
                m_start_pos_x = pos.x();
                m_ww_wl = true;

            /*
            if ( !m_ww_wl ) // going into new mode
            {
                //std::cout << "Installing filter" << std::endl;
                qApp->installEventFilter(this);
                // record where we are in pixels
                QPoint pos = QCursor::pos();
                m_start_pos_y = pos.y();
                m_start_pos_x = pos.x();
                m_ww_wl = true;
            }
            else
            {
                //std::cout << "Removing filter" << std::endl;
                qApp->removeEventFilter(this);

                m_session->image.generate_slices();
                m_ww_wl = false;
            }*/
        }
    }

    else if( mouseButtons == Qt::RightButton && !m_session->data_preprocessed && key_mod == Qt::ShiftModifier )
    {
        QGraphicsPolygonItem *r = dynamic_cast<QGraphicsPolygonItem*>(m_view->itemAt(m_view->mapFrom(this,event->pos())));
        if (r != NULL)
        {
            coord current = m_session->GetVoxel();
            
            int sel_row = r->data(0).toInt();
            int sel_col = r->data(1).toInt();
            int sel_slice = r->data(2).toInt();
            coord selected(sel_row, sel_col, sel_slice);
            
            int max_row, min_row;
            if ( selected.row > current.row )
            {
                max_row = selected.row;
                min_row = current.row;
            }
            else
            {
                max_row = current.row;
                min_row = selected.row;
            }
            int max_col, min_col;
            if ( selected.col > current.col )
            {
                max_col = selected.col;
                min_col = current.col;
            }
            else
            {
                max_col = current.col;
                min_col = selected.col;
            }
            int max_slice, min_slice;
            if ( selected.slice > current.slice )
            {
                max_slice = selected.slice;
                min_slice = current.slice;
            }
            else
            {
                max_slice = current.slice;
                min_slice = selected.slice;
            }
            
            tarquin::Workspace& workspace = m_session->GetWorkspace();
            tarquin::Options& options = workspace.GetOptions();
            std::vector<tarquin::coord>& fit_list = options.GetFitList();

            for ( int slice = min_slice; slice < max_slice + 1; slice++ )
            {
                for ( int col = min_col; col < max_col + 1; col++ )
                {
                    for ( int row = min_row; row < max_row +1; row++ )
                    {
                        coord temp(row, col, slice);
                        std::vector<tarquin::coord>::iterator it;
                        it = std::find(fit_list.begin(), fit_list.end(), temp);
                        if ( it == fit_list.end() )
                            fit_list.push_back(temp);
                    }
                }
            }
            m_ui.check_box_fit->setChecked(true);

            std::sort(fit_list.begin(), fit_list.end());
            UpdateGeom();
            OnVoxelChange();
        }

    }

    else if( mouseButtons == Qt::RightButton && !m_session->data_preprocessed )
    {
        QGraphicsPolygonItem *r = dynamic_cast<QGraphicsPolygonItem*>(m_view->itemAt(m_view->mapFrom(this,event->pos())));
        if (r != NULL)
        {
            /*
            r->setPen(m_view->selected_pen);
            r->setZValue(2);
            */
            int row = r->data(0).toInt();
            int col = r->data(1).toInt();
            int slice = r->data(2).toInt();

            coord selected(row, col, slice); 

            tarquin::Workspace& workspace = m_session->GetWorkspace();
            tarquin::Options& options = workspace.GetOptions();
            std::vector<tarquin::coord>& fit_list = options.GetFitList();

            std::vector<tarquin::coord>::iterator it;
            it = std::find(fit_list.begin(), fit_list.end(), selected);

            if ( it != fit_list.end() )
            {
                // voxel is in the fit list so better remove it
                //int ind = it - fit_list.begin();
                fit_list.erase(it);

                if ( selected == m_session->GetVoxel() )
                {
                    m_ui.check_box_fit->setChecked(false);
                }
                else
                {
                    r->setPen(m_view->std_vox_pen);
                    r->setZValue(1);
                }

            }
            else
            {
                // voxel is not in the fit list so better append it
                fit_list.push_back(selected);

                // sort fit list
                std::sort(fit_list.begin(), fit_list.end());

                if ( selected == m_session->GetVoxel() )
                {
                    m_ui.check_box_fit->setChecked(true);
                }
                else
                {
                    r->setPen(m_view->in_fit_list_pen);
                    r->setZValue(2);
                }
            }
        }
    }
}


double MainWindow::interpolate( double val, double y0, double x0, double y1, double x1 ) {
    return (val-x0)*(y1-y0)/(x1-x0) + y0;
}

double MainWindow::base( double val ) {
    val = 2*(val - 0.5);
    if ( val <= -0.75 ) return 0;
    else if ( val <= -0.25 ) return interpolate( val, 0.0, -0.75, 1.0, -0.25 );
    else if ( val <= 0.25 ) return 1.0;
    else if ( val <= 0.75 ) return interpolate( val, 1.0, 0.25, 0.0, 0.75 );
    else return 0.0;
}

double MainWindow::red( double gray ) {
    return base( gray - 0.5 );
}
double MainWindow::green( double gray ) {
    return base( gray );
}
double MainWindow::blue( double gray ) {
    return base( gray + 0.5 );
}

void MainWindow::UpdateGeom()
{

    // Clear previous scene
    qDeleteAll(m_scene->items());

    // Draw geom data

    treal mrs_row_dim;
    treal mrs_col_dim;

    if ( m_session->GetWorkspace().GetFID().IsKnownVoxelDim() )
    {
        mrs_row_dim = m_session->GetWorkspace().GetFID().GetVoxelDim()[0];
        mrs_col_dim = m_session->GetWorkspace().GetFID().GetVoxelDim()[1];
    }
    else
    {
        mrs_row_dim = 10;
        mrs_col_dim = 10;
    }

    int mrs_rows = m_session->GetWorkspace().GetFID().GetRows();
    int mrs_cols = m_session->GetWorkspace().GetFID().GetCols();
    int mrs_slices = m_session->GetWorkspace().GetFID().GetSlices();

    if ( m_session->mri_loaded )
    {
        m_scene->setSceneRect(0,0,m_session->image.get_col_fov(),m_session->image.get_row_fov());
        QGraphicsRectItem *r = new QGraphicsRectItem(QRectF(0,0,m_session->image.get_col_fov(),m_session->image.get_row_fov()));
        //r->setBrush(QBrush(Qt::black));
        r->setBrush(QBrush(Qt::NoBrush));
        QPen box = QPen(Qt::white,1.0);
        box.setJoinStyle(Qt::MiterJoin);
        r->setPen(box);
        r->setZValue(0);
        m_scene->addItem(r);
    }
    else
        m_scene->setSceneRect(0,0,mrs_cols*mrs_col_dim,mrs_rows*mrs_row_dim);

    /*
    QGraphicsRectItem *r = new QGraphicsRectItem(QRectF(0,0, cols*col_dim, rows*row_dim));
    //r->setBrush(QBrush(Qt::black));
    r->setBrush(QBrush(Qt::NoBrush));
    QPen box = QPen(Qt::white,2.5);
    box.setJoinStyle(Qt::MiterJoin);
    r->setPen(box);
    m_scene->addItem(r);
    */


	tarquin::Workspace& workspace = m_session->GetWorkspace();
	tarquin::Options& options = workspace.GetOptions();
    std::vector<tarquin::coord>& fit_list = options.GetFitList();
    
    cvm::rvector map;
    
    if ( m_session->data_fitted )
    {
        //m_session->map;

        map.resize(m_session->map.size());
        map = m_session->map; 
        /*
        const CBasis& basis     = workspace.GetBasis();
	    const rvec_stdvec ahat  = workspace.GetAmplitudesNormalised();
	    std::vector < std::string > signal_names = basis.GetSignalNames();	
        map.resize(ahat.size());
        for ( int n = 0; n < ahat.size(); n++ )
            map(n+1) = ahat[n](22) + ahat[n](23); // get a metabolite
            */
    }
    else
    {
	    if( options.GetFilenameWater() == "" )
        {
            cvm::rvector temp_map = m_session->GetWorkspace().GetFIDRaw().GetRawMap();
            map.resize(temp_map.size());
            map = temp_map;
        }
        else
        {
            tarquin::CFID water = m_session->GetWorkspace().GetFIDWater();
            
            map.resize(water.GetVoxelCount());

            int rows = water.GetRows();
            int cols = water.GetCols();
            int slices = water.GetSlices();
            
            int n = 1;
            for ( int slice = 0; slice < slices; slice++ )
            {
                for ( int col = 0; col < cols; col++ )
                {
                    for ( int row = 0; row < rows; row++ )
                    {
                        coord current(row+1, col+1, slice+1); 
                        map(n) = ComputeWaterNormalisation(current, water, GetLog());
                        n++;
                    }
                }
            }


            //cvm::rvector temp_map = m_session->GetWorkspace().GetFIDWater().GetRawMap();
        }
    }

    double min = map(map.indofmin());
    cvm::rvector min_vec(map.size());
    min_vec.set(min);
    cvm::rvector img_map = map - min_vec;
    double max = img_map(img_map.indofmax());

    if ( max == 0 )
        max = 1;

    img_map = img_map/max;

    for ( int slice = 0; slice < mrs_slices; slice++ )
    {
        for ( int col = 0; col < mrs_cols; col++ )
        {
            for ( int row = 0; row < mrs_rows; row++ )
            {
                if ( m_session->mri_loaded && m_session->GetWorkspace().GetFID().IsKnownPos() )
                {
                    size_t this_slice = m_ui.slice_slider->value();
                    // need to check voxel intersects
                    double di = m_session->image.get_col_dim();
                    double dj = m_session->image.get_row_dim();
                    double dk = m_session->image.get_slice_dim();

                    Eigen::Vector3d X = m_session->image.get_row_dirn(this_slice);
                    X.normalize();
                    Eigen::Vector3d Y = m_session->image.get_col_dirn(this_slice);
                    Y.normalize();

                    Eigen::Vector3d IPP = m_session->image.get_pos(this_slice);
                    Eigen::Vector3d mri_row_dir = X; 
                    Eigen::Vector3d mri_col_dir = Y; 
                    Eigen::Vector3d mri_slice_dir = mri_row_dir.cross(mri_col_dir);
                    mri_slice_dir.normalize();

                    // slice correction
                    //IPP = IPP + di*mri_row_dir - dj*mri_col_dir;  // not sure why this is needed but does seem to 
                                                                  // help geom agreement with Philips software
                    Eigen::Vector3d mri_pos = IPP; 


                    Eigen::Matrix4d M;
                    /*M << X(0)*di, Y(0)*dj, 0, IPP(0),
                             X(1)*di, Y(1)*dj, 0, IPP(1), 
                             X(2)*di, Y(2)*dj, 0, IPP(2),
                             0, 0, 0, 1;*/

                    M << X(0), Y(0), 0, IPP(0),
                             X(1), Y(1), 0, IPP(1), 
                             X(2), Y(2), 0, IPP(2),
                             0, 0, 0, 1;

                    Eigen::Matrix3d M_cut;
                    // dj and di is accounted for by using proper scene coords and rescaling MRI
                    //M_cut << X(0)*di, Y(0)*dj, IPP(0), X(1)*di, Y(1)*dj, IPP(1), X(2)*di, Y(2)*dj, IPP(2);
                    //M_cut << X(0), Y(0), IPP(0), X(1), Y(1), IPP(1), X(2), Y(2), IPP(2);
                    //std::cout << std::endl << M_cut << std::endl;

                    std::vector<double> rdir = m_session->GetWorkspace().GetFID().GetRowDirn();
                    Eigen::Vector3d mrs_row_dir;
                    mrs_row_dir << rdir[0], rdir[1], rdir[2];
                    mrs_row_dir.normalize();

                    std::vector<double> cdir = m_session->GetWorkspace().GetFID().GetColDirn();
                    Eigen::Vector3d mrs_col_dir;
                    mrs_col_dir << cdir[0], cdir[1], cdir[2];
                    mrs_col_dir.normalize();

                    Eigen::Vector3d mrs_slice_dir = mrs_row_dir.cross(mrs_col_dir);
                    mrs_slice_dir.normalize();

                    std::vector<double> pos = m_session->GetWorkspace().GetFID().GetPos();
                    Eigen::Vector3d mrs_pos;
                    mrs_pos << pos[0], pos[1], pos[2];

                    std::vector<double> vox_dim = m_session->GetWorkspace().GetFID().GetVoxelDim();
                    double mrs_row_dim = vox_dim[0];
                    double mrs_col_dim = vox_dim[1];
                    double mrs_slice_dim = vox_dim[2];

                    Eigen::Vector3d slice_norm = mri_slice_dir;
                    Eigen::Vector3d slice_point = mri_pos;
                    
                    Eigen::Vector3d offset;
                    offset << 0, 0, 0;

                    offset =  col*mrs_row_dim*mrs_row_dir + row*mrs_col_dim*mrs_col_dir + slice*mrs_slice_dim*mrs_slice_dir;

                    std::vector<Point> vox_points;

                    // SLICE wise
                    // point intersecting with cube edge
                    Eigen::Vector3d line_norm = mrs_slice_dir; 
                    {
                        Eigen::Vector3d line_point = mrs_pos + offset + mrs_row_dir*mrs_row_dim/2.0 + mrs_col_dir*mrs_col_dim/2.0 - mrs_slice_dir*mrs_slice_dim/2.0;
                        // add intersecting line
                        double numer = (slice_point - line_point).dot(slice_norm);
                        double denom = line_norm.dot(slice_norm);
                        if ((numer != 0) && (denom != 0))
                        {
                            double d = numer / denom;
                            if ((d > 0) && (d < mrs_slice_dim))
                            {
                                Point vox_pt = calc_pt(d, line_norm, line_point, M);
                                vox_points.push_back(vox_pt);
                            }
                        }
                    }
                    {
                        Eigen::Vector3d line_point = mrs_pos + offset + mrs_row_dir*mrs_row_dim/2.0 - mrs_col_dir*mrs_col_dim/2.0 - mrs_slice_dir*mrs_slice_dim/2.0;
                        // add intersecting line
                        double numer = (slice_point - line_point).dot(slice_norm);
                        double denom = line_norm.dot(slice_norm);
                        if ((numer != 0) && (denom != 0))
                        {
                            double d = numer / denom;
                            if ((d > 0) && (d < mrs_slice_dim))
                            {
                                Point vox_pt = calc_pt(d, line_norm, line_point, M);
                                vox_points.push_back(vox_pt);
                            }
                        }
                    }
                    {
                        Eigen::Vector3d line_point = mrs_pos + offset - mrs_row_dir*mrs_row_dim/2.0 + mrs_col_dir*mrs_col_dim/2.0 - mrs_slice_dir*mrs_slice_dim/2.0;
                        // add intersecting line
                        double numer = (slice_point - line_point).dot(slice_norm);
                        double denom = line_norm.dot(slice_norm);
                        if ((numer != 0) && (denom != 0))
                        {
                            double d = numer / denom;
                            if ((d > 0) && (d < mrs_slice_dim))
                            {
                                Point vox_pt = calc_pt(d, line_norm, line_point, M);
                                vox_points.push_back(vox_pt);
                            }
                        }
                    }
                    {
                        Eigen::Vector3d line_point = mrs_pos + offset - mrs_row_dir*mrs_row_dim/2.0 - mrs_col_dir*mrs_col_dim/2.0 - mrs_slice_dir*mrs_slice_dim/2.0;
                        // add intersecting line
                        double numer = (slice_point - line_point).dot(slice_norm);
                        double denom = line_norm.dot(slice_norm);
                        if ((numer != 0) && (denom != 0))
                        {
                            double d = numer / denom;
                            if ((d > 0) && (d < mrs_slice_dim))
                            {
                                Point vox_pt = calc_pt(d, line_norm, line_point, M);
                                vox_points.push_back(vox_pt);
                            }
                        }
                    }

                    // ROW wise
                    // point intersecting with cube edge
                    line_norm = mrs_row_dir; 
                    {
                        Eigen::Vector3d line_point = mrs_pos + offset - mrs_row_dir*mrs_row_dim/2.0 + mrs_col_dir*mrs_col_dim/2.0 + mrs_slice_dir*mrs_slice_dim/2.0;
                        // add intersecting line
                        double numer = (slice_point - line_point).dot(slice_norm);
                        double denom = line_norm.dot(slice_norm);
                        if ((numer != 0) && (denom != 0))
                        {
                            double d = numer / denom;
                            if ((d > 0) && (d < mrs_row_dim))
                            {
                                Point vox_pt = calc_pt(d, line_norm, line_point, M);
                                vox_points.push_back(vox_pt);
                            }
                        }
                    }
                    {
                        Eigen::Vector3d line_point = mrs_pos + offset - mrs_row_dir*mrs_row_dim/2.0 + mrs_col_dir*mrs_col_dim/2.0 - mrs_slice_dir*mrs_slice_dim/2.0;
                        // add intersecting line
                        double numer = (slice_point - line_point).dot(slice_norm);
                        double denom = line_norm.dot(slice_norm);
                        if ((numer != 0) && (denom != 0))
                        {
                            double d = numer / denom;
                            if ((d > 0) && (d < mrs_row_dim))
                            {
                                Point vox_pt = calc_pt(d, line_norm, line_point, M);
                                vox_points.push_back(vox_pt);
                            }
                        }
                    }
                    {
                        Eigen::Vector3d line_point = mrs_pos + offset - mrs_row_dir*mrs_row_dim/2.0 - mrs_col_dir*mrs_col_dim/2.0 + mrs_slice_dir*mrs_slice_dim/2.0;
                        // add intersecting line
                        double numer = (slice_point - line_point).dot(slice_norm);
                        double denom = line_norm.dot(slice_norm);
                        if ((numer != 0) && (denom != 0))
                        {
                            double d = numer / denom;
                            if ((d > 0) && (d < mrs_row_dim))
                            {
                                Point vox_pt = calc_pt(d, line_norm, line_point, M);
                                vox_points.push_back(vox_pt);
                            }
                        }
                    }
                    {
                        Eigen::Vector3d line_point = mrs_pos + offset - mrs_row_dir*mrs_row_dim/2.0 - mrs_col_dir*mrs_col_dim/2.0 - mrs_slice_dir*mrs_slice_dim/2.0;
                        // add intersecting line
                        double numer = (slice_point - line_point).dot(slice_norm);
                        double denom = line_norm.dot(slice_norm);
                        if ((numer != 0) && (denom != 0))
                        {
                            double d = numer / denom;
                            if ((d > 0) && (d < mrs_row_dim))
                            {
                                Point vox_pt = calc_pt(d, line_norm, line_point, M);
                                vox_points.push_back(vox_pt);
                            }
                        }
                    }
                    
                    // COL wise
                    // point intersecting with cube edge
                    line_norm = mrs_col_dir; 
                    {
                        Eigen::Vector3d line_point = mrs_pos + offset + mrs_row_dir*mrs_row_dim/2.0 - mrs_col_dir*mrs_col_dim/2.0 + mrs_slice_dir*mrs_slice_dim/2.0;
                        // add intersecting line
                        double numer = (slice_point - line_point).dot(slice_norm);
                        double denom = line_norm.dot(slice_norm);
                        if ((numer != 0) && (denom != 0))
                        {
                            double d = numer / denom;
                            if ((d > 0) && (d < mrs_col_dim))
                            {
                                Point vox_pt = calc_pt(d, line_norm, line_point, M);
                                vox_points.push_back(vox_pt);
                            }
                        }
                    }
                    {
                        Eigen::Vector3d line_point = mrs_pos + offset + mrs_row_dir*mrs_row_dim/2.0 - mrs_col_dir*mrs_col_dim/2.0 - mrs_slice_dir*mrs_slice_dim/2.0;
                        // add intersecting line
                        double numer = (slice_point - line_point).dot(slice_norm);
                        double denom = line_norm.dot(slice_norm);
                        if ((numer != 0) && (denom != 0))
                        {
                            double d = numer / denom;
                            if ((d > 0) && (d < mrs_col_dim))
                            {
                                Point vox_pt = calc_pt(d, line_norm, line_point, M);
                                vox_points.push_back(vox_pt);
                            }
                        }
                    }
                    {
                        Eigen::Vector3d line_point = mrs_pos + offset - mrs_row_dir*mrs_row_dim/2.0 - mrs_col_dir*mrs_col_dim/2.0 + mrs_slice_dir*mrs_slice_dim/2.0;
                        // add intersecting line
                        double numer = (slice_point - line_point).dot(slice_norm);
                        double denom = line_norm.dot(slice_norm);
                        if ((numer != 0) && (denom != 0))
                        {
                            double d = numer / denom;
                            if ((d > 0) && (d < mrs_col_dim))
                            {
                                Point vox_pt = calc_pt(d, line_norm, line_point, M);
                                vox_points.push_back(vox_pt);
                            }
                        }
                    }
                    {
                        Eigen::Vector3d line_point = mrs_pos + offset - mrs_row_dir*mrs_row_dim/2.0 - mrs_col_dir*mrs_col_dim/2.0 - mrs_slice_dir*mrs_slice_dim/2.0;
                        // add intersecting line
                        double numer = (slice_point - line_point).dot(slice_norm);
                        double denom = line_norm.dot(slice_norm);
                        if ((numer != 0) && (denom != 0))
                        {
                            double d = numer / denom;
                            if ((d > 0) && (d < mrs_col_dim))
                            {
                                Point vox_pt = calc_pt(d, line_norm, line_point, M);
                                vox_points.push_back(vox_pt);
                            }
                        }
                    }

                    if ( vox_points.size() > 2 )
                    {
                        std::vector<Point> vox_points_sorted = convex_hull(vox_points);

                        QPolygonF voxel_path;
                        for ( size_t n = 0; n < vox_points.size(); n++ )
                        {
                            //std::cout << vox_points_sorted[n].x << "," << vox_points_sorted[n].y << std::endl;
                            voxel_path << QPointF( vox_points_sorted[n].x, vox_points_sorted[n].y );
                        }
                        
                        QGraphicsPolygonItem *poly = new QGraphicsPolygonItem(voxel_path);

                        int n = 1+row + mrs_rows*col + mrs_rows*mrs_cols*slice;

                        if ( !m_session->data_fitted )
                        {
                            poly->setBrush(QBrush(GetColMap(img_map(n))));
                        }

                        poly->setZValue(1);
                        poly->setData(0,row+1); // row
                        poly->setData(1,col+1); // col
                        poly->setData(2,slice+1); // slice

                        coord current(row+1, col+1, slice+1); 
                        std::vector<tarquin::coord>::iterator it;
                        it = std::find(fit_list.begin(), fit_list.end(), current);

                        if ( it != fit_list.end() )
                        {
                            poly->setPen(m_view->in_fit_list_pen);
                            poly->setZValue(2);

                            if ( m_session->data_fitted )
                            {
                                int ind = it - fit_list.begin() + 1;
                                poly->setBrush(QBrush(GetColMap(img_map(ind))));
                            }
                        }
                        else
                        {
                            poly->setPen(m_view->std_vox_pen);
                            poly->setZValue(1);

                            if ( m_session->data_fitted )
                            {
                                poly->setBrush(QBrush(Qt::NoBrush));
                            }

                        }
                        if ( !m_ui.hide_grid->isChecked() )
                            m_scene->addItem(poly);

                    }



                }
                else
                {

                    if ( m_ui.spinSlices->value() == slice+1 ) // if we're on the correct slice
                    {
                    QGraphicsPolygonItem *rect = new QGraphicsPolygonItem(QRectF(col*mrs_row_dim,row*mrs_col_dim, mrs_row_dim, mrs_col_dim));

                    int n = 1+row + mrs_rows*col + mrs_rows*mrs_cols*slice;

                    //rect->setBrush(QBrush(QColor(255*red(img_map(n)),255*green(img_map(n)),255*blue(img_map(n)),255)));

                    //std::cout << slice <<std::endl;

                    if ( !m_session->data_fitted )
                    {
                        //rect->setBrush(QBrush(GetColMap(img_map(n+mrs_rows*mrs_cols*slice))));
                        rect->setBrush(QBrush(GetColMap(img_map(n))));
                    }

                    rect->setZValue(1);
                    rect->setData(0,row+1); // row
                    rect->setData(1,col+1); // col
                    rect->setData(2,slice+1); // slice

                    coord current(row+1, col+1, slice+1); 
                    std::vector<tarquin::coord>::iterator it;
                    it = std::find(fit_list.begin(), fit_list.end(), current);

                    if ( it != fit_list.end() )
                    {
                        rect->setPen(m_view->in_fit_list_pen);
                        rect->setZValue(2);

                        if ( m_session->data_fitted )
                        {
                            int ind = it - fit_list.begin() + 1;
                            rect->setBrush(QBrush(GetColMap(img_map(ind))));
                        }
                    }
                    else
                    {
                        rect->setPen(m_view->std_vox_pen);
                        rect->setZValue(1);

                        if ( m_session->data_fitted )
                        {
                            //rect->setBrush(QBrush(QColor(Qt::black)));
                            rect->setBrush(QBrush(Qt::NoBrush));
                        }

                    }

                    if ( !m_ui.hide_grid->isChecked() )
                        m_scene->addItem(rect);
                    }
                }
            }
        }
    }

    if ( m_session->mri_loaded )
    {
        QImage image_slice = m_session->image.get_slice(m_ui.slice_slider->value());

        if ( !m_ui.hide_mri->isChecked() )
        {
            QGraphicsItem *slice = m_scene->addPixmap( QPixmap::fromImage(image_slice).scaled(m_session->image.get_col_fov(),m_session->image.get_row_fov()) );
            slice->setZValue(0);
        }

        //slice->setFlag(QGraphicsItem::ItemIsMovable, true);

        m_session->mri_loaded = true;
    }
}

void MainWindow::OnSpinChange()
{
    size_t slice = m_ui.slice_spin->value();  
    m_ui.slice_slider->setValue(slice);
    m_ui.slice_spin->setFocus();
}

QColor MainWindow::GetColMap(double val)
{
    //std::cout << val << std::endl;
    //return QColor(255,255,0,240*val);
    //return QColor(255*red(val),255*green(val),255*blue(val),m_session->grid_trans/100.0*255);
    COLOUR c = GetColour(val,0,1);
    return QColor(255*c.r,255*c.g,255*c.b,255-m_session->grid_trans/100.0*255);
}

void MainWindow::OnSliderChange()
{
    size_t slice = m_ui.slice_slider->value();  

    disconnect(m_ui.slice_spin, SIGNAL(valueChanged(int)), this, SLOT(OnSpinChange()));
    m_ui.slice_spin->setValue(slice);
    connect(m_ui.slice_spin, SIGNAL(valueChanged(int)), this, SLOT(OnSpinChange()));

    /*
    std::cout << slice << std::endl;
    std::cout << "Pos      : " << m_session->image.get_pos(slice) << std::endl;
    std::cout << "Row dirn : " << m_session->image.get_row_dirn(slice) << std::endl;
    std::cout << "Col dirn : " << m_session->image.get_col_dirn(slice) << std::endl;
    */
    
    //m_session->image.print_paras();

    UpdateGeom();
    OnVoxelChange(false);

    m_ui.slice_slider->setFocus();
}



// http://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain#C.2B.2B
// Implementation of Andrew's monotone chain 2D convex hull algorithm.
// Asymptotic complexity: O(n log n).
// Practical performance: 0.5-1.0 seconds for n=1000000 on a 1GHz machine.
// #include <algorithm>
// #include <vector>
// using namespace std;
/* 
typedef double coord_t;         // coordinate type
typedef double coord2_t;  // must be big enough to hold 2*max(|coordinate|)^2
 
struct Point {
        coord_t x, y;
 
        bool operator <(const Point &p) const {
                return x < p.x || (x == p.x && y < p.y);
        }
};
*/
 
// 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
// Returns a positive value, if OAB makes a counter-clockwise turn,
// negative for clockwise turn, and zero if the points are collinear.
coord2_t MainWindow::cross(const Point &O, const Point &A, const Point &B)
{
        return (A.x - O.x) * (B.y - O.y) - (A.y - O.y) * (B.x - O.x);
}
 
// Returns a list of points on the convex hull in counter-clockwise order.
// Note: the last point in the returned list is the same as the first one.
std::vector<Point> MainWindow::convex_hull(std::vector<Point> P)
{
        int n = P.size(), k = 0;
        std::vector<Point> H(2*n);
 
        // Sort points lexicographically
        sort(P.begin(), P.end());
 
        // Build lower hull
        for (int i = 0; i < n; i++) {
                while (k >= 2 && cross(H[k-2], H[k-1], P[i]) <= 0) k--;
                H[k++] = P[i];
        }
 
        // Build upper hull
        for (int i = n-2, t = k+1; i >= 0; i--) {
                while (k >= t && cross(H[k-2], H[k-1], P[i]) <= 0) k--;
                H[k++] = P[i];
        }
 
        H.resize(k);
        return H;
}

Point MainWindow::calc_pt(double d, Eigen::Vector3d line_norm, Eigen::Vector3d line_point, Eigen::Matrix3d M_cut)
{
	Eigen::Vector3d intersect_pt = d*line_norm + line_point;

	/*Eigen::JacobiSVD<Eigen::Matrix3d> svd(M_cut,Eigen::ComputeThinU|Eigen::ComputeThinV);
	Eigen::Matrix3d result;
	svd.pinv(result);
	Eigen::Vector3d pt = result*intersect_pt;
	*/

	Eigen::Vector3d pt = M_cut.inverse()*intersect_pt;
	//std::cout << std::endl << pt << std::endl;
	Point vox_pt;
	vox_pt.x = pt(0);
	vox_pt.y = pt(1);
	return vox_pt;
}

Point MainWindow::calc_pt(double d, Eigen::Vector3d line_norm, Eigen::Vector3d line_point, const Eigen::Matrix4d& M)
{
	Eigen::Vector3d intersect_pt = d*line_norm + line_point;

	Eigen::Vector4d intersect_pt_long;
	intersect_pt_long(0) = intersect_pt(0);
	intersect_pt_long(1) = intersect_pt(1);
	intersect_pt_long(2) = intersect_pt(2);
	intersect_pt_long(3) = 1;

	Eigen::JacobiSVD<Eigen::MatrixXd> svd(M,Eigen::ComputeThinU|Eigen::ComputeThinV);
	//Eigen::JacobiSVD<Eigen::Matrix4d> svd(M);
	Eigen::MatrixXd result;
	svd.pinv(result);
	Eigen::Vector4d pt = result*intersect_pt_long;

	//Eigen::Vector3d pt = M_cut.inverse()*intersect_pt;
	//std::cout << std::endl << pt << std::endl;
	Point vox_pt;
	vox_pt.x = pt(0);
	vox_pt.y = pt(1);
	return vox_pt;
}
