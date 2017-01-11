
#include "proto_buf.hpp"
#include "session.h"
#include "exception.hpp"
#include "inputdlg.h"
#include "guicommon.h"
#include "voxel_select.h"

#include <boost/filesystem.hpp>
#include <QtGui>


InputDlg::InputDlg(QWidget* parent, Session* session) :
	QDialog(parent),
	m_session(session),
	m_loaded_ws(false),
	m_loaded_basis(false),
	m_start_fit(false)
{
	assert( session );

	m_ui.setupUi(this);

	// add the valid options to the combo box
	m_ui.cmbFormat->addItem("DICOM",                 QVariant(tarquin::DCM));
	m_ui.cmbFormat->addItem("Siemens IMA",           QVariant(tarquin::SIEMENS));
	m_ui.cmbFormat->addItem("Siemens RDA",           QVariant(tarquin::RDA));
	m_ui.cmbFormat->addItem("Philips private DICOM", QVariant(tarquin::PHILIPS_DCM));
	m_ui.cmbFormat->addItem("Philips SDAT/SPAR",     QVariant(tarquin::PHILIPS));
	m_ui.cmbFormat->addItem("GE P file",             QVariant(tarquin::GE));
	m_ui.cmbFormat->addItem("GE SHF and P file",     QVariant(tarquin::SHF));
	m_ui.cmbFormat->addItem("Varian FID",            QVariant(tarquin::VARIAN));
	m_ui.cmbFormat->addItem("Bruker FID",            QVariant(tarquin::BRUKER));
	m_ui.cmbFormat->addItem("LCModel RAW",           QVariant(tarquin::LCM));
	m_ui.cmbFormat->addItem("Dangerplot",            QVariant(tarquin::DANGER));
	m_ui.cmbFormat->addItem("JMRUI txt",             QVariant(tarquin::JMRUI_TXT));

    m_ui.cmbRefMode->addItem("1H NAA Cr Cho Lip",    QVariant(tarquin::PROTON_NAA_CR_CHO_LIP));
    m_ui.cmbRefMode->addItem("1H NAA Cr Cho",        QVariant(tarquin::PROTON_NAA_CR_CHO));
    m_ui.cmbRefMode->addItem("1H NAA Cho",           QVariant(tarquin::PROTON_NAA_CHO));
    m_ui.cmbRefMode->addItem("1H Cr Cho",            QVariant(tarquin::PROTON_CR_CHO));
    m_ui.cmbRefMode->addItem("1H NAA",               QVariant(tarquin::PROTON_NAA));
    m_ui.cmbRefMode->addItem("1H Cr",                QVariant(tarquin::PROTON_CR));
    m_ui.cmbRefMode->addItem("1H Cho",               QVariant(tarquin::PROTON_CHO));
    m_ui.cmbRefMode->addItem("1H H2O",               QVariant(tarquin::PROTON_H2O));
    m_ui.cmbRefMode->addItem("1H Lip",               QVariant(tarquin::PROTON_LIP));
    m_ui.cmbRefMode->addItem("31P PCr",              QVariant(tarquin::PHOSPH_PCR));
    m_ui.cmbRefMode->addItem("31P PCr gamma-ATP",    QVariant(tarquin::PHOSPH_PCR_GAMMAATP));

    m_ui.cmbDynRefMode->addItem("1H NAA Cr Cho Lip",    QVariant(tarquin::PROTON_NAA_CR_CHO_LIP));
    m_ui.cmbDynRefMode->addItem("1H NAA Cr Cho",        QVariant(tarquin::PROTON_NAA_CR_CHO));
    m_ui.cmbDynRefMode->addItem("1H NAA Cho",           QVariant(tarquin::PROTON_NAA_CHO));
    m_ui.cmbDynRefMode->addItem("1H Cr Cho",            QVariant(tarquin::PROTON_CR_CHO));
    m_ui.cmbDynRefMode->addItem("1H NAA",               QVariant(tarquin::PROTON_NAA));
    m_ui.cmbDynRefMode->addItem("1H Cr",                QVariant(tarquin::PROTON_CR));
    m_ui.cmbDynRefMode->addItem("1H Cho",               QVariant(tarquin::PROTON_CHO));
    m_ui.cmbDynRefMode->addItem("1H H2O",               QVariant(tarquin::PROTON_H2O));
    m_ui.cmbDynRefMode->addItem("1H Lip",               QVariant(tarquin::PROTON_LIP));
    m_ui.cmbDynRefMode->addItem("31P PCr",              QVariant(tarquin::PHOSPH_PCR));
    m_ui.cmbDynRefMode->addItem("31P PCr gamma-ATP",    QVariant(tarquin::PHOSPH_PCR_GAMMAATP));

    m_ui.cmbIntBasisSet->addItem("1H brain",                QVariant(tarquin::PROTON_BRAIN));
    m_ui.cmbIntBasisSet->addItem("1H brain + Glth",                QVariant(tarquin::PROTON_BRAIN_GLTH));
    m_ui.cmbIntBasisSet->addItem("1H brain + Gly, Glth",                QVariant(tarquin::PROTON_BRAIN_GLY_GLTH));
    m_ui.cmbIntBasisSet->addItem("1H brain + Gly, Cit, Glth",    QVariant(tarquin::PROTON_BRAIN_GLY_CIT_GLTH));
    m_ui.cmbIntBasisSet->addItem("1H brain full",                QVariant(tarquin::PROTON_BRAIN_FULL));
    m_ui.cmbIntBasisSet->addItem("1H brain, long echo",     QVariant(tarquin::PROTON_BRAIN_LE));
    m_ui.cmbIntBasisSet->addItem("1H brain no PCr",                QVariant(tarquin::PROTON_BRAIN_NO_PCR));
    m_ui.cmbIntBasisSet->addItem("1H MEGA-PRESS GABA",      QVariant(tarquin::PROTON_MEGAPRESS_GABA));
    m_ui.cmbIntBasisSet->addItem("1H BRAINO phantom",      QVariant(tarquin::PROTON_BRAINO));
    m_ui.cmbIntBasisSet->addItem("1H brain + Glth + exp MM",      QVariant(tarquin::PROTON_BRAIN_MMEXP));
    m_ui.cmbIntBasisSet->addItem("1H brain + Glth - Lip/MM",      QVariant(tarquin::PROTON_BRAIN_METAB_ONLY));
    m_ui.cmbIntBasisSet->addItem("1H brain LCModel",      QVariant(tarquin::PROTON_BRAIN_LCM));
    m_ui.cmbIntBasisSet->addItem("31P brain, 1H decoupled", QVariant(tarquin::PHOSPH_BRAIN_DECOUP));

    m_ui.cmbDynAv->addItem("Default",                         QVariant(tarquin::DEFAULT));
    m_ui.cmbDynAv->addItem("No averaging",                    QVariant(tarquin::NONE));
    m_ui.cmbDynAv->addItem("Average all scans",               QVariant(tarquin::ALL));
    m_ui.cmbDynAv->addItem("Subtract even from odd",          QVariant(tarquin::SUBTRACT));
    m_ui.cmbDynAv->addItem("Average odd scans only",          QVariant(tarquin::ODD));
    m_ui.cmbDynAv->addItem("Average even scans only",         QVariant(tarquin::EVEN));

    m_ui.cmbDynAv_W->addItem("Default",                         QVariant(tarquin::DEFAULT));
    m_ui.cmbDynAv_W->addItem("No averaging",                    QVariant(tarquin::NONE));
    m_ui.cmbDynAv_W->addItem("Average all scans",               QVariant(tarquin::ALL));
    m_ui.cmbDynAv_W->addItem("Subtract even from odd",          QVariant(tarquin::SUBTRACT));
    m_ui.cmbDynAv_W->addItem("Average odd scans only",          QVariant(tarquin::ODD));
    m_ui.cmbDynAv_W->addItem("Average even scans only",         QVariant(tarquin::EVEN));
    
    m_ui.cmbPS->addItem("PRESS",             QVariant(tarquin::PRESS));
    m_ui.cmbPS->addItem("STEAM",             QVariant(tarquin::STEAM));
    m_ui.cmbPS->addItem("sLASER",             QVariant(tarquin::SEMI_LASER));
    m_ui.cmbPS->addItem("LASER",             QVariant(tarquin::LASER));
    m_ui.cmbPS->addItem("Pulse-acquire",     QVariant(tarquin::PULSE_ACQUIRE));
    m_ui.cmbPS->addItem("CPMG",              QVariant(tarquin::CPMG));
    m_ui.cmbPS->addItem("Spin-echo",         QVariant(tarquin::SPIN_ECHO));
    m_ui.cmbPS->addItem("MEGA-PRESS",        QVariant(tarquin::MEGA_PRESS)); // only used in preprocessor for adjuting phase

	// disable the input controls until a file has been loaded
	m_ui.txtSF->setEnabled(false);
	m_ui.txtTF->setEnabled(false);
	m_ui.txtET->setEnabled(false);
	m_ui.txtTE1->setEnabled(false);
	m_ui.txtTM->setEnabled(false);
	m_ui.txtCN->setEnabled(false);
	m_ui.cmbPS->setEnabled(false);
	m_ui.txtNP->setEnabled(false);
	m_ui.txtRO->setEnabled(false);
	m_ui.txtParaFile->setEnabled(false);
	m_ui.btnOpenParaFile->setEnabled(false);
	m_ui.btnOpenWU->setEnabled(false);
	m_ui.txtWU->setEnabled(false);

	m_ui.txtZF->setEnabled(false);
	m_ui.txtPPMstart->setEnabled(false);
	m_ui.txtPPMend->setEnabled(false);
	m_ui.txtBL->setEnabled(false);
	m_ui.txtLB->setEnabled(false);

	m_ui.btnOpenXML->setEnabled(false);
	m_ui.btnOpenCSV->setEnabled(false);
	m_ui.btnSaveXML->setEnabled(false);
	m_ui.btnSaveLCM->setEnabled(false);
	
    m_ui.txtBasisCSV->setEnabled(false);
    m_ui.txtBasisXML->setEnabled(false);
    m_ui.txtSaveBasis->setEnabled(false);
    m_ui.txtSaveBasisLCM->setEnabled(false);

    m_ui.txtSP->setEnabled(false);
	m_ui.txtEP->setEnabled(false);
	m_ui.txtWC->setEnabled(false);
	m_ui.txtCW->setEnabled(false);
	m_ui.txtLD->setEnabled(false);
	m_ui.txtDR->setEnabled(false);
	m_ui.txtZFKS->setEnabled(false);
	m_ui.txtIM->setEnabled(false);
	m_ui.txtMI->setEnabled(false);
	m_ui.txtIB->setEnabled(false);
	m_ui.txtWConc->setEnabled(false);
	m_ui.txtWA->setEnabled(false);
	m_ui.cbCP->setEnabled(false);
	m_ui.cbLF->setEnabled(false);
	m_ui.cbKPWFS->setEnabled(false);
    m_ui.cbKF->setEnabled(false);
    m_ui.cbAP->setEnabled(false);
	m_ui.cbAR->setEnabled(false);
	m_ui.cbECC->setEnabled(false);
    m_ui.txtMMS->setEnabled(false);
	m_ui.txtMBS->setEnabled(false);
	m_ui.cbFE->setEnabled(false);
	m_ui.cbDFC->setEnabled(false);
	m_ui.cbSRC->setEnabled(false);
	m_ui.cmbRefMode->setEnabled(false);
	m_ui.cmbDynRefMode->setEnabled(false);
	m_ui.cmbIntBasisSet->setEnabled(false);
	m_ui.cmbDynAv->setEnabled(false);
	m_ui.cmbDynAv_W->setEnabled(false);

}

QString InputDlg::GetPathAndSave(
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

void InputDlg::OnBtnOpenWS()
{
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

void InputDlg::OnBtnOpenWU()
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

void InputDlg::LoadFID(QString filename, fid_type_e fid_type)
{
	// retrieve the format of this kind of FID from the combo box
	tarquin::fid_format_e format = 
		static_cast<tarquin::fid_format_e>
		(m_ui.cmbFormat->itemData(m_ui.cmbFormat->currentIndex(), Qt::UserRole).toInt());

	// set the options on the session
	tarquin::Options& opts = m_session->GetWorkspace().GetOptions();
	opts.SetFormat(format);

	try
	{
		if( WATER_SUPPRESSED_FID == fid_type )
		{
			tarquin::CFID& fid = m_session->GetWorkspace().GetFID();

			fid.Load(filename.toStdString(), opts, m_session->GetWorkspace(), GetLog());


			// enable all the buttons
			m_ui.txtSF->setEnabled(true);
			m_ui.txtTF->setEnabled(true);
			m_ui.txtET->setEnabled(true);
			m_ui.txtTE1->setEnabled(true);
			m_ui.txtTM->setEnabled(true);
			m_ui.txtCN->setEnabled(true);
	        m_ui.cmbPS->setEnabled(true);
			m_ui.txtRO->setEnabled(true);
            m_ui.txtParaFile->setEnabled(true);
            m_ui.btnOpenParaFile->setEnabled(true);
            m_ui.btnOpenWU->setEnabled(true);
            m_ui.txtWU->setEnabled(true);
            
            m_ui.txtZF->setEnabled(true);
            m_ui.txtPPMstart->setEnabled(true);
            m_ui.txtPPMend->setEnabled(true);
            m_ui.txtBL->setEnabled(true);
            m_ui.txtLB->setEnabled(true);

            m_ui.btnOpenXML->setEnabled(true);
            m_ui.btnOpenCSV->setEnabled(true);
            m_ui.btnSaveXML->setEnabled(true);
            m_ui.btnSaveLCM->setEnabled(true);

            m_ui.txtBasisCSV->setEnabled(true);
            m_ui.txtBasisXML->setEnabled(true);
            m_ui.txtSaveBasis->setEnabled(true);
            m_ui.txtSaveBasisLCM->setEnabled(true);

            m_ui.txtSP->setEnabled(true);
            m_ui.txtEP->setEnabled(true);
            m_ui.txtWC->setEnabled(true);
            m_ui.txtCW->setEnabled(true);
            m_ui.txtLD->setEnabled(true);
	        m_ui.txtDR->setEnabled(true);
            m_ui.txtZFKS->setEnabled(true);
            m_ui.txtIM->setEnabled(true);
            m_ui.txtMI->setEnabled(true);
            m_ui.txtIB->setEnabled(true);
            m_ui.txtWConc->setEnabled(true);
            m_ui.txtWA->setEnabled(true);
            m_ui.cbCP->setEnabled(true);
	        m_ui.cbLF->setEnabled(true);
	        m_ui.cbKPWFS->setEnabled(true);
            m_ui.cbKF->setEnabled(true);
            m_ui.cbAP->setEnabled(true);
            m_ui.cbAR->setEnabled(true);
	        m_ui.cbECC->setEnabled(true);
            m_ui.txtMMS->setEnabled(true);
            m_ui.txtMBS->setEnabled(true);
            m_ui.cbFE->setEnabled(true);
	        m_ui.cbDFC->setEnabled(true);
	        m_ui.cbSRC->setEnabled(true);
	        m_ui.cmbRefMode->setEnabled(true);
	        m_ui.cmbDynRefMode->setEnabled(true);
            m_ui.cmbIntBasisSet->setEnabled(true);
            m_ui.cmbDynAv->setEnabled(true);
            m_ui.cmbDynAv_W->setEnabled(true);

			// did the FID contain the water data as well?
			if( fid.GetCWF() )
			{
				m_ui.txtWU->setText(filename);
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

        //
        // transfer the acquisition parameters to the dialog
        //
        UpdateDlg();

		m_ui.btnOpenWS->setEnabled(false);
	}
	catch( const tarquin::Exception& e )
	{
		ErrorDialog(this, tr("File I/O Error"), tr("The exact message was: ") + e.what());

		if( WATER_SUPPRESSED_FID == fid_type )
			m_ui.txtWS->setText(QString());
		else
			m_ui.txtWU->setText(QString());
	}
}

void InputDlg::UpdateDlg()
{
    tarquin::CFID& fid = m_session->GetWorkspace().GetFID();
	tarquin::Options& opts = m_session->GetWorkspace().GetOptions();

    if( fid.IsKnownSamplingFrequency() )
        m_ui.txtSF->setText( QString::number(fid.GetSamplingFrequency(),'g',12) );

    if( fid.IsKnownTransmitterFrequency() )
        m_ui.txtTF->setText( QString::number(fid.GetTransmitterFrequency(),'g',12) );

    if( fid.IsKnownEchoTime() )
        m_ui.txtET->setText( QString::number(fid.GetEchoTime(),'g',12) );

    m_ui.txtTE1->setText( QString::number(opts.GetPRESS_TE1(),'g',12) );
    m_ui.txtTM->setText( QString::number(opts.GetSTEAM_TM(),'g',12) );
    m_ui.txtCN->setText( QString::number(opts.GetCPMG_N()) );
    
    m_ui.cmbRefMode->setCurrentIndex( opts.GetRefSignals() );
    m_ui.cmbDynRefMode->setCurrentIndex( opts.GetDynRefSignals() );
    m_ui.cmbIntBasisSet->setCurrentIndex( opts.GetIntBasisSet() );
    m_ui.cmbPS->setCurrentIndex( opts.GetPulSeq() );
    m_ui.cmbDynAv->setCurrentIndex( opts.GetDynAv() );
    m_ui.cmbDynAv_W->setCurrentIndex( opts.GetDynAvW() );

    m_ui.txtNP->setText( QString::number(fid.GetNumberOfPoints()) );
    m_ui.txtRO->setText( QString::number(fid.GetPPMRef(0)) ); // TODO
    m_ui.txtSP->setText( QString::number(opts.GetRangeStart()) );
    m_ui.txtEP->setText( QString::number(opts.GetRangeEnd()) );
    m_ui.txtWC->setText( QString::number(opts.GetWaterWindow()) );
    m_ui.txtCW->setText( QString::number(opts.GetConvWindowWidth()) );
    m_ui.txtLD->setText( QString::number(opts.GetLambda()) );
    m_ui.txtDR->setText( QString::number(opts.GetMaxDRef()) );
	m_ui.txtZFKS->setText( QString::number(opts.GetZfillKspace()) );
    m_ui.txtIM->setText( QString::number(opts.GetNUMERICAL_TOL_INIT_MU()) );
    m_ui.txtMI->setText( QString::number(opts.GetMaxIters()) );
    m_ui.txtIB->setText( QString::number(opts.GetInitBeta()) );
    m_ui.txtWConc->setText( QString::number(opts.GetWConc()) );
    m_ui.txtWA->setText( QString::number(opts.GetWAtt()) );
    m_ui.txtMMS->setText( QString::number( opts.GetMaxMetabShift() ));
    m_ui.txtMBS->setText( QString::number( opts.GetMaxBroadShift() ));
    m_ui.cbCP->setChecked( opts.GetCombinePreproc() );
    m_ui.cbLF->setChecked( opts.GetLipidFilter() );
    m_ui.cbKPWFS->setChecked( opts.GetKeepPreWsShift() );

    m_ui.cbKF->setChecked( opts.GetFilterKspace() );
    m_ui.cbAP->setChecked( opts.GetAutoPhase() );
    m_ui.cbAR->setChecked( opts.GetAutoReference() );
    m_ui.cbECC->setChecked( opts.GetWaterEddy() );

    m_ui.cbFE->setChecked( opts.GetFullEcho() );
    m_ui.cbDFC->setChecked( opts.GetDynFreqCorr() );
    m_ui.cbSRC->setChecked( opts.GetSwapRowCol() );

    m_ui.txtZF->setText( QString::number(opts.GetZF()));
    m_ui.txtPPMstart->setText( QString::number(opts.GetPPMstart()));
    m_ui.txtPPMend->setText( QString::number(opts.GetPPMend()));
    m_ui.txtBL->setText( QString::number(opts.GetBL()));
    m_ui.txtLB->setText( QString::number(opts.Getlb()));

    if ( opts.GetUsePrecompiled() )
    {
        m_ui.txtBasisXML->setText( QString::fromStdString(opts.GetBasisPath()) );
        m_ui.txtBasisCSV->setText( QString::fromStdString("") );
    }
    else
    {
        m_ui.txtBasisCSV->setText( QString::fromStdString(opts.GetBasisPath()) );
        m_ui.txtBasisXML->setText( QString::fromStdString("") );
    }
	m_ui.txtSaveBasis->setText( QString::fromStdString(opts.GetBasisSaveFile()) );
	m_ui.txtSaveBasisLCM->setText( QString::fromStdString(opts.GetBasisSaveFileLCM()) );
}

/*!
 * Select a precompiled basis header file and attempt to load it.
 */
void InputDlg::OnBtnOpenXML()
{
	// the FID must be loaded before we can simulate the basis
	if( !m_loaded_ws )
	{
		InfoDialog(this, tr("Wrong Order"), tr("You must load a FID before selecting the basis."));
		return;
	}
    
    // check the options used
    if( !CheckDlg() )
		return;

	// retrieve the path that was last used from the settings
	QSettings settings;

	QString path = QFileDialog::getOpenFileName(this, 
			tr("Open Basis File (XML/LCModel BASIS)"), 
			settings.value("basis_xml_input_dir", QString()).toString(),
			tr("XML/basis Files (*.xml *.basis);;all files (*.*)"));

	// nothing selected, they must have pushed cancel
	if( !path.size() )
		return;

	// they have loaded the file, so keep the directory for next time
	GetPathAndSave(path, "basis_xml_input_dir");
	tarquin::Options& opts = m_session->GetWorkspace().GetOptions();
	
	// try and load the basis
	try
	{
		tarquin::CBasis& basis = m_session->GetWorkspace().GetBasis();
            
        int len = path.toStdString().size();
        std::string str_end;
        if ( len > 6 )
            str_end = path.toStdString().substr(len-6,6);
        else
            str_end = "";

        if (( str_end == ".basis" ) || ( str_end == ".BASIS" ))
        {
            if( !basis.ReadLCMBasis(path.toStdString(), m_session->GetWorkspace().GetFIDRaw(), opts, GetLog()) )
                throw std::runtime_error("The ReadLCMBasis function generally failed.");
        }
        else
        {
            if( !load_basis(path.toStdString(), basis) )
                throw std::runtime_error("The load_basis function generally failed.");
        }
		
		// check that the basis matches the FID parameters
		if( !basis.check(m_session->GetWorkspace().GetFIDRaw(), GetLog()) )
		{
			ErrorDialog(this, tr("Parameter Mismatch"), 
					tr("The parameters in the basis file do not match those of the FID undergoing analysis.") +
					tr("\nPlease use a different basis file or create one to match your FID parameters."));
			return;
		}

		// set the file on the form
		m_ui.txtBasisXML->setText(path);

		// set the flag
		m_loaded_basis = true;

		opts.SetUsePrecompiled(true);
	}
	catch( const std::exception& e )
	{
		ErrorDialog(this, tr("Error Loading Basis"), tr("The exact message was: ") + e.what());
	}
	
}

void InputDlg::OnBtnOpenCSV()
{
	// retrieve the path that was last used from the settings
	QSettings settings;

	QString path = QFileDialog::getExistingDirectory(this, 
			tr("Choose TARQUIN Basis Directory"), 
			settings.value("basis_csv_input_dir", QString()).toString());

	// nothing selected, they must have pushed cancel
	if( !path.size() )
		return;

	// they have loaded the file, so keep the directory for next time
	settings.setValue("basis_csv_input_dir", path);

	// set the file on the form
	m_ui.txtBasisCSV->setText(path);

	tarquin::Options& opts = m_session->GetWorkspace().GetOptions();
	opts.SetBasisPath(path.toStdString());	
	opts.SetUsePrecompiled(false);
	m_loaded_basis = true;
}

// save the basis set if it has been simulated
void InputDlg::OnBtnSaveXML()
{
	// retrieve the path that was last used from the settings
	QSettings settings;
	
	QString path = QFileDialog::getSaveFileName(this, 
			tr("Save TARQUIN Basis File (.xml)"), 
			settings.value("save_basis", QString()).toString(),
			tr("XML Files (*.xml)"));
	
	// nothing selected, they must have pushed cancel
	if( !path.size() )
		return;
	
	// set the file on the form
	m_ui.txtSaveBasis->setText(path);
	
	// they have specified a file, so keep the directory for next time
	settings.setValue("save_basis", path);
}

// save the basis set if it has been simulated
void InputDlg::OnBtnSaveLCM()
{
	// retrieve the path that was last used from the settings
	QSettings settings;
	
	QString path = QFileDialog::getSaveFileName(this, 
			tr("Save TARQUIN Basis File (.basis)"), 
			settings.value("save_basis_lcm", QString()).toString(),
			tr("basis Files (*.basis)"));
	
	// nothing selected, they must have pushed cancel
	if( !path.size() )
		return;
	
	// set the file on the form
	m_ui.txtSaveBasisLCM->setText(path);
	
	// they have specified a file, so keep the directory for next time
	settings.setValue("save_basis_lcm", path);
}




/*!
 * Select a paramter file
 */
void InputDlg::OnBtnOpenParaFile()
{
	// retrieve the path that was last used from the settings
	QSettings settings;

	QString path = QFileDialog::getOpenFileName(this, 
			tr("Open TARQUIN Parameter File"), 
			settings.value("para_input_dir", QString()).toString(),
			tr("txt Files (*.txt);;all files (*.*)"));

	// nothing selected, they must have pushed cancel
	if( !path.size() )
		return;
	
	// set the file on the form
	m_ui.txtParaFile->setText(path);
	
	// they have loaded the file, so keep the directory for next time
	settings.setValue("para_input_dir", path);

	tarquin::Options& opts = m_session->GetWorkspace().GetOptions();
	tarquin::CFID& fidraw = m_session->GetWorkspace().GetFIDRaw();

	// TODO what happens if this fails?
	QByteArray path_str_by = path.toAscii();
	char* path_str = path_str_by.data();
	int argc = 3;
	char* argv[3];
	argv[0] = (char *)"tarquingui";
	argv[1] = (char *)"--para_file";
	argv[2] = path_str;
	
	if ( !ParseCommandLine(argc, argv, opts, fidraw) )
	{
		QMessageBox::information(this, tr("Error in para file"), tr("Error, the parameter file could not be parsed correctly."));
		return;
	}

    UpdateDlg();
}

bool InputDlg::CheckDlg()
{
    bool ok = false;

	QString water_suppressed   = m_ui.txtWS->text();
	QString water_unsuppressed = m_ui.txtWU->text();

	if( !m_loaded_ws )
	{
		QMessageBox::information(this, tr("Missing FID"), tr("You must provide a path for the water suppressed FID data."));
		return false;
	}

	/*if( !m_loaded_basis )
	{
		QMessageBox::information(this, tr("Missing Basis"), tr("You must load or simulate a basis before continuing."));
		return;
	}*/

	if( ( m_ui.txtBasisCSV->text().toStdString() != "" ) 
		&& ( m_ui.txtBasisXML->text().toStdString() != "" ) )
	{
		QMessageBox::information(this, tr("Basis Set Problem"), tr("Warning : You have specified to simulate a basis and chosen a pre-simulated basis set.  The pre-simulated basis set will be used."));
		tarquin::Options& opts = m_session->GetWorkspace().GetOptions();
		opts.SetUsePrecompiled(true);
	}

	//
	// extract the options from the dialog 
	//
	double fs = m_ui.txtSF->text().toDouble(&ok);
	if( !ok || fs < 0 )
	{
		QMessageBox::information(this, tr("Field Error"), tr("The sampling frequency must be a positive real number."));
		return false;
	}

	double ft = m_ui.txtTF->text().toDouble(&ok);
	if( !ok || ft < 0 )
	{
		QMessageBox::information(this, tr("Field Error"), tr("The transmitter frequency must be a positive real number."));
		return false;
	}

	double te = m_ui.txtET->text().toDouble(&ok);
	if( !ok || te < 0 )
	{
		QMessageBox::information(this, tr("Field Error"), tr("The echo time must be a non-negative real number."));
		return false;
	}

	double te1 = m_ui.txtTE1->text().toDouble(&ok);
	if( !ok || te1 < 0 )
	{
		QMessageBox::information(this, tr("Field Error"), tr("TE1 must be a non-negative real number."));
		return false;
	}

	double tm = m_ui.txtTM->text().toDouble(&ok);
	if( !ok || tm < 0 )
	{
		QMessageBox::information(this, tr("Field Error"), tr("TM must be a non-negative real number."));
		return false;
	}
    
    int cn = m_ui.txtCN->text().toInt(&ok);
	if( !ok || cn < 1 )
	{
		QMessageBox::information(this, tr("Field Error"), tr("TM must be a real integer greater than 0."));
		return false;
	}

	double ppmref = m_ui.txtRO->text().toDouble(&ok);
	if( !ok )
	{
		QMessageBox::information(this, tr("Field Error"), tr("The PPM reference offset real number, probably 4.65."));
		return false;
	}
    
    int n_start = m_ui.txtSP->text().toInt(&ok);
	if( !ok || n_start < 1 )
	{
		QMessageBox::information(this, tr("Field Error"), tr("The start point must be an integer greater than 0."));
		return false;
	}
    
    int n_end = m_ui.txtEP->text().toInt(&ok);
	if( !ok || n_end < 1 )
	{
		QMessageBox::information(this, tr("Field Error"), tr("The end point must be an integer greater than 0."));
		return false;
	}
    
    double water_cutoff = m_ui.txtWC->text().toDouble(&ok);
	if( !ok )
	{
		QMessageBox::information(this, tr("Field Error"), tr("The cutoff must be numeric."));
		return false;
	}

    int conv_width = m_ui.txtCW->text().toInt(&ok);
	if( !ok || conv_width < 0 ) // TODO
	{
		QMessageBox::information(this, tr("Field Error"), tr("The convolution width must be an integer greater or equal to 0."));
		return false;
	}

    double lambda = m_ui.txtLD->text().toDouble(&ok);
	if( !ok || lambda < 0 )
	{
		QMessageBox::information(this, tr("Field Error"), tr("Lambda must me a non-negative numeric value."));
		return false;
	}

    double max_dref = m_ui.txtDR->text().toDouble(&ok);
	if( !ok || max_dref < 0 )
	{
		QMessageBox::information(this, tr("Field Error"), tr("Max dref must me a non-negative numeric value."));
		return false;
	}

    int zfks = m_ui.txtZFKS->text().toInt(&ok);
	if( !ok || zfks <= 0 )
	{
		QMessageBox::information(this, tr("Field Error"), tr("The k-space zero filling factor must be an integer greater than 0."));
		return false;
	}

    double init_mu = m_ui.txtIM->text().toDouble(&ok);
	if( !ok || init_mu < 0 )
	{
		QMessageBox::information(this, tr("Field Error"), tr("Init mu must me a non-negative numeric value."));
		return false;
	}
    
    int max_iters = m_ui.txtMI->text().toInt(&ok);
	if( !ok || max_iters < 0 )
	{
		QMessageBox::information(this, tr("Field Error"), tr("Maximum iterations must me a non-negative numeric value."));
		return false;
	}

    double init_beta = m_ui.txtIB->text().toDouble(&ok);
	if( !ok || init_beta < 0 )
	{
		QMessageBox::information(this, tr("Field Error"), tr("Init beta must be a non-negative numeric value."));
		return false;
	}
    
    double water_conc = m_ui.txtWConc->text().toDouble(&ok);
	if( !ok || water_conc < 0 )
	{
		QMessageBox::information(this, tr("Field Error"), tr("Water conc must be a non-negative numeric value."));
		return false;
	}

    double water_att = m_ui.txtWA->text().toDouble(&ok);
	if( !ok || water_att < 0 )
	{
		QMessageBox::information(this, tr("Field Error"), tr("Water att must be a non-negative numeric value."));
		return false;
	}

    double max_metab_shift = m_ui.txtMMS->text().toDouble(&ok);
	if( !ok || max_metab_shift < 0 )
	{
		QMessageBox::information(this, tr("Field Error"), tr("Max metabolite shift must be a non-negative numeric value."));
		return false;
	}

    double max_broad_shift = m_ui.txtMBS->text().toDouble(&ok);
	if( !ok || max_broad_shift < 0 )
	{
		QMessageBox::information(this, tr("Field Error"), tr("Max broad shift must be a non-negative numeric value."));
		return false;
	}

    int zf = m_ui.txtZF->text().toInt(&ok);
    if( !ok || zf <= 0 )
	{
		QMessageBox::information(this, tr("Field Error"), tr("The zero filling factor must be an integer greater than 0."));
		return false;
	}

    double ppm_start = m_ui.txtPPMstart->text().toDouble(&ok);
	if( !ok )
	{
		QMessageBox::information(this, tr("Field Error"), tr("PPM start must be a numeric value."));
		return false;
	}

    double ppm_end = m_ui.txtPPMend->text().toDouble(&ok);
	if( !ok || ppm_end < ppm_start )
	{
		QMessageBox::information(this, tr("Field Error"), tr("PPM end must be a numeric value greater than PPM start."));
		return false;
	}

    int bl = m_ui.txtBL->text().toInt(&ok);
    if( !ok || bl <= 2 )
	{
		QMessageBox::information(this, tr("Field Error"), tr("The baseline points must be an intger values greater than 2."));
		return false;
	}
    
    double lb = m_ui.txtLB->text().toDouble(&ok);
	if( !ok || lb < 0 )
	{
		QMessageBox::information(this, tr("Field Error"), tr("Line-broadening must be a non-negative numeric value."));
		return false;
	}

	//
	// transfer the data to the FID
	//
	tarquin::CFID& fidws = m_session->GetWorkspace().GetFID();
	tarquin::CFID& fidwu = m_session->GetWorkspace().GetFIDWater();

	fidws.SetSamplingFrequency(fs);
	fidws.SetTransmitterFrequency(ft);
	fidws.SetEchoTime(te);
	fidws.SetPPMRef(ppmref);
	fidws.SetFilename(water_suppressed.toStdString());

	// FIXME: the water reference should be coupled to the FID data somehow
	// do we even need to set this?
	if( water_unsuppressed.size() )
	{
		fidwu.SetSamplingFrequency(fs);
		fidwu.SetTransmitterFrequency(ft);
		fidwu.SetEchoTime(te);
		fidwu.SetPPMRef(ppmref);
		fidwu.SetFilename(water_unsuppressed.toStdString());
	}

	// options that may be useful later
	tarquin::Options& opts = m_session->GetWorkspace().GetOptions();
    
    opts.SetRef(ppmref);
    // this flag stopps the reference from being overwritten when 
    // reloading the fid
    opts.SetRefSpec(true);

	opts.SetPRESS_TE1(te1);
	opts.SetSTEAM_TM(tm);
	opts.SetCPMG_N(cn);
	opts.SetFilenameWater(water_unsuppressed.toStdString());
	opts.SetFilename(water_suppressed.toStdString());
	opts.SetBasisSaveFile(m_ui.txtSaveBasis->text().toStdString());
	opts.SetBasisSaveFileLCM(m_ui.txtSaveBasisLCM->text().toStdString());
    opts.SetRangeStart(n_start);
    opts.SetRangeEnd(n_end);
    opts.SetWaterWindow(water_cutoff);
    opts.SetConvWindowWidth(conv_width);
    opts.SetMaxDRef(max_dref);
    opts.SetLambda(lambda);
    opts.SetInitMu(init_mu);
	opts.SetZfillKspace(zfks);
    opts.SetMaxIters(max_iters);
    opts.SetInitBeta(init_beta);
    opts.SetWConc(water_conc);
    opts.SetWAtt(water_att);
    opts.SetCombinePreproc( m_ui.cbCP->isChecked() );
    opts.SetLipidFilter( m_ui.cbLF->isChecked() );
    opts.SetKeepPreWsShift( m_ui.cbKPWFS->isChecked() );
    opts.SetFullEcho( m_ui.cbFE->isChecked() );
    
    opts.SetFilterKspace( m_ui.cbKF->isChecked() );
    opts.SetAutoPhase( m_ui.cbAP->isChecked() );
    opts.SetAutoReference( m_ui.cbAR->isChecked() );
    opts.SetWaterEddy( m_ui.cbECC->isChecked() );

    if ( opts.GetSwapRowCol() != m_ui.cbSRC->isChecked() )
    {
        fidws.swap_row_col();
        /*
        int f_rows = opts.GetFitRows();
        int f_cols = opts.GetFitCols();
        opts.SetFitRows(f_cols);
        opts.SetFitCols(f_rows);
        */

        if ( opts.GetFilenameWater() != "" || m_session->GetWorkspace().GetFIDRaw().GetCWF() )
            fidwu.swap_row_col();
    }
    opts.SetSwapRowCol( m_ui.cbSRC->isChecked() );
    opts.SetDynFreqCorr( m_ui.cbDFC->isChecked() );
    opts.SetMaxMetabShift(max_metab_shift);
    opts.SetMaxBroadShift(max_broad_shift);

	opts.SetZF(zf);
	opts.SetPPMstart(ppm_start);
	opts.SetPPMend(ppm_end);
	opts.SetBL(bl);
	opts.Setlb(lb);

    // find ref signals
	tarquin::chem_shift_ref_e ref_signals = static_cast<tarquin::chem_shift_ref_e> (m_ui.cmbRefMode->itemData(m_ui.cmbRefMode->currentIndex(), Qt::UserRole).toInt());
    opts.SetRefSignals(ref_signals);

    // find ref signals
	tarquin::chem_shift_ref_e dyn_ref_signals = static_cast<tarquin::chem_shift_ref_e> (m_ui.cmbRefMode->itemData(m_ui.cmbDynRefMode->currentIndex(), Qt::UserRole).toInt());
    opts.SetDynRefSignals(dyn_ref_signals);
    
    // set basis set
    tarquin::basis_set_e int_basis = static_cast<tarquin::basis_set_e> (m_ui.cmbIntBasisSet->itemData(m_ui.cmbIntBasisSet->currentIndex(), Qt::UserRole).toInt());
    opts.SetIntBasisSet(int_basis);
    
    // set pul seq
    tarquin::pul_seq_e pul_seq = static_cast<tarquin::pul_seq_e> (m_ui.cmbPS->itemData(m_ui.cmbPS->currentIndex(), Qt::UserRole).toInt());
    opts.SetPulSeq(pul_seq);

    // set dynamic averaging scheme
    tarquin::dyn_av_e dyn_av = static_cast<tarquin::dyn_av_e> (m_ui.cmbDynAv->itemData(m_ui.cmbDynAv->currentIndex(), Qt::UserRole).toInt());
    opts.SetDynAv(dyn_av);

    // set dynamic averaging scheme (water ref)
    tarquin::dyn_av_e dyn_av_w = static_cast<tarquin::dyn_av_e> (m_ui.cmbDynAv_W->itemData(m_ui.cmbDynAv_W->currentIndex(), Qt::UserRole).toInt());
    opts.SetDynAvW(dyn_av_w);

    // reload the fids to account for dynamic scans or full echo processing
    // or k-space zero-filling
    fidws.Load(water_suppressed.toStdString(), opts, m_session->GetWorkspace(), GetLog()); // TODO???
	fidws.SetSamplingFrequency(fs);
	fidws.SetTransmitterFrequency(ft);
    fidws.SetEchoTime(te);
    if ( water_unsuppressed.toStdString() != "" )
    {
        if ( m_session->GetWorkspace().GetFIDRaw().GetCWF() )
            fidwu.LoadW(water_unsuppressed.toStdString(), opts, GetLog());
        else
            fidwu.Load(water_unsuppressed.toStdString(), opts, m_session->GetWorkspace(), GetLog());

        fidwu.SetSamplingFrequency(fs);
        fidwu.SetTransmitterFrequency(ft);
		fidwu.SetEchoTime(te);
    }

	return true;
}

void InputDlg::OnBtnFit()
{
    if( CheckDlg() )
    {
        // if CSI ask about voxels
	    tarquin::CFID& fidws = m_session->GetWorkspace().GetFID();
        if ( fidws.GetVoxelCount() > 1 )
        {
            voxel_select voxel_select_dlg(this, m_session);
            if( QDialog::Accepted == voxel_select_dlg.exec() )
            {
                m_session->data_loaded = true;
                m_start_fit = true;
                QDialog::accept();
            }
        }
        else
        {
            m_session->data_loaded = true;
            m_start_fit = true;
            QDialog::accept();
        }
    }
}

void InputDlg::OnBtnLoad()
{
    if( CheckDlg() )
    {
        tarquin::CFID& fidws = m_session->GetWorkspace().GetFID();
        if ( fidws.GetVoxelCount() > 1 )
        {
            voxel_select voxel_select_dlg(this, m_session);
            if( QDialog::Accepted == voxel_select_dlg.exec() )
            {
                m_session->data_loaded = true;
                QDialog::accept();
            }
        }
        else
        {
                m_session->data_loaded = true;
                QDialog::accept();
        }
    }
}

bool InputDlg::startFit() const
{
	return m_start_fit;
}
