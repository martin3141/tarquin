#include "session.h"
#include "exception.hpp"
#include "formatdlg.h"
#include "guicommon.h"

#include <boost/filesystem.hpp>
#include <QtWidgets>


formatdlg::formatdlg(QWidget *parent, Session* session) :
    QDialog(parent),
	m_session(session)
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
}

void formatdlg::accept()
{
	// retrieve the format of this kind of FID from the combo box
	tarquin::fid_format_e format = 
		static_cast<tarquin::fid_format_e>
		(m_ui.cmbFormat->itemData(m_ui.cmbFormat->currentIndex(), Qt::UserRole).toInt());

	// set the options on the session
	tarquin::Options& opts = m_session->GetWorkspace().GetOptions();
	opts.SetFormat(format);
	QDialog::accept();
}
