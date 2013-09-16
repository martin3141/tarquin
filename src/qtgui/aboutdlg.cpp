#include "aboutdlg.h"
#include "version/version.h"

AboutDlg::AboutDlg(QWidget* parent) :
	QDialog(parent)
{
	m_ui.setupUi(this);

	// set the text for the version and copyright
	m_ui.lblCopyright->setText( QString::fromStdString(tarquin::version::copyright()) );
	m_ui.lblVersion->setText( tr("Version: ") + QString::fromStdString(tarquin::version::version_string()) );
}

