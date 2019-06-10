#include "Workspace.hpp"
#include "Options.hpp"
#include "CFID.hpp"
#include "resultsdlg.h"
#include <QtWidgets>
#include "session.h"

ResultsDlg::ResultsDlg(QWidget* parent, Session* session) :
	QDialog(parent),
    m_session(session),
	m_model(NULL)
{
	assert( session );
	m_ui.setupUi(this);

	// create the empty model
	m_model = new QStandardItemModel(0, 4, this);
	m_ui.treeView->setModel(m_model);

	// aesthetics
	m_ui.treeView->setRootIsDecorated(false);
	m_ui.treeView->setAlternatingRowColors(true);
	m_ui.treeView->header()->setResizeMode(QHeaderView::ResizeToContents);
}

namespace
{

struct ResultRow
{
	int         idx;   //!< index into original data arrays
	std::string name;  //!< name of the basis vector
	float       ampl;  //!< the amplitude to be displayed
	float       crlb;  //!< the Cramer Rao Lower Bound
	float       error; //!< error as a percentage of amplitude

	ResultRow() :
		idx(0), ampl(0.0f), crlb(0.0f), error(0.0f) { }
};

// this needs a nice, unambiguous prototype
int to_lower(int c)
{
	return tolower(c);
}

bool operator< (const ResultRow& lhs, const ResultRow& rhs)
{
	std::string left  = lhs.name;
	std::string right = rhs.name;

	std::transform(left.begin(),  left.end(),  left.begin(),  to_lower);
	std::transform(right.begin(), right.end(), right.begin(), to_lower);

	return left > right;
}

} // namespace

void ResultsDlg::InitFromWorkspace(const tarquin::Workspace& workspace)
{
	// some useful locals
	const tarquin::CBasis&  basis      = workspace.GetBasis();
	const tarquin::rvec_stdvec&	ampl_norm  = workspace.GetAmplitudesNormalised();
	const tarquin::rvec_stdvec&    crlb_norm  = workspace.GetCRLBsNormalised();
	const tarquin::Options& options    = workspace.GetOptions();
	int fit_no = m_session->m_fit_number;

    if ( fit_no > ampl_norm.size() || fit_no < 0 )
    {
        std::cout << "Error, not sure why this happens - but the fit number is invalid." << std::endl;
        std::cout << "Current fit number       : " << fit_no << std::endl;
        std::cout << "Max allowable fit number : " << ampl_norm.size() << std::endl;
        std::cout << "I have set fit number back to 0 to prevent a crash." << std::endl;
        fit_no = 0;
        m_session->m_fit_number = 0;
    }

	assert( ampl_norm[fit_no].size() == crlb_norm[fit_no].size() );

	// these column have the same heading, regardless of scaling
	m_model->setHeaderData(0, Qt::Horizontal, tr("Basis Vector"));
	m_model->setHeaderData(3, Qt::Horizontal, tr("% std. dev."));

	// we are not doing water scaling, so arbitary units
	if( options.GetFilenameWater() == "" )
	{
		m_model->setHeaderData(1, Qt::Horizontal, tr("Amplitude (a.u.)"));
		m_model->setHeaderData(2, Qt::Horizontal, tr("std. dev. (a.u.)"));
	}
	// we are doing water scaling
	else
	{
		m_model->setHeaderData(1, Qt::Horizontal, tr("Amplitude (mM)"));
		m_model->setHeaderData(2, Qt::Horizontal, tr("std. dev. (mM)"));
	}

	std::vector<ResultRow> results;
	// for each basis vector
	for( int i = 0; i < ampl_norm[fit_no].size(); ++i )
	{
		ResultRow row;
		row.idx   = i;
		row.name  = basis.GetSignalName(i);
		row.ampl  = ampl_norm[fit_no][i+1];
		row.crlb  = crlb_norm[fit_no][i+1];
		row.error = row.crlb/row.ampl * 100.0;

		results.push_back(row);
	}

	// sort the rows
	std::sort(results.begin(), results.end());

	// add the sorted data to the model
	for( size_t i = 0; i < results.size(); ++i )
	{
		m_model->insertRow(0);
		m_model->setData(m_model->index(0, 0), QString::fromStdString(results[i].name)  );
		m_model->setData(m_model->index(0, 1), QString::number(results[i].ampl,  'g', 4));
		m_model->setData(m_model->index(0, 2), QString::number(results[i].crlb,  'g', 4));
		m_model->setData(m_model->index(0, 3), QString::number(results[i].error, 'g', 4));
	}

    if ( workspace.GetAmplitudesNormalisedComb().size() > 0 )
    {
        for(int n = 1; n < workspace.GetAmplitudesNormalisedComb()[fit_no].size()+1; n++)
        {
            m_model->insertRow(results.size()+n-1);
            m_model->setData(m_model->index(results.size()+n-1, 0), QString::fromStdString(workspace.GetMetabNamesComb()[n-1])  );
            m_model->setData(m_model->index(results.size()+n-1, 1), QString::number(workspace.GetAmplitudesNormalisedComb()[fit_no][n],  'g', 4));
            m_model->setData(m_model->index(results.size()+n-1, 2), QString::number(workspace.GetCRLBsNormalisedComb()[fit_no][n],  'g', 4));
            m_model->setData(m_model->index(results.size()+n-1, 3), QString::number(workspace.GetCRLBsNormalisedComb()[fit_no][n]/workspace.GetAmplitudesNormalisedComb()[fit_no][n]*100, 'g', 4));

        }
    }

}
