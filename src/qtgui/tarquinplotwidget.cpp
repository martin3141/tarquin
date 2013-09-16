#include "tarquinplotwidget.h"
#include <QtGui>

#include <qwt_plot.h>
#include <qwt_plot_marker.h>
#include <qwt_plot_curve.h>
#include <qwt_legend.h>
#include <qwt_data.h>
#include <qwt_text.h>
#include <qwt_math.h>
#include <qwt_plot_zoomer.h>
#include <qwt_plot_magnifier.h>
#include <qwt_plot_panner.h>
#include <qwt_scale_draw.h>
#include <qwt_scale_widget.h>
#include <qwt_plot_layout.h>
#include <qwt_legend_item.h>

TarquinPlotWidget::TarquinPlotWidget(QWidget* parent) :
	QwtPlot(parent),
	m_zoomer(NULL),
	m_panner(NULL),
	m_picker(NULL)
{
	// displayed items may be toggled
	connect(this, SIGNAL(legendChecked(QwtPlotItem*, bool)), this, SLOT(ToggleCurve(QwtPlotItem*, bool))); 
	
	// set the plot background color	
	this->setCanvasBackground(Qt::white);

	// may be useful for changing the color of various parts of the plot
	QPalette pal = this->palette();
	pal.setColor( QPalette::Window, Qt::white );
	pal.setColor( QPalette::WindowText, Qt::black );
	pal.setColor( QPalette::Text, Qt::black );
	this->setPalette( pal );
	this->setAutoFillBackground( true );
	this->setCanvasLineWidth( 0 );
    
    // fonts and sizes
    //this->setAxisFont(QwtPlot::yLeft, QFont("Helvetica", 18));
    //this->setAxisFont(QwtPlot::xBottom, QFont("Helvetica", 18));
    //this->setAxisTitleFont(QwtPlot::yLeft, QFont("Helvetica", 18));
    //this->setAxisTitleFont(QwtPlot::xBottom, QFont("Helvetica", 18));
    //this->setTitleFont(QFont("Helvetica", 18));

}

/*!
 * This needs to be called every time the range of data changes.
 *
 * It sets up the following behaviour:
 *
 *  left button         - zoom rectangle 
 *  middle button       - panning 
 *  right button        - go one back up zoom stack
 *  ctrl + right button - go all the way back up zoom stack
 */
void TarquinPlotWidget::InitZoomer()
{
	if( m_zoomer )
	{
		delete m_zoomer;
		m_zoomer = NULL;

		delete m_panner;
		m_panner = NULL;

		delete m_picker;
		m_picker = NULL;
	}

	// create the zoomer control
    m_zoomer = new QwtPlotZoomer(canvas());
    m_zoomer->setMousePattern(QwtEventPattern::MouseSelect2, Qt::RightButton, Qt::ControlModifier);
    m_zoomer->setMousePattern(QwtEventPattern::MouseSelect3, Qt::RightButton);

	// panning control
    m_panner = new QwtPlotPanner(canvas());
    m_panner->setAxisEnabled(QwtPlot::yRight, false);
    m_panner->setMouseButton(Qt::MidButton);

    // picker
    m_picker = new QwtPlotPicker(QwtPlot::xBottom, QwtPlot::yLeft, canvas());
    m_picker->setTrackerMode(QwtPlotPicker::ActiveOnly);
    // emit the position of clicks on widget
    m_picker->setSelectionFlags(QwtPlotPicker::PointSelection | QwtPlotPicker::ClickSelection);

	// Avoid jumping when labels with more/less digits
	// appear/disappear when scrolling vertically
	const QFontMetrics fm(axisWidget(QwtPlot::yLeft)->font());
	QwtScaleDraw *sd = axisScaleDraw(QwtPlot::yLeft);
	sd->setMinimumExtent( fm.width("100.00") );

	// aesthetics of selection rectangle
    const QColor c(Qt::darkBlue);
    m_zoomer->setRubberBandPen(c);
    m_zoomer->setTrackerPen(c);
}

void TarquinPlotWidget::DisableZoomer()
{
    m_zoomer->setEnabled(false);
}

void TarquinPlotWidget::EnableZoomer()
{
    m_zoomer->setEnabled(true);
}

void TarquinPlotWidget::InitLegend()
{
	QwtLegend* legend = new QwtLegend(this);
	legend->setFont(QFont("Helvetica", 8));
    legend->setItemMode(QwtLegend::CheckableItem);
    insertLegend(legend, QwtPlot::RightLegend);

	// check each item
	QList<QWidget*> items = legend->legendItems();

	for( int i = 0; i < items.size(); ++i )
	{
		((QwtLegendItem*)items[i])->setChecked(true);
	}
}

void TarquinPlotWidget::ToggleCurve(QwtPlotItem* item, bool on)
{
	item->setVisible(on);
	QWidget* w = legend()->find(item);

	if( w && w->inherits("QwtLegendItem") )
		((QwtLegendItem *)w)->setChecked(on);

    replot();	
}

