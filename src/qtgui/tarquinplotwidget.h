#ifndef Tarquin_PLOT_WIDGET_INCLUDED
#define Tarquin_PLOT_WIDGET_INCLUDED

#include <qwt_plot.h>

class QwtPlotZoomer;
class QwtPlotPanner;
class QwtPlotPicker;

class TarquinPlotWidget : public QwtPlot
{
	Q_OBJECT 

	public:

		TarquinPlotWidget(QWidget* parent);

		void InitZoomer();

		void InitLegend();
		
        void DisableZoomer();

        void EnableZoomer();

	private slots:

		void ToggleCurve(QwtPlotItem* item, bool on);

	private:

		//! The zoomer we are using.
		QwtPlotZoomer* m_zoomer;

		//! The panner we are using.
		QwtPlotPanner* m_panner;

	public:

		//! The picker we are using.
        QwtPlotPicker* m_picker;
		
};

#endif // Tarquin_PLOT_WIDGET_INCLUDED
