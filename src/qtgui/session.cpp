#include "proto_buf.hpp"
#include "session.h"
#include "mainwindow.h"
#include "tarquinplotwidget.h"
#include "resultsdlg.h"
#include "td_conv_ws.hpp"
#include "common.hpp"
#include "CFID.hpp"
#include <sstream> 
#include <QtGui>
#include "guicommon.h"
#include <qwt_plot_curve.h>
#include <qwt_plot_marker.h>
#include <qwt_plot_zoomer.h>
#include <qwt_scale_engine.h>
#include <qwt_legend.h>
#include <qwt_legend_item.h>

#include "common.hpp"

// session constructor
Session::Session(MainWindow* parent) :
	QObject(parent),
	m_parent(parent)
{
    m_fit_number = 0;
    data_loaded = false;
    data_preprocessed = false;
    data_fitted = false;
    mri_loaded = false;
    grid_trans = 70;
}


void Session::add_timedomain_signal(
		const QString&           label, //!< label to plot with FID
		const QColor&            colour,//!< colour to draw line in
		const float              dt,    //!< time interval (recip. of sampling frequency)
		const int                ns,    //!< start point (inclusive)
		const int                ne,    //!< end point (inclusive),
		const cvm::cvector&      y,      //!< the signal to plot
		show_curve_type          is_visible
		)
{
	// get the widget we are adding the plot to
	QwtPlot* plot = m_parent->GetPlot();

	// no attributes on scale
	QwtScaleEngine* engine = plot->axisScaleEngine(QwtPlot::xBottom);
	engine->setAttribute(QwtScaleEngine::Inverted, false);

    
	//
	// add the time domain data
	//
	
	// create the curve
	QwtPlotCurve* curve_td = new QwtPlotCurve(label);
	curve_td->setPaintAttribute(QwtPlotCurve::ClipPolygons);

    if ( m_show_flags.AA )
        curve_td->setRenderHint(QwtPlotItem::RenderAntialiased);

	QPen pen = QPen(colour);
	pen.setWidthF(1.0);
	curve_td->setPen(QPen(colour));
	curve_td->setYAxis(QwtPlot::yLeft);
	curve_td->attach(plot);

	// the length of the subsignal we are plotting
	const int signal_length = ne-ns+1;

	// generate the signal
	cvm::rvector signal(signal_length);
	if( REAL == m_show_flags.mode )
	{
		//plot->setTitle(label + QObject::tr("(real part)"));

		for( int n = ns; n <= ne; ++n )
			signal[n+1] = y[n+1].real();
	}
	else if( IMAG == m_show_flags.mode )
	{
		//plot->setTitle(label + QObject::tr("(imaginary part)"));

		for( int n = ns; n <= ne; ++n )
			signal[n+1] = y[n+1].imag();
	}
	else if( ABS == m_show_flags.mode )
    {
		for( int n = ns; n <= ne; ++n )
			signal[n+1] = pow(pow(y[n+1].real(),2.0) + pow(y[n+1].imag(),2.0),0.5)  ;
    }

	// create the plot data
    
    // generate the time axis
    if ( SEC == m_show_flags.units_td )
    {
        cvm::rvector xaxis(signal_length);
        for( int n = ns; n <= ne; ++n )
            xaxis[n+1] = n*dt;
        curve_td->setData(&xaxis[1], &signal[1], signal_length);
    }
    else if ( TD_PTS == m_show_flags.units_td )
    {
        cvm::rvector xaxis(signal_length);
        for( int n = ns; n <= ne; ++n )
            xaxis[n+1] = n+1;
        curve_td->setData(&xaxis[1], &signal[1], signal_length);
    }

	//curve_td->setVisible( DISPLAY_CURVE == is_visible );
	curve_td->setVisible(true);
}


void Session::add_frequencydomain_signal(
		const QString&           label,    //!< label to use for the series
		const QColor&            colour,   //!< colour to draw the series in
		cvm::rvector&            hz,       //!< Hz scale for x axis, may not be used
		cvm::rvector&            ppm,      //!< ppm scale for x axis, may not be used
		const cvm::cvector&      y,        //!< the data to take the DFT of
		show_curve_type          is_visible
		)
{
	assert( y.size() == ppm.size() );
	assert( y.size() == hz.size() );

	// take the FFT
	cvm::cvector Y(y.size());
	tarquin::fft(y, Y);

	// 0Hz is in the middle
	Y = tarquin::fftshift(Y);

	// take the real part (we want a stride of 1 element so do it manually)
	cvm::rvector YR(Y.size());
	cvm::rvector PTS(Y.size());
	for( int n = 1; n <= YR.size(); ++n )
    {
        if( REAL == m_show_flags.mode )
            YR[n] = Y[n].real();
        if( IMAG == m_show_flags.mode )
            YR[n] = Y[n].imag();
        if( ABS == m_show_flags.mode )
            YR[n] = pow(pow(Y[n].real(),2.0) + pow(Y[n].imag(),2.0),0.5);
        PTS[n] = n;
    }

	// get the widget we are adding the plot to
	QwtPlot* plot = m_parent->GetPlot();

	// no reverse scaling unless we turn it on later
	QwtScaleEngine* engine = plot->axisScaleEngine(QwtPlot::xBottom);
	engine->setAttribute(QwtScaleEngine::Inverted, false);

	// create the curve
	QwtPlotCurve* curve_fd = new QwtPlotCurve(label);
	curve_fd->setPaintAttribute(QwtPlotCurve::ClipPolygons);

    if ( m_show_flags.AA )
        curve_fd->setRenderHint(QwtPlotItem::RenderAntialiased);

	QPen pen = QPen(colour);
	pen.setWidthF(1.0);
	curve_fd->setPen(pen);
	curve_fd->setYAxis(QwtPlot::yLeft);
	curve_fd->attach(plot);


	// PPM on the x axis?
	if( PPM == m_show_flags.units_fd )
	{
		curve_fd->setData(&ppm[1], &YR[1], YR.size());
		plot->setAxisTitle(QwtPlot::xBottom, QObject::tr("Chemical Shift (ppm)"));

		// flip sign of axis for retarded convention
		engine->setAttribute(QwtScaleEngine::Inverted, true);
	}
	// Hz on the x axis?
	else if( HERTZ == m_show_flags.units_fd )
 
	{
		curve_fd->setData(&hz[1], &YR[1], YR.size());
		plot->setAxisTitle(QwtPlot::xBottom, QObject::tr("Frequency (Hz)"));
	}
	else if( FD_PTS == m_show_flags.units_fd )
 
	{
		curve_fd->setData(&PTS[1], &YR[1], YR.size());
		plot->setAxisTitle(QwtPlot::xBottom, QObject::tr("Frequency / (points)"));
	}
	
	curve_fd->setVisible(true);
	//curve_fd->setVisible( DISPLAY_CURVE == is_visible );
}

/*!
 * This function is called after the state has been changed, it will show whatever signals
 * have currently been set to show.
 */
void Session::Update(bool all/*=true*/)
{
	TarquinPlotWidget* plot = m_parent->GetPlot();
	plot->clear();

	tarquin::Options& options = m_workspace.GetOptions();
	const std::vector<tarquin::coord>& fit_list = options.GetFitList();

    // TODO, a better way?
    if ( m_show_flags.WUS && ( options.GetFilenameWater() == "" && !m_workspace.GetFIDRaw().GetCWF() ) )
    {
        m_show_flags.WUS = false;
		ErrorDialog(m_parent, tr("Warning"), tr("Warning. Water unsuppressed file not loaded."));
    }

	//
	// are we drawing in the time domain?
	//
	if( TIME_DOMAIN == m_show_flags.domain )
	{
		assert( m_workspace.GetFIDRaw().GetVectorFID().size() );

		// what is the time step?
        float dt;
        if ( !m_show_flags.WUS )
            dt = 1.0f / m_workspace.GetFIDRaw().GetSamplingFrequency();
        else
            dt = 1.0f / m_workspace.GetFIDWater().GetSamplingFrequency();

		// input signal
		if ( !data_preprocessed || m_show_flags.show_input || m_show_flags.WUS || !m_in_fit_list )
        {
            tarquin::CFID fid;
            if ( !m_show_flags.WUS )
            {
                fid = m_workspace.GetFIDRaw();
            }
            else
            {
                fid = m_workspace.GetFIDWater();
            }

			cvm::cvector         y   = fid.GetVectorFID(m_voxel);

            lb(y,m_workspace.GetFIDRaw(),options.Getlb());
            tarquin::ZeroPad(y,options.GetZF());
			const int            N   = y.size();

			// add it to the plot
			add_timedomain_signal(QObject::tr("original signal"), Qt::blue, dt, 0, N-1, y, m_show_flags.show_input ? DISPLAY_CURVE : HIDE_CURVE);
		}

		// draw the preprocessed signal
		if( m_show_flags.show_preproc && !m_show_flags.WUS && m_in_fit_list )
		{
			const tarquin::CFID& fid = m_workspace.GetFIDProc();
			cvm::cvector  y   = fid.GetVectorFID(m_voxel);
            lb(y,m_workspace.GetFIDProc(),options.Getlb());
            tarquin::ZeroPad(y,options.GetZF());
			const int            N   = y.size();

            //lb(y,m_workspace.GetFIDRaw(),m_show_flags.LB);
            //ZeroPad(y);
			
            // add it to the plot
			add_timedomain_signal(QObject::tr("preprocessed signal"), Qt::red, dt, 0, N-1, y);
		}

		// draw the estimate of the signal (the model)
		if( m_show_flags.show_model && !m_show_flags.WUS && m_in_fit_list )
		{
			tarquin::cvec_stdvec yhat_vec = m_workspace.GetSignalEstimate();
			cvm::cvector yhat = yhat_vec[m_fit_number];
            
            lb(yhat,m_workspace.GetFIDRaw(),options.Getlb());
            tarquin::ZeroPad(yhat,options.GetZF());
			const int    N    = yhat.size();

			// add it to the plot
			add_timedomain_signal(QObject::tr("model signal"), Qt::magenta, dt, 0, N-1, yhat);
		}

		// for each of the basis vectors to display
		const tarquin::cmat_stdvec& S_vec = m_workspace.GetBasisMatrix();
		const cvm::cmatrix& S = S_vec[m_fit_number];

        const tarquin::rvec_stdvec& a_vec = m_workspace.GetAmplitudes();
        const cvm::rvector&         a     = a_vec[m_fit_number];

		for( size_t i = 0; i < m_show_flags.basis_indices.size(); ++i )
		{
			const int    idx = m_show_flags.basis_indices[i];
			cvm::cvector v   = S(idx+1)*a(idx+1);
            lb(v,m_workspace.GetFIDRaw(),options.Getlb());
            tarquin::ZeroPad(v,options.GetZF());
			const int    N   = v.size();

			QString name = QString::fromStdString( m_workspace.GetBasis().GetSignalName(idx) );

			// add it to the plot
			add_timedomain_signal(name, Qt::green, dt, 0, N-1, v);
		}

		// set some titles
        if ( SEC == m_show_flags.units_td )
            plot->setAxisTitle(QwtPlot::xBottom, "time / s");
        else if ( TD_PTS == m_show_flags.units_td )
            plot->setAxisTitle(QwtPlot::xBottom, "time / data points");

		plot->setAxisTitle(QwtPlot::yLeft,   "amplitude / a.u.");
	}
	//
	// are we drawing in the frequency domain?
	//
	else if( FREQUENCY_DOMAIN == m_show_flags.domain )
	{
		// all the signals will be the same length, so we only need
		// to compute the PPM and Hz scales once
        //std::cout << m_fit_number << std::endl;
        //std::cout << m_fit_number << std::endl;
        //std::cout << std::endl << fit_list.size() << std::endl << std::flush;
        //std::cout << m_workspace.GetFIDRaw().GetVectorFID(m_fit_number).size() << std::endl;

		cvm::cvector y   = m_workspace.GetFIDRaw().GetVectorFID(m_voxel);

        // line broaden
        lb(y,m_workspace.GetFIDRaw(),options.Getlb());
        // zero fill
		tarquin::ZeroPad(y,options.GetZF());
		cvm::rvector hz  = m_workspace.GetFID().GetFreqScale(options.GetZF());

		// draw the input signal
		if ( !data_preprocessed || m_show_flags.show_input || m_show_flags.WUS || !m_in_fit_list )
		{
			// zero pad and generate the x axis
            if ( !m_show_flags.WUS )
            {
                cvm::cvector y = m_workspace.GetFIDRaw().GetVectorFID(m_voxel);
		        cvm::rvector ppm = m_workspace.GetFIDRaw().GetPPMScale( m_voxel, options.GetZF() );
                lb(y,m_workspace.GetFIDRaw(),options.Getlb());
			    tarquin::ZeroPad(y,options.GetZF());
			    // add it to the plot
			    add_frequencydomain_signal(QObject::tr("original signal"), Qt::blue, hz, ppm, y, m_show_flags.show_input ? DISPLAY_CURVE : HIDE_CURVE);
            }
            else
            {
                cvm::cvector y = m_workspace.GetFIDWater().GetVectorFID(m_voxel);
		        cvm::rvector ppm = m_workspace.GetFIDRaw().GetPPMScale( m_voxel, options.GetZF() );
                //tarquin::plot(y);
                lb(y,m_workspace.GetFIDWater(),options.Getlb());
			    tarquin::ZeroPad(y,options.GetZF());
                // add it to the plot
			    add_frequencydomain_signal(QObject::tr("original signal"), Qt::blue, hz, ppm, y, m_show_flags.show_input ? DISPLAY_CURVE : HIDE_CURVE);
            }
		}

		// draw the input signal
		if( m_show_flags.show_preproc && !m_show_flags.WUS && m_in_fit_list )
		{
			// zero pad and generate the x axis
		    cvm::rvector ppm = m_workspace.GetFID().GetPPMScale( m_voxel, options.GetZF() );
			cvm::cvector yp = m_workspace.GetFID().GetVectorFID( m_voxel );
            lb(yp,m_workspace.GetFID(),options.Getlb());
			tarquin::ZeroPad(yp,options.GetZF());

			// add it to the plot
			add_frequencydomain_signal(QObject::tr("processed signal"), Qt::black, hz, ppm, yp);
		}

		// draw the model signal
		cvm::cvector baseline(y.size());
		if( m_show_flags.show_model && !m_show_flags.WUS && m_in_fit_list )
		{
			// zero pad and generate the x axis
			tarquin::cvec_stdvec yhat_vec = m_workspace.GetSignalEstimate();
			cvm::cvector yhat = yhat_vec[m_fit_number];
            lb(yhat,m_workspace.GetFID(),options.Getlb());
			tarquin::ZeroPad(yhat,options.GetZF());

			// draw the residual and baseline
			cvm::cvector yp = m_workspace.GetFID().GetVectorFID( m_voxel );
            lb(yp,m_workspace.GetFID(),options.Getlb());
			tarquin::ZeroPad(yp,options.GetZF());
			//add_frequencydomain_signal(QObject::tr("processed signal"), Qt::red, hz, ppm, yp);

			cvm::cvector YHAT(yhat.size());
			//ZeroPad(yhat);
			tarquin::fft(yhat, YHAT);
			// 0Hz is in the middle
			YHAT = tarquin::fftshift(YHAT);

			cvm::cvector YP(yp.size());
			tarquin::fft(yp, YP);
			// 0Hz is in the middle
			YP = tarquin::fftshift(YP);
	
			// estimate the baseline
			cvm::cvector RES = YP - YHAT;
			cvm::cvector BASELINE;
            if ( options.GetBL() > YHAT.size() - 11 )
                options.SetBL(YHAT.size() - 11);

			tarquin::td_conv_ws( RES, BASELINE, options.GetBL()*options.GetZF(), 10);	

			// FIXME FIXME FIXME
			// This code does not belong here, it is doing calculations, i.e. results. It should be 
			// part of the results workspace. This routine should be fast, in case we have to redraw the
			// plot.
			//

			// adjust the offset for RES
			// find points corresponding to 0.2 and 4.0 ppm

		    cvm::rvector ppm = m_workspace.GetFID().GetPPMScale( m_voxel, options.GetZF() );

			int left = 1, right = ppm.size();
            // note, n doesn't go to the end since n+1 points are required
			for ( int n = 1; n < ppm.size(); n++ ) 
			{
				if ( ( ppm(n) >  options.GetPPMend() ) && ( ppm(n+1) <= options.GetPPMend() ) )
					left = n;
				if ( ( ppm(n) >  options.GetPPMstart() ) && ( ppm(n+1) <= options.GetPPMend() ) )
					right = n;
			}

			double Ymax = 0;
			// find max of y for plotting residual	
			for ( int n = left; n < right+1; n++ ) 
			{
                if ( m_show_flags.SB == true )
                {
                    if ( REAL == m_show_flags.mode )
                    {
                        if ( Ymax < YP(n).real() - BASELINE(n).real())
                            Ymax = YP(n).real() - BASELINE(n).real();
                        if ( Ymax < YHAT(n).real() )
                            Ymax = YHAT(n).real();
                    }
                    else if ( IMAG == m_show_flags.mode )
                    {
                        if ( Ymax < YP(n).imag() - BASELINE(n).imag())
                            Ymax = YP(n).imag() - BASELINE(n).imag();
                        if ( Ymax < YHAT(n).imag() )
                            Ymax = YHAT(n).imag();
                    }
                    else if ( ABS == m_show_flags.mode )
                    {
                        if ( Ymax < abs(YP(n)) - abs(BASELINE(n)))
                            Ymax = abs(YP(n)) - abs(BASELINE(n));
                        if ( Ymax < abs(YHAT(n)) )
                            Ymax = abs(YHAT(n));
                    }
                }
                else
                {
                    if ( REAL == m_show_flags.mode )
                    {
                        if ( Ymax < YP(n).real() )
                            Ymax = YP(n).real();
                        if ( Ymax < YHAT(n).real()+BASELINE(n).real() )
                            Ymax = YHAT(n).real()+BASELINE(n).real();
                    }
                    else if ( IMAG == m_show_flags.mode )
                    {
                        if ( Ymax < YP(n).imag() )
                            Ymax = YP(n).imag();
                        if ( Ymax < YHAT(n).imag()+BASELINE(n).imag() )
                            Ymax = YHAT(n).imag()+BASELINE(n).imag();
                    }
                    else if ( ABS == m_show_flags.mode )
                    {
                        if ( Ymax < abs(YP(n)))
                            Ymax = abs(YP(n));
                        if ( Ymax < abs(YHAT(n)) + abs(BASELINE(n)) )
                            Ymax = abs(YHAT(n)) + abs(BASELINE(n));
                    }
                }
			}

			double res_min = 0;
			// find min of residual for plotting residual	
			for ( int n = left; n < (right+1); n++ ) 
            {
                if ( REAL == m_show_flags.mode )
                {
                    if ( res_min > RES(n).real() - BASELINE(n).real() )
                        res_min = RES(n).real() - BASELINE(n).real() ;
                }
                else if ( IMAG == m_show_flags.mode )
                {
                    if ( res_min > RES(n).imag() - BASELINE(n).imag() )
                        res_min = RES(n).imag() - BASELINE(n).imag() ;
                }
                else if ( ABS == m_show_flags.mode )
                {
                    if ( res_min > abs(RES(n)) - abs(BASELINE(n)))
                        res_min = abs(RES(n)) - abs(BASELINE(n));
                }
            }

            
			for ( int n = 1; n < (ppm.size()+1); n++ ) 
            {
                if ( REAL == m_show_flags.mode )
                    RES(n) = tarquin::tcomplex(RES(n).real() - BASELINE(n).real() + Ymax - res_min, 0);
                else if ( IMAG == m_show_flags.mode )
                    RES(n) = tarquin::tcomplex(0,RES(n).imag() - BASELINE(n).imag() + Ymax - res_min);
                else if ( ABS == m_show_flags.mode )
                    RES(n) = tarquin::tcomplex(abs(RES(n)) - abs(BASELINE(n)) + Ymax - res_min,0);
            }
			
			baseline(BASELINE.size());
			BASELINE = tarquin::fftshift(BASELINE);
			tarquin::ifft(BASELINE, baseline);

            
            if ( m_show_flags.SB == true && !m_show_flags.WUS && m_in_fit_list )
            {
                // add baseline corrected signal to the plot
                add_frequencydomain_signal(QObject::tr("bc signal"), Qt::black, hz, ppm, yp-baseline);
            }

			// add yhat to the plot
            if ( m_show_flags.SB == true )
                add_frequencydomain_signal(QObject::tr("model signal"), Qt::red, hz, ppm, yhat);
            else
                add_frequencydomain_signal(QObject::tr("model signal"), Qt::red, hz, ppm, yhat+baseline);

			cvm::cvector res(RES.size());
			RES = tarquin::fftshift(RES);
			tarquin::ifft(RES, res);

			add_frequencydomain_signal(QObject::tr("residual"), Qt::black, hz, ppm, res);
		    
            add_frequencydomain_signal(QObject::tr("baseline"), Qt::green, hz, ppm, baseline);
		
			// for each of the basis vectors to display
			const tarquin::cmat_stdvec& S_vec = m_workspace.GetBasisMatrix();
			const cvm::cmatrix&         S     = S_vec[m_fit_number];
			const tarquin::rvec_stdvec& a_vec = m_workspace.GetAmplitudes();
			const cvm::rvector&         a     = a_vec[m_fit_number];

			for( size_t i = 0; i < m_show_flags.basis_indices.size(); ++i )
			{
				const int    idx = m_show_flags.basis_indices[i];
				cvm::cvector v   = S(idx+1) * a[idx+1];
				lb(v,m_workspace.GetFID(),options.Getlb());
				tarquin::ZeroPad(v,options.GetZF());

				QString name = QString::fromStdString( m_workspace.GetBasis().GetSignalName(idx) );

				// add it to the plot
				if ( m_show_flags.SB == true )
					add_frequencydomain_signal(name, Qt::green, hz, ppm, v);
				else
					add_frequencydomain_signal(name, Qt::green, hz, ppm, v+baseline);
			}

		}

		// set some titles
		//plot->setTitle("DFT (real part)");
		plot->setAxisTitle(QwtPlot::yLeft, "Intensity (a.u.)");
	}
    
    /*
    if ( m_workspace.GetFIDRaw().GetVoxelCount() > 1 )
    {
        std::ostringstream ss;
        ss << "Row : " << m_voxel.row << ", Col : " << m_voxel.col << ", Slice : " << m_voxel.slice;
        std::string title = ss.str();
        plot->setTitle(QString::fromStdString(title));
    }
    else
        plot->setTitle("");
    */

	// put axes back to automatic scaling
	//plot->setAxisAutoScale(QwtPlot::xBottom);
	//plot->setAxisAutoScale(QwtPlot::yLeft);
    
    if ( all )
    {
        // configure the zoomer
        plot->InitZoomer();
    }

    // configure the legend
    plot->InitLegend();

	// refresh the display
	plot->replot();
}


/*!
 * This function is called to bring up the results dialog for the current set of results.
 */
void Session::ShowResults()
{
	ResultsDlg dlg(m_parent, this);
	dlg.InitFromWorkspace(m_workspace);
	dlg.exec();
}

/*!
 * Is fitted data available?
 */
bool Session::FitAvailable()
{
	return m_workspace.GetAmplitudes().size() > 0;
}

/*!
 * Save the workspace as an XML file.
 */
void Session::Save(const QString& filename)
{
	if( !save_results(filename.toStdString().c_str(), m_workspace) )
		throw std::runtime_error("failed to save file");
}

/*!
 * Load the workspace as an XML file.
 */
void Session::Load(const QString& filename)
{
	if( !load_results(filename.toStdString().c_str(), m_workspace) )
		throw std::runtime_error("failed to load file");

	// now turn on all the flags
	m_show_flags.show_input   = true;
	m_show_flags.show_preproc = true;
	m_show_flags.show_model   = true;
	m_show_flags.domain       = Session::FREQUENCY_DOMAIN;
	m_show_flags.units_fd     = Session::PPM;
	// set a default fit number
	m_fit_number = 0;

	// add all the vectors to display
	m_show_flags.basis_indices.clear();
	const tarquin::cmat_stdvec& S_vec = GetWorkspace().GetBasisMatrix();
	const cvm::cmatrix& S = S_vec[m_fit_number];
	for( int i = 0; i < S.nsize(); ++i )
		m_show_flags.basis_indices.push_back(i);

	// and redraw
	Update();

	// hide the processed signal	
	//QwtPlot* plot = m_parent->GetPlot();
	//QwtLegend* legend = plot->legend();
	//QList<QWidget*> items = legend->legendItems();
	//((QwtLegendItem*)items[0])->setChecked(false);
	// does not work - ((QwtPlotItem*)items[0])->setVisible(false);
    //plot->replot();	
	//legend->updateLegend();
	//greenCurve.detach();
	
}

void Session::SetVoxel(const tarquin::coord& voxel)
{
    m_voxel = voxel;
}

const tarquin::coord Session::GetVoxel()
{
    return m_voxel;
}
