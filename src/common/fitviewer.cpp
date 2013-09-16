#include <iostream>
#include "CBasis.hpp"
#include "Workspace.hpp"
#include "tarquin.hpp"
#include "cvm_util.hpp"
#include "preprocess.hpp"
#include "common.hpp"
#include "td_conv_ws.hpp"
#include "fitviewer.hpp"
#include "proto_buf.hpp"
#include "export_data.hpp"

using namespace std;

int tarquin::FitView(
		std::string strFilename,  //!< name of XML file to plot
		std::string strPlotSigs, 
		treal ppm_start, 
		treal ppm_end, 
		treal lb,
		bool pause 
		)
{
	// create a new workspace
	Workspace workspace;		

	// deserialise 
	/*if( false == load_xml(workspace, strFilename) ) 
	{
		std::cerr << "\nerror: couldn't load file: \"" << strFilename << "\"\n" << std::endl;
		return -1;
	}*/
    
    if( false == load_results(strFilename, workspace) ) 
	{
		std::cerr << "\nerror: couldn't load file: \"" << strFilename << "\"\n" << std::endl;
		return -1;
	}

	//workspace.GenerateNonSTLstruc();

	return FitView(workspace, strPlotSigs, ppm_start, ppm_end, lb, pause, "plot.eps");
}

int tarquin::FitView(
		const Workspace& workspace, 
		std::string strPlotSigs, 
		treal ppm_start, 
		treal ppm_end, 
		treal lb,
		bool pause_on_screen,
		std::string output_filename
		)
{
	if ( ppm_start > ppm_end ) 
	{
		treal temp;
		temp = ppm_start;
		ppm_start = ppm_end;
		ppm_end = temp;
	}

	// create a new options 
	Options options = workspace.GetOptions();

	CFID yfid = workspace.GetFID();
	coord_vec fit_list = options.GetFitList();
	cvm::cvector y = yfid.GetVectorFID(fit_list[0]); //TODO
	cvec_stdvec yhat_vec = workspace.GetSignalEstimate();
	cvm::cvector yhat = yhat_vec[0]; // TODO
    
    // test code... 
    /*
    cvm::cvector linesh(yhat.size());
    treal frac;
    for ( int n = 1; n < yhat.size() + 1; n++ ) 
    {
        frac = (treal) n/yhat.size();
        linesh(n) = y(n)/(yhat(n)*0.999+y(n)*0.001);
    }
    
    cvm::cvector linesh_smo(yhat.size());
    td_conv_ws( linesh, linesh_smo, yfid.GetSamplingFrequency()/10, 10);	
    //plot(linesh);
    //plot(linesh_smo);
    
    for ( int n = 1; n < yhat.size() + 1; n++ ) 
        yhat(n) = yhat(n) * linesh_smo(n);
    */
            
    int zf = 1;
    // zero fill if less than 4096 points
    if ( y.size() < 4096 )
        zf = int(round(4096/y.size()));

    y.resize(y.size()*zf);
    yhat.resize(yhat.size()*zf);


	//apply broadening to y and yhat
	int N = y.size();
	//int N = yfid.GetNumberOfPoints();

	// the sampling interval (time step)
	treal dt = 1.0 / yfid.GetSamplingFrequency(); 

	// generate time signal
	cvm::rvector t(N);
	for(integer n = 0; n < N; n++)
		t(n+1) = n*dt;

	// treal lb = 5;
	// treal lb = options.Getlb();
	// std::cout << lb << std::endl;
	treal alpha = -lb*3.141593;

	for(integer n = 0; n < N; n++) {
		y(n+1) = y(n+1)*exp(alpha*t(n+1));
		yhat(n+1) = yhat(n+1)*exp(alpha*t(n+1));
	}



	cvm::cvector Y = fft(y);
	cvm::cvector YHAT = fft(yhat);



	Y = fftshift(Y);
	YHAT = fftshift(YHAT);

	cvm::cvector residual = y-yhat;
	cvm::cvector baseline = residual;

	//	MakeBaseline(yfid, yhat, baseline);
	//	exit(-1);
	/* old method
	// set baseline points not considered in fitting to be zero	
	// for ( int n = options.GetRangeStart(); n < baseline.size()+1-11; n ++ )
	//for ( int n = 1; n < baseline.size()+1; n ++ )
	for ( int n = options.GetRangeStart(); n < baseline.size()+1-11; n ++ )
	baseline(n) = 0;

	cvm::cvector BASELINE = fft(baseline);
	BASELINE = fftshift(BASELINE);
	*/

	cvm::cvector RESIDUAL = fft(residual);
	RESIDUAL = fftshift(RESIDUAL);

	cvm::cvector BASELINE;
	//td_conv_ws( RESIDUAL, BASELINE, 100, 10);	
	td_conv_ws( RESIDUAL, BASELINE, 50, 10);	

	cmat_stdvec s_vec = workspace.GetBasisMatrix(); // TODO should be a reference?
	cvm::cmatrix s = s_vec[0]; // TODO

    s.resize(s.msize()*zf, s.nsize());


	CBasis basis = workspace.GetBasis();
	std::vector<std::string> sig_names = basis.GetSignalNames();
	//std::cout << strPlotSigs << std::endl;
	//for ( size_t n = 0; n < sig_names.size(); n++ )
	//   std::cout << sig_names[n] << std::endl;


	std::vector<int> plot_index;

	// if the argument is a star plot all signals 
	if ( strPlotSigs == "*" ) {
		for ( size_t n = 0; n < sig_names.size(); n++ )
			plot_index.push_back(n);
	}

	else {
		// convert strPlotSigs into a vector of strings
		std::string buf; // buffer string
		stringstream ss(strPlotSigs);
		// for each token
		while (ss >> buf) {
			// search through basis signals
			for ( size_t n = 0; n < sig_names.size(); n++ ){
				// until a match is found
				if ( buf == sig_names[n] ) {
					//std::cout << "Found match " << sig_names[n] << " " << buf << " " << n << std::endl;
					plot_index.push_back(n);
				}
			}
		}
	}

	//cvm::cmatrix& s = workspace.GetGroupMatrix();
	rvec_stdvec ahat_vec = workspace.GetAmplitudes();
	cvm::rvector ahat = ahat_vec[0]; // TODO

	//y(n+1) = y(n+1)*exp(alpha*t(n+1));

	// apply amplitude and lb paras
	for ( int n = 1; (n < ahat.size() + 1); n++) {
		for ( int m = 1; (m < s.msize() + 1); m++) {
			s(m,n) = s(m,n)*ahat(n)*exp(alpha*t(m));
		}
	}

	cvm::cmatrix S = fft(s);
	S = fftshift(S);
	coord voxel(1, 1, 1); // TODO
	cvm::rvector freq_scale = yfid.GetPPMScale(voxel, zf);

	cvm::cmatrix all(S.msize(),S.nsize()+7);
	all.assign(1,8,S);

	// find points corresponding to 0.2 and 4.0 ppm
	int left = 1, right = 1;
	for ( int n = 1; n < (freq_scale.size()); n++ ) {
		if ( ( freq_scale(n) > ppm_end ) && ( freq_scale(n+1) < ppm_end ) )
			left = n;
		if ( ( freq_scale(n) > ppm_start ) && ( freq_scale(n+1) < ppm_start ) )
			right = n;
	}
	// std::cout << left << std::endl;
	// std::cout << right << std::endl;

	double Ymax = 0;
	// find max of y for plotting residual	
	for ( int n = left; n < right+1; n++ ) {
		if ( Ymax < Y(n).real() )
			Ymax = Y(n).real();
		if ( Ymax < YHAT(n).real()+BASELINE(n).real() )
			Ymax = YHAT(n).real()+BASELINE(n).real();
	}

	double res_min = 0;
	// find min of residual for plotting residual	
	for ( int n = left; n < (right+1); n++ ) 
		if ( res_min > RESIDUAL(n).real() - BASELINE(n).real() )
			res_min = RESIDUAL(n).real() - BASELINE(n).real() ;

	double res_max = 0;
	// find max of residual for plotting residual	
	for ( int n = left; n < (right+1); n++ ) 
		if ( res_max < RESIDUAL(n).real() - BASELINE(n).real() )
			res_max = RESIDUAL(n).real() - BASELINE(n).real() ;

	for ( int m = 1; m < (S.msize()+1); m++ ) {	
		all(m,1) = Y(m);
		all(m,2) = YHAT(m)+BASELINE(m);
		all(m,3) = RESIDUAL(m) - BASELINE(m) + Ymax - res_min;
		all(m,4) = Ymax-res_min;
		all(m,5) = Ymax;
		all(m,6) = Ymax-res_min+res_max;
		all(m,7) = BASELINE(m);
		for ( int n = 8; n < all.nsize() + 1 ; n++ )
			all(m,n) = all(m,n) + BASELINE(m);
	}

	// output some convergence information first
	std::vector < std::vector < double > >  info = workspace.GetLMinfo();
    // TODO update [0]'s below
	std::cout << "\nl2 norm of error at initial p  = " << info[0][0];
	std::cout << "\nl2 norm of error at final p    = " << info[0][1];
	std::cout << "\nl2 norm of J.'*e at final p    = " << info[0][2];
	std::cout << "\nl2 norm of D*p at final p      = " << info[0][3];
	std::cout << "\nnumber of iterations           = " << info[0][5];
	std::cout << "\nnumber of function evaluations = " << info[0][7];
	std::cout << "\nnumber of Jacobian evaluations = " << info[0][8];

	if( 1 == info[0][6] )
		std::cout << "\nstopped by small gradient\n";
	else if( 2 == info[0][6] )
		std::cout << "\nstopped by small D*p\n";
	else if( 3 == info[0][6] )
		std::cout << "\nstopped by iteration limit\n";
	else if( 4 == info[0][6] )
		std::cout << "\nstopped by singular matrix - restart with bigger mu";
	else if( 5 == info[0][6] )
		std::cout << "\nno further reduction possible - restart with bigger mu\n";
	else if( 6 == info[0][6] )
		std::cout << "\nstopping by small l2 norm of error\n";


    ExportCsvFit("fitview.csv", workspace);

    // plot the results 
    plotfit(freq_scale, all, plot_index, ppm_start, ppm_end, output_filename, pause_on_screen);

	return 0;
}

void tarquin::plotfit(
		cvm::rvector& scale, 
		cvm::cmatrix& freq_sig, 
		std::vector<int>& plot_index, 
		treal ppm_start, 
		treal ppm_end,
		std::string output_filename,
		bool pause_on_screen
		)
{	
	// write data to file 
	// TODO: use mkstemp instead
	std::ofstream data_outfile("plot.txt");
	for(int m = 1; m < freq_sig.msize() + 1; m++) {
		data_outfile << scale(m);
		for(int n = 1; n < freq_sig.nsize() + 1; n++) {
			data_outfile << " " << freq_sig(m,n).real();
		}
		data_outfile << std::endl;
	}

	// write gnuplot script to file 
	std::ofstream gnuplot_outfile("gnuplot.txt");
	//gnuplot_outfile << "set term " << terminal << std::endl;
	gnuplot_outfile << "set grid" << std::endl;
	gnuplot_outfile << "unset ytics" << std::endl;
	gnuplot_outfile << "set key off" << std::endl;
	gnuplot_outfile << "set xtics nomirror" << std::endl;
	//gnuplot_outfile << "set xtics 0.2" << std::endl;
	//gnuplot_outfile << "set mxtics 2" << std::endl;
	//gnuplot_outfile << "set grid xtics mxtics lw 0.5, lw 0.5" << std::endl;
	//gnuplot_outfile << "set grid xtics mxtics ls 0, ls 9" << std::endl;
	gnuplot_outfile << "set xtic out" << std::endl;
	//gnuplot_outfile << "set grid xtics mxtics" << std::endl;
	gnuplot_outfile << "set xlabel \"Chemical Shift (ppm)\"" << std::endl;
	//gnuplot_outfile << "set xrange [0.0 : 6.0] reverse" << std::endl;
	//gnuplot_outfile << "set xrange [-0.2 : 6.0] reverse" << std::endl;
	gnuplot_outfile << "set xrange [" << ppm_start << " : " << ppm_end << "] reverse" << std::endl;
	//gnuplot_outfile	<< "plot 'plot.txt' using 1:3 with lines lt 1 lw 3 title \"fit\", ";
	gnuplot_outfile	<< "plot 'plot.txt' using 1:2 with lines lt -1 title \"data\", ";
	gnuplot_outfile	<< "'plot.txt' using 1:3 with lines lt -1 lc rgb 'red' title \"fit\", ";
	gnuplot_outfile	<< "'plot.txt' using 1:4 with lines lt -1 title \"residual\", ";
	gnuplot_outfile	<< "'plot.txt' using 1:5 with lines lt -1 title \"zero res\", ";
	gnuplot_outfile	<< "'plot.txt' using 1:6 with lines lt -1, ";
	gnuplot_outfile	<< "'plot.txt' using 1:7 with lines lt -1, ";
	gnuplot_outfile	<< "'plot.txt' using 1:8 with lines lt -1";


	for ( size_t n = 0; n < plot_index.size(); n++)
		gnuplot_outfile	<< ", 'plot.txt' using 1:" << plot_index[n]+9 << " with lines lt -1 lc rgb 'blue' title \"metab\"";

	gnuplot_outfile << std::endl;

	//for (int n = 9; n < freq_sig.nsize()+1; n++) 
	//    gnuplot_outfile	<< ", 'plot.txt' using 1:" << n+1 << " with lines lt -1 title \"metab\"";

	gnuplot_outfile << std::endl;

	if( true == pause_on_screen )
	{
		gnuplot_outfile << "pause -1" << std::endl;
	}

	gnuplot_outfile << "set term " << epsterminal << std::endl;
	gnuplot_outfile << "set output '" << output_filename << "'" << std::endl;
	gnuplot_outfile << "replot" << std::endl;

	gnuplot_outfile << std::endl;

	// run script
	std::string strRunMe;
	strRunMe = g_strGnuPlot + " gnuplot.txt";
	int ret = system(strRunMe.c_str());
    if ( ret != -1 )
        assert(0);
}


