#ifndef __COMMON__
#define __COMMON__

#include <stdlib.h>
#include <complex>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <math.h>
#include <fftw3.h>
#include <assert.h>
#include <string.h>
#include "version/version.h"


#define BOOST_SYSTEM_NO_LIB 1

#include <cvm.h>
#include "compat.hpp"

#ifdef _MSC_VER
// this isn't available on windows
inline double round(double x)
{
	return double(int(x + 0.5));
}
#endif


namespace tarquin 
{
    // name of gnuplot binary
    extern std::string g_strGnuPlot;

    // if we are using a toy compiler
#ifdef _MSC_VER

    // terminal type for gnuplot
    const std::string terminal = "windows";	
    const std::string pngterminal = "png truecolor font arial 10";	
    const std::string epsterminal = "postscript eps enhanced color solid";
    // directory / file separator
    // for some reason, this is better of as a forward slash rather than backslash
    const std::string filesep = "/";

#ifndef M_PI
	const double M_PI = 3.14159265358979323846264338327950288;
#endif

#ifndef isnan
#include <float.h>
#define isnan _isnan
#endif

    // no, we are probably using a man's compiler (tm)
#else

    // terminal type for gnuplot
    //const std::string terminal = "wxt size 1024,768";	
    const std::string terminal = "wxt";	
    //const std::string terminal = "x11";	
    const std::string pngterminal = "png truecolor font Vera 10";	
    const std::string epsterminal = "postscript eps enhanced color";

    // directory / file separator
    const std::string filesep = "/";

#endif

    //! the type used to represent real numbers
    typedef cvm::treal treal;

    //! the type used to represent complex numbers
    typedef cvm::tcomplex tcomplex;

    //! the type used to represent integers 
    typedef int integer;

    //! The tolerance to use for determining what is a zero in NNLS calculations.
    #define NUMERICAL_TOL_NNLS 1e-11

    #define NUMERICAL_TOL_INF_JE     1e-30   // stopping threshold for infinity norm of J.'*e
    #define NUMERICAL_TOL_L2_RES     1e-30  // stopping threshold for l2 norm of e
    #define NUMERICAL_TOL_DIFF_DELTA 1e-30   // only appears to be used for numerical Jacobian


    //
    // Common small functions.
    //	

    //! Get extension of file, empty string if not present.
    inline std::string GetFilenameExtension(std::string strFilename)
    {
	std::string strExt;

	if( strFilename.rfind(".", strFilename.size()) != std::string::npos ) 
	    strExt = strFilename.substr(strFilename.rfind(".",strFilename.size()));

	return strExt;
    }

	inline std::size_t GetFileSize(std::string strFilename)
	{
		std::size_t nStart, nEnd;

		std::ifstream file(strFilename.c_str());

		if( true == file.fail() )
			return 0;

		nStart = size_t(file.tellg());
		file.seekg(0, std::ios::end);
		nEnd = size_t(file.tellg());

		return nEnd-nStart;
	}

    //! Get part of file name up to extension.
    inline std::string GetFilenameBase(std::string strFilename)
    {
	std::string strExt = GetFilenameExtension(strFilename);

	return strFilename.substr(0, strFilename.size() - strExt.size());
    }

    //! Loose trailing / leading white space from a string.
    inline void LooseWhiteSpace(std::string& str)
    {
	std::string strNew;

	for(std::size_t n = 0; n < str.size(); n++)
	    if( (str.at(n) != ' ') && (str.at(n) != '\t') )
		strNew += str.at(n);

	str = strNew;
    }

    //! Returns true if the file with the fully pathed name exists.
    inline bool FileExists(const std::string& strFilename)
    {
	std::ifstream file(strFilename.c_str());
	file.close();

	return !file.fail(); 
    }


#define BYTESWAP(x) ByteSwap((unsigned char *) &x,sizeof(x))

    inline void ByteSwap(unsigned char * b, int n)
    {
	register int i = 0;
	register int j = n-1;
	while (i<j)
	{
	    std::swap(b[i], b[j]);
	    i++, j--;
	}
    }
    

    inline treal stdev(const cvm::rvector& x)
    {
        int N = x.size();
        // find the mean of x
        treal mean = 0;
        for ( int n = 1; n < N + 1; n++)
            mean += x(n);
        mean = mean / N;
        // find the standard deviation of x
        treal stdev = 0;
        for ( int n = 1; n < N + 1; n++)
            stdev += ( x(n) - mean ) * ( x(n) - mean );
        stdev = sqrt(stdev/N);
        return stdev;
    }
    
    inline treal stdev(const cvm::rvector& x, int start, int end)
    {
        int N = end - start;
        // find the mean of x
        treal mean = 0;
        for ( int n = start; n < end + 1; n++)
            mean += x(n);
        mean = mean / N;
        // find the standard deviation of x
        treal stdev = 0;
        for ( int n = start; n < end + 1; n++)
            stdev += ( x(n) - mean ) * ( x(n) - mean );
        stdev = sqrt(stdev/N);
        return stdev;
    }

    inline void plot(cvm::rvector& freq_sig)
    {	
	// write data to file 
	std::ofstream data_outfile("plot.txt");

	for(int m = 1; m < freq_sig.size() + 1; m++) {
	    data_outfile << freq_sig(m);
	    data_outfile << std::endl;
	}

	// write gnuplot script to file 
	std::ofstream gnuplot_outfile("gnuplot.txt");
	//gnuplot_outfile << "set term " << terminal << std::endl;
	gnuplot_outfile << "bind \"r\" \"set xrange [] reverse;replot\"" << std::endl;
	gnuplot_outfile << "set xrange [] " << std::endl;
	gnuplot_outfile	<< "plot 'plot.txt' using 1 with lines";
	gnuplot_outfile << std::endl;
	gnuplot_outfile << "pause -1";
	gnuplot_outfile << std::endl;

	// run script
	std::string strRunMe;
	strRunMe = g_strGnuPlot + " gnuplot.txt";
	int ret = system(strRunMe.c_str());
    (void)ret;
    /*if ( ret != -1 )
        assert(0);
        */
    }
    
    inline void plot(cvm::cvector& freq_sig)
    {	
	// write data to file 
	std::ofstream data_outfile("plot.txt");

	for(int m = 1; m < freq_sig.size() + 1; m++) {
	    data_outfile << freq_sig(m).real();
	    data_outfile << std::endl;
	}

	// write gnuplot script to file 
	std::ofstream gnuplot_outfile("gnuplot.txt");
	//gnuplot_outfile << "set term " << terminal << std::endl;
	gnuplot_outfile << "bind \"r\" \"set xrange [] reverse;replot\"" << std::endl;
	gnuplot_outfile << "set xrange [] " << std::endl;
	gnuplot_outfile	<< "plot 'plot.txt' using 1 with lines";
	gnuplot_outfile << std::endl;
	gnuplot_outfile << "pause -1";
	gnuplot_outfile << std::endl;

	// run script
	std::string strRunMe;
	strRunMe = g_strGnuPlot + " gnuplot.txt";
	int ret = system(strRunMe.c_str());
    (void)ret;
    /*if ( ret != -1 )
        assert(0);
        */
    }


    inline void plot(cvm::rvector& scale, cvm::cvector& freq_sig)
    {	
	// write data to file 
	std::ofstream data_outfile("plot.txt");

	for(int m = 1; m < freq_sig.size() + 1; m++) {
	    data_outfile << scale(m);
	    data_outfile << " " << freq_sig(m).real();
	    data_outfile << std::endl;
	}

	// write gnuplot script to file 
	std::ofstream gnuplot_outfile("gnuplot.txt");
	//gnuplot_outfile << "set term " << terminal << std::endl;
	gnuplot_outfile << "bind \"r\" \"set xrange [] reverse;replot\"" << std::endl;
	gnuplot_outfile << "set xrange [] reverse" << std::endl;
	gnuplot_outfile	<< "plot 'plot.txt' using 1:2 with lines";
	gnuplot_outfile << std::endl;
	gnuplot_outfile << "pause -1";
	gnuplot_outfile << std::endl;

	// run script
	std::string strRunMe;
	strRunMe = g_strGnuPlot + " gnuplot.txt";
	int ret = system(strRunMe.c_str());
    (void)ret;
    //if ( ret != -1 )
        //assert(0);
    }
    
    inline void plot(cvm::rvector& scale, cvm::rvector& freq_sig)
    {	
	// write data to file 
	std::ofstream data_outfile("plot.txt");

	for(int m = 1; m < freq_sig.size() + 1; m++) {
	    data_outfile << scale(m);
	    data_outfile << " " << freq_sig(m);
	    data_outfile << std::endl;
	}

	// write gnuplot script to file 
	std::ofstream gnuplot_outfile("gnuplot.txt");
	//gnuplot_outfile << "set term " << terminal << std::endl;
	gnuplot_outfile << "bind \"r\" \"set xrange [] reverse;replot\"" << std::endl;
	gnuplot_outfile << "set xrange [] reverse" << std::endl;
	gnuplot_outfile	<< "plot 'plot.txt' using 1:2 with lines";
	gnuplot_outfile << std::endl;
	gnuplot_outfile << "pause -1";
	gnuplot_outfile << std::endl;

	// run script
	std::string strRunMe;
	strRunMe = g_strGnuPlot + " gnuplot.txt";
	int ret = system(strRunMe.c_str());
    (void)ret;
    //if ( ret != -1 )
        //assert(0);
    }


    inline void plot(cvm::rvector& scale, cvm::cmatrix& freq_sig, std::vector<std::string> labels)
    {	
	// write data to file 
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
	gnuplot_outfile << "bind \"r\" \"set xrange [] reverse;replot\"" << std::endl;
	//gnuplot_outfile << "set xrange [] reverse" << std::endl;
	gnuplot_outfile << "set xrange [] reverse" << std::endl;
	//gnuplot_outfile << "set xrange [] writeback" << std::endl;
	//gnuplot_outfile << "set noautoscale x" << std::endl;
	gnuplot_outfile	<< "plot 'plot.txt' using 1:2 with lines title '";
	gnuplot_outfile << labels[0] << "'" ;
	for (int n = 2; n < freq_sig.nsize()+1; n++) 
	    gnuplot_outfile	<< ", 'plot.txt' using 1:" << n+1 << " with lines title '" << labels[n-1] << "'";
	gnuplot_outfile << std::endl;
	gnuplot_outfile << "pause -1";
	gnuplot_outfile << std::endl;

	// run script
	std::string strRunMe;
	strRunMe = g_strGnuPlot + " gnuplot.txt";
	int ret = system(strRunMe.c_str());
    if ( ret != -1 )
        assert(0);
    }

    inline void plot(cvm::rvector& scale, cvm::cmatrix& freq_sig)
    {	
	// write data to file 
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
	gnuplot_outfile << "set xrange [] reverse" << std::endl;
	gnuplot_outfile << "bind \"r\" \"set xrange [] reverse;replot\"" << std::endl;
	gnuplot_outfile	<< "plot 'plot.txt' using 1:2 with lines";
	for (int n = 2; n < freq_sig.nsize()+1; n++) 
	    gnuplot_outfile	<< ", 'plot.txt' using 1:" << n+1 << " with lines";
	gnuplot_outfile << std::endl;
	gnuplot_outfile << "pause -1";
	gnuplot_outfile << std::endl;

	// run script
	std::string strRunMe;
	strRunMe = g_strGnuPlot + " gnuplot.txt";
	int ret = system(strRunMe.c_str());
    if ( ret != -1 )
        assert(0);
    }

    inline void plot(cvm::cmatrix& freq_sig)
    {	
	// write data to file 
	std::ofstream data_outfile("plot.txt");
	for(int m = 1; m < freq_sig.msize() + 1; m++) {
	    data_outfile << m;
	    for(int n = 1; n < freq_sig.nsize() + 1; n++) {
		data_outfile << " " << freq_sig(m,n).real();
	    }
	    data_outfile << std::endl;
	}

	// write gnuplot script to file 
	std::ofstream gnuplot_outfile("gnuplot.txt");
	//gnuplot_outfile << "set term " << terminal << std::endl;
	//gnuplot_outfile << "set xrange [] reverse" << std::endl;
	gnuplot_outfile << "bind \"r\" \"set xrange [] reverse;replot\"" << std::endl;
	gnuplot_outfile	<< "plot 'plot.txt' using 1:2 with lines";
	for (int n = 2; n < freq_sig.nsize()+1; n++) 
	    gnuplot_outfile	<< ", 'plot.txt' using 1:" << n+1 << " with lines";
	gnuplot_outfile << std::endl;
	gnuplot_outfile << "pause -1";
	gnuplot_outfile << std::endl;

	// run script
	std::string strRunMe;
	strRunMe = g_strGnuPlot + " gnuplot.txt";
	int ret = system(strRunMe.c_str());
    if ( ret != -1 )
        assert(0);
    }

    inline void plot(cvm::rmatrix& freq_sig)
    {	
	// write data to file 
	std::ofstream data_outfile("plot.txt");
	for(int m = 1; m < freq_sig.msize() + 1; m++) {
	    data_outfile << m;
	    for(int n = 1; n < freq_sig.nsize() + 1; n++) {
		data_outfile << " " << freq_sig(m,n);
	    }
	    data_outfile << std::endl;
	}

	// write gnuplot script to file 
	std::ofstream gnuplot_outfile("gnuplot.txt");
	//gnuplot_outfile << "set term " << terminal << std::endl;
	gnuplot_outfile << "set xrange [:] reverse" << std::endl;
	gnuplot_outfile << "bind \"r\" \"set xrange [] reverse;replot\"" << std::endl;
	gnuplot_outfile	<< "plot 'plot.txt' using 1:2 with lines";
	for (int n = 2; n < freq_sig.nsize()+1; n++) 
	    gnuplot_outfile	<< ", 'plot.txt' using 1:" << n+1 << " with lines";
	gnuplot_outfile << std::endl;
	gnuplot_outfile << "pause -1";
	gnuplot_outfile << std::endl;

	// run script
	std::string strRunMe;
	strRunMe = g_strGnuPlot + " gnuplot.txt";
	int ret = system(strRunMe.c_str());
	(void)ret;
    }

    /*! Plot two FIDs on the same graph.
     * \param fidA is the first FID.
     * \param fidB is the second FID.
     */
    class CFID;
    void PlotFIDs(const CFID& fidA, const CFID& fidB);

    inline float GeneratePosGauss(float mean, float stdev) 
    {
	float x1, x2, w, y;
	do {
	    do {
		x1 = 2.0f * rand()/RAND_MAX - 1.0f;
		x2 = 2.0f * rand()/RAND_MAX - 1.0f;
		w = x1 * x1 + x2 * x2;
	    } while ( w >= 1.0f );
	    w = std::sqrt( (-2.0f * std::log( w ) ) / w );
	    y = x1 * w * stdev + mean;
	} while ( y <= 0);
	return y;	
    }

    inline float GenerateGauss(float mean, float stdev) 
    {
	float x1, x2, w, y;
	do {
	    x1 = 2.0f * rand()/RAND_MAX - 1.0f;
	    x2 = 2.0f * rand()/RAND_MAX - 1.0f;
	    w = x1 * x1 + x2 * x2;
	} while ( w >= 1.0f );
	w = std::sqrt( (-2.0f * std::log( w ) ) / w );
	y = x1 * w * stdev + mean;
	return y;
    }

    inline cvm::cvector GenerateNoise(int size, float stdev) 
    {
	cvm::cvector noise(size);
	std::complex<float> j(0,1);
	for ( int n = 1; n < size + 1; n++ )
	    noise(n) = GenerateGauss(0, stdev) + j*GenerateGauss(0, stdev);

	return noise;
    }

    // for deserialisation
    /*template <typename T> bool
	inline load_xml(T &s, std::string strFilename) 
	{
	    // open the archive
	    std::ifstream ifs(strFilename.c_str());

	    if( 0 == ifs.good() )
		return false;

	    //boost::archive::xml_iarchive ia(ifs);
	    boost::archive::text_iarchive ia(ifs);
	    //boost::archive::binary_iarchive ia(ifs);

	    // restore from the archive
	    ia >> BOOST_SERIALIZATION_NVP(s);		

	    return true;
	}*/

    // for serialisation
    /*template <class T>
	inline bool save_xml(const T &s, const char * filename)
	{
	    // make an archive
	    std::ofstream ofs(filename);

	    if( 0 == ofs.good() )
		return false;

	    //boost::archive::xml_oarchive oa(ofs);
	    boost::archive::text_oarchive oa(ofs);
	    //boost::archive::binary_oarchive oa(ofs);
	    oa << BOOST_SERIALIZATION_NVP(s);

	    return true;
	}
	*/

    void plotfit(
			cvm::rvector& scale, 
			cvm::cmatrix& freq_sig, 
			std::vector<int>& plot_index, 
			treal ppm_start, 
			treal ppm_end,
			std::string output_filename="plot.eps",
			bool pause_on_screen=true
			);

    inline void savefitimag(cvm::rvector& scale, cvm::cmatrix& freq_sig, std::string strFilename)
    {	
	// write data to file 
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
	gnuplot_outfile << "set term pdf size 12.7cm, 8.9cm" << std::endl;
	//gnuplot_outfile << "set term " << epsterminal << std::endl;
	//gnuplot_outfile << "set term " << pngterminal << std::endl;
	gnuplot_outfile << "set output '" << strFilename << "'" << std::endl;
	gnuplot_outfile << "set grid" << std::endl;
	gnuplot_outfile << "unset ytics" << std::endl;
	gnuplot_outfile << "set key off" << std::endl;
	gnuplot_outfile << "set xtics nomirror" << std::endl;
	//gnuplot_outfile << "set border 1" << std::endl;
	//gnuplot_outfile << "set xtics 0.2" << std::endl;
	//gnuplot_outfile << "set mxtics 2" << std::endl;
	//gnuplot_outfile << "set grid xtics mxtics lw 0.5, lw 0.5" << std::endl;
	//gnuplot_outfile << "set grid xtics mxtics ls 0, ls 9" << std::endl;
	gnuplot_outfile << "set xtic out" << std::endl;
	//gnuplot_outfile << "set grid xtics mxtics" << std::endl;
	gnuplot_outfile << "set xlabel \"Chemical Shift (ppm)\"" << std::endl;
	//gnuplot_outfile << "set xrange [] reverse" << std::endl;
	//gnuplot_outfile << "set xrange [0.0 : 6.0] reverse" << std::endl;
	//gnuplot_outfile << "set xrange [-0.2 : 6.0] reverse" << std::endl;
	gnuplot_outfile << "set xrange [0.2 : 4.0] reverse" << std::endl;
    
    //gnuplot_outfile	<< "plot 'plot.txt' using 1:3 with lines lw 6 lc rgb 'red' title \"fit\", ";
    gnuplot_outfile	<< "plot 'plot.txt' using 1:3 with lines lw 6 lc rgb 'dark-grey' title \"fit\", ";
	gnuplot_outfile	<< "'plot.txt' using 1:2 with lines lw 2 lc rgb 'black' title \"data\", ";
	gnuplot_outfile	<< "'plot.txt' using 1:4 with lines lw 2 lc rgb 'black' title \"residual\", ";
	gnuplot_outfile	<< "'plot.txt' using 1:5 with lines lw 1 lc rgb 'black' title \"zero res\", ";
	gnuplot_outfile	<< "'plot.txt' using 1:6 with lines lw 1 lc rgb 'black', ";
	gnuplot_outfile	<< "'plot.txt' using 1:7 with lines lw 1 lc rgb 'black', ";
	gnuplot_outfile	<< "'plot.txt' using 1:8 with lines lw 2 lc rgb 'black'";

	//for (int n = 5; n < freq_sig.nsize()+1; n++) 
	//    gnuplot_outfile	<< ", 'plot.txt' using 1:" << n+1 << " with lines lt -1 title \"metab\"";

	gnuplot_outfile << std::endl;
	// run script
	std::string strRunMe;
	strRunMe = g_strGnuPlot + " gnuplot.txt";
	int ret = system(strRunMe.c_str());
    if ( ret != -1 )
        assert(0);
    }

    inline void savepdffit(cvm::rvector& scale, cvm::cmatrix& freq_sig, std::string strFilename, std::ostringstream& table, std::string title, const std::vector<std::string>& names, bool ext_output, double ppm_start, double ppm_end, double cex) 
    {	
    
	// write data to file 
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
	gnuplot_outfile << "set term pdfcairo size 29.7cm,21.0cm font \"Arial," << 6*cex << "\"" << std::endl;
	gnuplot_outfile << "set bmargin at screen 0.1" << std::endl;
	gnuplot_outfile << "set lmargin at screen 0.05" << std::endl;
	gnuplot_outfile << "set tmargin at screen 0.9" << std::endl;
	gnuplot_outfile << "set rmargin at screen 0.72" << std::endl;
	gnuplot_outfile << "set output '" << strFilename << "'" << std::endl;
	gnuplot_outfile << "set grid" << std::endl;
	gnuplot_outfile << "unset ytics" << std::endl;
	gnuplot_outfile << "set key off" << std::endl;
	gnuplot_outfile << "set xtics nomirror" << std::endl;
	gnuplot_outfile << "set xtic out" << std::endl;
	gnuplot_outfile << "set xtics 0.2" << std::endl;
	gnuplot_outfile << "set xlabel \"Chemical Shift (ppm)\"" << std::endl;
    gnuplot_outfile << "set label \"" << table.str() << "\" at graph(1.02),graph(0.99) font \"Courier," << 6*cex << "\"" << std::endl;

    gnuplot_outfile << "set label \"TARQUIN version " << version::version_string() << "\" at graph(0.0),graph(1.03)" << std::endl;
    
    gnuplot_outfile << "set label \"" << title << "\" at screen(0.5),graph(1.05) center font \"Arial,"<< 11*cex <<"\"" << std::endl;

    gnuplot_outfile << "set timestamp top" << std::endl;
    gnuplot_outfile << "set format x '%1.1f'" << std::endl;

	//gnuplot_outfile << "set xrange [" << ppm_start << ":" << ppm_end << "] reverse" << std::endl;
    gnuplot_outfile << "set xrange [" << ppm_end << ":" << ppm_start << "]" << std::endl;
	gnuplot_outfile	<< "plot 'plot.txt' using 1:3 with lines lw 6 lc rgb 'red' title \"fit\", ";
	gnuplot_outfile	<< "'plot.txt' using 1:2 with lines lw 2 lc rgb 'black' title \"data\", ";
	gnuplot_outfile	<< "'plot.txt' using 1:4 with lines lw 2 lc rgb 'black' title \"residual\", ";
	gnuplot_outfile	<< "'plot.txt' using 1:5 with lines lw 1 lc rgb 'black' title \"zero res\", ";
	gnuplot_outfile	<< "'plot.txt' using 1:6 with lines lw 1 lc rgb 'black', ";
	gnuplot_outfile	<< "'plot.txt' using 1:7 with lines lw 1 lc rgb 'black', ";
	gnuplot_outfile	<< "'plot.txt' using 1:8 with lines lw 2 lc rgb 'black'";

    gnuplot_outfile << std::endl;
    if ( ext_output )
    {
        for (int n = 8; n < freq_sig.nsize()+1; n++) 
        {
            gnuplot_outfile << "reset" << std::endl;
            gnuplot_outfile << "set bmargin at screen 0.1" << std::endl;
            gnuplot_outfile << "set lmargin at screen 0.05" << std::endl;
            gnuplot_outfile << "set tmargin at screen 0.9" << std::endl;
            gnuplot_outfile << "set rmargin at screen 0.72" << std::endl;
            gnuplot_outfile << "set grid" << std::endl;
            gnuplot_outfile << "unset ytics" << std::endl;
            gnuplot_outfile << "set key off" << std::endl;
            gnuplot_outfile << "set xtics nomirror" << std::endl;
            gnuplot_outfile << "set xtic out" << std::endl;
            gnuplot_outfile << "set xtics 0.2" << std::endl;
            gnuplot_outfile << "set xlabel \"Chemical Shift (ppm)\"" << std::endl;
            gnuplot_outfile << "set label \"" << table.str() << "\" at graph(1.02),graph(0.99) font \"Courier," << 6*cex << "\"" << std::endl;

            gnuplot_outfile << "set label \"TARQUIN version " << version::version_string() << "\" at graph(0.0),graph(1.03)" << std::endl;

            gnuplot_outfile << "set label \"" << names[n-8] << "\" at screen(0.5),graph(1.05) center font \"Arial," << 11*cex << "\"" << std::endl;

            gnuplot_outfile << "set timestamp top" << std::endl;
            gnuplot_outfile << "set format x '%1.1f'" << std::endl;

            //gnuplot_outfile << "set xrange [" << ppm_start << ":" << ppm_end << "] reverse" << std::endl;
            gnuplot_outfile << "set xrange [" << ppm_end << ":" << ppm_start << "]" << std::endl;

            gnuplot_outfile	<< "plot 'plot.txt' using 1:2 with lines lw 2 lc rgb 'black' title \"data\", ";
            gnuplot_outfile	<< "'plot.txt' using 1:" << n+1 << " with lines lw 6 lc rgb 'red' title \"metab\", ";
	        gnuplot_outfile	<< "'plot.txt' using 1:8 with lines lw 2 lc rgb 'black'";
            gnuplot_outfile << std::endl;
        }
    }

	gnuplot_outfile << std::endl;
	// run script
	std::string strRunMe;
	strRunMe = g_strGnuPlot + " gnuplot.txt";
	int ret = system(strRunMe.c_str());
    (void)ret;
    //if ( ret != -1 )
    //    assert(0);
    }

    //void ExportToTXT(std::string& strFile, Workspace& workspace);

    inline void str2rvec(const std::string& str_in, std::vector<double>& vec_out)
    {
        vec_out.clear();
        std::string chunk;
        std::stringstream stream(str_in);
        while (getline(stream, chunk, '\\' ))
        {
            std::istringstream istr_row(chunk);
            double tmp;
            istr_row >> tmp;
            vec_out.push_back(tmp);
        }
    }

    inline void rvec2str(const std::vector<double>& vec_in, std::string& str_out)
    {
        std::ostringstream szero;
        szero << vec_in[0];
        str_out = szero.str();
        for ( size_t n = 1; n < vec_in.size(); n++ )
        {
            std::ostringstream sn;
            sn << vec_in[n];
            str_out += "\\" + sn.str();
        }
    }

    void ZeroPad(cvm::cvector& y, int factor);
    void lb(cvm::cvector& y, const CFID& fid, double lb);
    
} // namespace tarquin






#endif
