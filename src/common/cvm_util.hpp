#ifndef __CVMUTIL__
#define __CVMUTIL__

#include <cvm.h>
#include <fstream>
#include <fftw3.h>
#include "common.hpp"

namespace tarquin 
{
	class fft_wrap
	{
		public:

			fft_wrap(const cvm::cvector& in, cvm::cvector& out) :
                m_out(out)
			{
				m_plan = fftw_plan_dft_1d(in.size(), (fftw_complex*)in.get(), (fftw_complex*)out.get(), FFTW_FORWARD, FFTW_ESTIMATE);
			}

			void operator() ()
			{
				fftw_execute(m_plan);
				
                double scale = 1.0; // m_out.size();
				m_out *= scale;
			}

			~fft_wrap()
			{
				fftw_destroy_plan(m_plan);
			}

			fftw_plan    m_plan;
			cvm::cvector& m_out;
	};
    
    class ifft_wrap
	{
		public:

			ifft_wrap(const cvm::cvector& in, cvm::cvector& out) :
				m_out(out)
			{
				m_plan = fftw_plan_dft_1d(in.size(), (fftw_complex*)in.get(), (fftw_complex*)out.get(), FFTW_BACKWARD, FFTW_ESTIMATE);
			}

			void operator() ()
			{
				fftw_execute(m_plan);

				double scale = 1.0 / m_out.size();
				m_out *= scale;
			}

			~ifft_wrap()
			{
				fftw_destroy_plan(m_plan);
			}

			fftw_plan     m_plan;
			cvm::cvector& m_out;
	};


    /*! Do the IFFT of a vector.
     * \param vIn is the in vector (frequency domain)
     * \param vOut is the out vector (time domain)
     */
	inline void ifft(const cvm::cvector& vIn, cvm::cvector& vOut)
	{
		if( vOut.size() != vIn.size() )
			vOut.resize( vIn.size() );

		ifft_wrap func(vIn, vOut);
		func();
	}

	/*! Do the FFT of a vector.
	 * \param vIn is the in vector (time domain)
	 * \param vOut is the out vector (frequency domain)
	 */
	inline void fft(const cvm::cvector& vIn, cvm::cvector& vOut)
	{
		if( vOut.size() != vIn.size() )
			vOut.resize( vIn.size() );

		// this cast may be dangerous on wierd platforms which is probably why static_cast doesn't work
		//fftw_plan plan = fftw_plan_dft_1d(vIn.size(), (fftw_complex*)vIn.get(), (fftw_complex*)vOut.get(), 
		//		FFTW_FORWARD, FFTW_ESTIMATE);

		//fftw_execute(plan);		
		//fftw_destroy_plan(plan);
	
		fft_wrap func(vIn, vOut);
		func();
	}

	inline cvm::cvector fft(const cvm::cvector& in)
	{
		cvm::cvector out;
		fft(in, out);
		return out;
	}
    
    inline cvm::cvector ifft(const cvm::cvector& in)
	{
		cvm::cvector out;
		ifft(in, out);
		return out;
	}

	
	inline cvm::cmatrix fft(cvm::cmatrix& time_sig) 
	{
		cvm::cmatrix freq_sig(time_sig.msize(), time_sig.nsize());

		for (int n = 1 ; n < time_sig.nsize()+1 ; n++) 
		{
			cvm::cvector time_sig_vec(time_sig(n));
			cvm::cvector freq_sig_vec(fft(time_sig_vec));
			// cvm::cvector freq_sig_shift_vec(fftshift(freq_sig_vec));
			// this is for non-fftshfting
			freq_sig(n) = freq_sig_vec;
			//freq_sig(n) = freq_sig_shift_vec;
		}
		return freq_sig;
	}

    inline cvm::cmatrix ifft(cvm::cmatrix& time_sig) 
	{
		cvm::cmatrix freq_sig(time_sig.msize(), time_sig.nsize());

		for (int n = 1 ; n < time_sig.nsize()+1 ; n++) 
		{
			cvm::cvector time_sig_vec(time_sig(n));
			cvm::cvector freq_sig_vec(ifft(time_sig_vec));
			// cvm::cvector freq_sig_shift_vec(fftshift(freq_sig_vec));
			// this is for non-fftshfting
			freq_sig(n) = freq_sig_vec;
			//freq_sig(n) = freq_sig_shift_vec;
		}
		return freq_sig;
	}

	

    
	template<typename vector_t> vector_t fftshift(const vector_t& vec_in) 
	{
		int p = 1;
		vector_t vec_out(vec_in.size());

		int cutoff = (int) ceil(vec_in.size()/2.0);

		for( int n = cutoff+1 ; n < vec_in.size()+1 ; ++n) 
		{
			vec_out(p) = vec_in(n);
			p++;
		}

		for( int n = 1 ; n < cutoff+1 ; ++n ) 
		{
			vec_out(p) = vec_in(n);
			p++;
		}

		return vec_out;
	}

    inline cvm::cmatrix fftshift(cvm::cmatrix& freq_sig) 
    {
		cvm::cmatrix freq_sig_shift(freq_sig.msize(), freq_sig.nsize());

		for( int n = 1 ; n <= freq_sig.nsize() ; ++n) 
		{

			cvm::cvector temp_vec(freq_sig(n));

			//cvm::cvector temp_vec(fftshift(temp_vec));
			cvm::cvector freq_sig_vec(fftshift(temp_vec));

			// this is for non-fftshfting
			freq_sig_shift(n) = freq_sig_vec;
		}

		return freq_sig_shift;
	}



    /*! Write a cvm matrix to a file.
     * \param strFilename is the file to write to.
     * \param mA is the matrix to write.
     */
    template<typename matrix_t>
	inline bool writeMatrixToFile(std::string strFilename, matrix_t& mA)
	{
	    std::ofstream fout(strFilename.c_str());

	    if( true == fout.bad() ) {
		std::cerr << "\nerror writing file: " << strFilename;
		return false;
	    }

	    fout << mA.msize();
	    fout << "\n" << mA.nsize();
	    fout << "\n" << mA;

	    return true;
	}

    /*! Read a cvm matrix from a file.
     * \param strFilename is the file in which the matrix is stored.
     * \param mA is the matrix to populate.
     */
    template<typename matrix_t>
	inline bool readMatrixFromFile(std::string strFilename, matrix_t& mA)
	{
	    std::ifstream fin(strFilename.c_str());

	    if( true == fin.bad() ) {
		std::cerr << "\nerror reading file: " << strFilename;
		return false;
	    }

	    integer nRows = -1;
	    integer nCols = -1;
	    fin >> nRows;
	    fin >> nCols;

	    if( (nRows <= 0) or (nCols <= 0) ) {
		std::cerr << "\nerror reading number of rows and cols.";
		return false;
	    }

	    mA.resize(nRows, nCols);

	    fin >> mA;

	    return true;
	}

	inline void diff(const cvm::rvector& vIn, cvm::rvector& vOut)
    {
		vOut.resize( vIn.size() - 1 );
        for ( int n = 1; n < vOut.size() + 1; n++ )
        {
            vOut(n) = -vIn(n)+vIn(n+1);
        }
    }

    // y = mx + c fit
	inline void lsqf(const cvm::rvector& x, const cvm::rvector& y, cvm::rvector& mc)
    {
        double x_sum = 0;
        double x_sq_sum = 0;
        double y_sum = 0;
        double xy_sum = 0;
        int N = y.size();

        for ( int n = 1; n < N + 1; n++ )
        {
            x_sum += x(n);
            x_sq_sum += x(n) * x(n);
            y_sum += y(n); 
            xy_sum += x(n) * y(n); 
        }

        double xy_covar = 1.0/(N - 1.0) * (xy_sum - x_sum * y_sum / N);
        double x_var = 1.0/(N - 1.0) * (x_sq_sum - x_sum * x_sum / N);

        double m = xy_covar / x_var;
        double c = y_sum/N - m * x_sum/N;

        mc.resize(2);
        mc(1) = m;
        mc(2) = c;
    }

	inline void get_fit(const cvm::rvector& x, cvm::rvector& y, const cvm::rvector& mc)
    {
        int N = x.size();
		y.resize( N );

        for ( int n = 1; n < N + 1; n++ )
            y(n) = mc(1) * x(n) + mc(2);
    }
}

#endif
