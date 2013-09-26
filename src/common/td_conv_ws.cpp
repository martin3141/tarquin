#include "td_conv_ws.hpp"
#include "cvm_util.hpp"
#include <complex> 
#include <cmath>

namespace tarquin {

void fft_conv(cvm::rvector f, cvm::rvector g, cvm::rvector& result)
    {
        // do some checks on f and g first
        if ( f.size() != g.size() )
            std::cout << "size mismatch" << std::endl;

        int N = f.size();
        result.resize(N*2);
        result.set(0);

        // make f complex
        f.resize(N*2);
        cvm::cvector f_comp(N*2);
        f_comp.set(0);
        f_comp.real() = f;

        // find the fft
        cvm::cvector F(N*2);
        fft(f_comp, F);

        // make g complex
        g.resize(N*2);
        cvm::cvector g_comp(N*2);
        g_comp.set(0);
        g_comp.real() = g;

        // find the fft
        cvm::cvector G(N*2);
        fft(g_comp, G);
        
        cvm::cvector result_ifft(N*2);
        cvm::cvector result_comp(N*2);

        G = G.conj();

        // apply the convolution
        for ( int n = 1; n < 2*N + 1; n++ )
            result_ifft(n) = F(n) * G(n);

        ifft(result_ifft, result_comp);
        
        result = fftshift(result_comp.real());
    }


void fft_conv_cyc(const cvm::rvector& f, const cvm::rvector& g, cvm::rvector& result)
    {
        // do some checks on f and g first
        if ( f.size() != g.size() )
            std::cout << "size mismatch" << std::endl;

        int N = f.size();
        result.resize(N);
        result.set(0);

        // make f complex
        cvm::cvector f_comp(N);
        f_comp.set(0);
        f_comp.real() = f;

        // find the fft
        cvm::cvector F(N);
        fft(f_comp, F);

        // make g complex
        cvm::cvector g_comp(N);
        g_comp.set(0);
        g_comp.real() = g;

        // find the fft
        cvm::cvector G(N);
        fft(g_comp, G);
        
        cvm::cvector result_ifft(N);
        cvm::cvector result_comp(N);

        // apply the convolution
        for ( int n = 1; n < N + 1; n++ )
            result_ifft(n) = F(n) * G(n);

        ifft(result_ifft, result_comp);
        
        result = result_comp.real();
    }
 
    void td_conv_ws_noext(const cvm::cvector& S, cvm::cvector& L, const int K) 
	{
		// initialise window function vector and set elements to equal zero
		cvm::rvector wind_fun( 2 * K + 1, 0 );

		// generate gaussian window function
		for ( int k = -K; k < K + 1; k++ ) {
			wind_fun(k+K+1) = exp(-4.0 * std::pow(treal(k),2.0) / std::pow(treal(K),2.0));
			//wind_fun(k+K+1) = cos(k*M_PI / (2.0*K+2.0));
		}

		// normalise window function	
		double wind_fun_sum = 0;
		for ( int n = 1; n < 2 * K + 2; n++ )
			wind_fun_sum += wind_fun(n);
		wind_fun = wind_fun/wind_fun_sum;
		
		// resize water estimate vector	and set elements to equal zero
		integer N = S.size();
		L.resize( N );
		L.set( 0 );
		
		// apply the convolution
		std::complex <double> prod;
		double prod_sum;
		std::complex <double> zero(0,0);
		for ( int n = 1; n < N + 1; n++ ) {
			prod = zero;
			prod_sum = 0;
			for (int m = n - K; m < n + K + 1; m++) {
                if (( m > 0 ) && ( m < N + 1 ))
                {
                    prod += wind_fun( m - (n - K) + 1 ) * S( m );
                    prod_sum += wind_fun( m - (n - K) + 1 );
                }
			}
			L(n) = prod/prod_sum;
		}

		}
    
  void td_conv_ws_noext(const cvm::rvector& S, cvm::rvector& L, const int K) 
	{
		// initialise window function vector and set elements to equal zero
		cvm::rvector wind_fun( 2 * K + 1, 0 );

		// generate gaussian window function
		for ( int k = -K; k < K + 1; k++ ) {
			wind_fun(k+K+1) = exp(-4.0 * std::pow(treal(k),2.0) / std::pow(treal(K),2.0));
			//wind_fun(k+K+1) = cos(k*M_PI / (2.0*K+2.0));
		}

		// normalise window function	
		double wind_fun_sum = 0;
		for ( int n = 1; n < 2 * K + 2; n++ )
			wind_fun_sum += wind_fun(n);
		wind_fun = wind_fun/wind_fun_sum;
		
		// resize water estimate vector	and set elements to equal zero
		integer N = S.size();
		L.resize( N );
		L.set( 0 );
		
		// apply the convolution
		double prod;
		double prod_sum;
		double zero(0);
		for ( int n = 1; n < N + 1; n++ ) {
			prod = zero;
			prod_sum = 0;
			for (int m = n - K; m < n + K + 1; m++) {
                if (( m > 0 ) && ( m < N + 1 ))
                {
                    prod += wind_fun( m - (n - K) + 1 ) * S( m );
                    prod_sum += wind_fun( m - (n - K) + 1 );
                }
			}
			L(n) = prod/prod_sum;
		}

		}


void td_conv_ws(const cvm::cvector& S, cvm::cvector& L, const int K, const int M) 
	{
		// initialise window function vector and set elements to equal zero
		cvm::rvector wind_fun( 2 * K + 1, 0 );

		// generate gaussian window function
		for ( int k = -K; k < K + 1; k++ ) {
            // gaussian
			//wind_fun(k+K+1) = exp(-4.0 * std::pow(treal(k),2.0) / std::pow(treal(K),2.0));
            // sin bell
			wind_fun(k+K+1) = cos(k*M_PI / (2.0*K+2.0));
            // hamming
			//wind_fun(k+K+1) = 0.54 + 0.46*cos(k*M_PI / (2.0*K));
		}

        //std::cout << wind_fun.size() << std::cout;

		// normalise window function	
		double wind_fun_sum = 0;
		for ( int n = 1; n < 2 * K + 2; n++ )
			wind_fun_sum += wind_fun(n);
		wind_fun = wind_fun/wind_fun_sum;

        //plot(wind_fun);
		
		// resize water estimate vector	and set elements to equal zero
		integer N = S.size();
		L.resize( N );
		L.set( 0 );
		
		// apply the convolution
		std::complex <double> prod;
		std::complex <double> zero(0,0);
		for ( int n = 1 + K; n < N - K + 1; n++ ) {
			prod = zero;
			for (int m = n - K; m < n + K + 1; m++) {
				prod += wind_fun( m - (n - K) + 1 ) * S( m );
			}
			L(n) = prod;
		}

		// go back to the start and extrapolate first K points
		std::complex <double> kc(0,0);
		std::complex <double> Mc(M,0);
		for ( int k = -K; k < 0; k++ ) {
			kc = k;
			L(1+K+k) = L(K+1) - kc*(L(K+1)-L(K+M+1))/Mc;
		}

		// go to the end and extrapolate last K points
		for ( int k = 1; k < K+1; k++ ) {
			kc = k;
			L(N-K+k) = L(N-K) + kc*(L(N-K)-L(N-K-M))/Mc;
		}
    }

	void td_conv_ws(const cvm::rvector& S, cvm::rvector& L, const int K, const int M) 
	{
		// initialise window function vector and set elements to equal zero
		cvm::rvector wind_fun( 2 * K + 1, 0 );

		// generate gaussian window function
		for ( int k = -K; k < K + 1; k++ ) {
			//wind_fun(k+K+1) = exp(-4 * pow(k,2) / pow(K,2));
			wind_fun(k+K+1) = cos(k*M_PI / (2.0*K+2.0));
		}

		// normalise window function	
		double wind_fun_sum = 0;
		for ( int n = 1; n < 2 * K + 2; n++ )
			wind_fun_sum += wind_fun(n);
		wind_fun = wind_fun/wind_fun_sum;
		
		// resize water estimate vector	and set elements to equal zero
		integer N = S.size();
		L.resize( N );
		L.set( 0 );
		
		// apply the convolution
		double prod;
		double zero(0);
		for ( int n = 1 + K; n < N - K + 1; n++ ) {
			prod = zero;
			for (int m = n - K; m < n + K + 1; m++) {
				prod += wind_fun( m - (n - K) + 1 ) * S( m );
			}
			L(n) = prod;
		}

		// go back to the start and extrapolate first K points
		double kc(0);
		double Mc(M);
		for ( int k = -K; k < 0; k++ ) {
			kc = k;
			L(1+K+k) = L(K+1) - kc*(L(K+1)-L(K+M+1))/Mc;
		}

		// go to the end and extrapolate last K points
		for ( int k = 1; k < K+1; k++ ) {
			kc = k;
			L(N-K+k) = L(N-K) + kc*(L(N-K)-L(N-K-M))/Mc;
		}
	}

	void td_conv_ws_fft(const cvm::rvector& S, cvm::rvector& L, const int K, const int M) 
	{
		// initialise window function vector and set elements to equal zero
		cvm::rvector wind_fun( 2 * K + 1, 0 );

		// generate gaussian window function
		for ( int k = -K; k < K + 1; k++ ) {
			//wind_fun(k+K+1) = exp(-4 * pow(k,2) / pow(K,2));
			wind_fun(k+K+1) = cos(k*M_PI / (2.0*K+2.0));
		}

		// normalise window function	
		double wind_fun_sum = 0;
		for ( int n = 1; n < 2 * K + 2; n++ )
			wind_fun_sum += wind_fun(n);
		wind_fun = wind_fun/wind_fun_sum;
      
        cvm::rvector wind_fun_full(S.size());
        wind_fun_full.set(0);
		
        for ( int n = 1; n < 2 * K + 1 + 1; n++ )
            wind_fun_full( S.size()/2 - K + n ) = wind_fun(n);

        wind_fun_full = fftshift(wind_fun_full);
        
        //plot(wind_fun_full);

		// resize water estimate vector	and set elements to equal zero
		integer N = S.size();
		L.resize( N );
		L.set( 0 );

        //plot(wind_fun_full);
        fft_conv_cyc(S, wind_fun_full, L);
        //plot(wind_fun_full);

        //plot(L);
		
		// go back to the start and extrapolate first K points
		double kc(0);
		double Mc(M);
		for ( int k = -K; k < 0; k++ ) {
			kc = k;
			L(1+K+k) = L(K+1) - kc*(L(K+1)-L(K+M+1))/Mc;
		}

		// go to the end and extrapolate last K points
		for ( int k = 1; k < K+1; k++ ) {
			kc = k;
			L(N-K+k) = L(N-K) + kc*(L(N-K)-L(N-K-M))/Mc;
		}

        //plot(L);
	}


}
