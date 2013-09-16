#ifndef __TDCONVWS__
#define __TDCONVWS__

#include "CFID.hpp"

namespace tarquin {

	void td_conv_ws(const cvm::cvector& y, cvm::cvector& yw, const int K, const int M);
	
    void td_conv_ws_noext(const cvm::cvector& y, cvm::cvector& yw, const int K);
    
    void td_conv_ws_noext(const cvm::rvector& y, cvm::rvector& yw, const int K);
	
	void td_conv_ws(const cvm::rvector& y, cvm::rvector& yw, const int K, const int M);
	
    void td_conv_ws_fft(const cvm::rvector& y, cvm::rvector& yw, const int K, const int M);
    
    void fft_conv_cyc(const cvm::rvector& f, const cvm::rvector& g, cvm::rvector& result);
    
    void fft_conv(cvm::rvector f,cvm::rvector g, cvm::rvector& result);

}

#endif
