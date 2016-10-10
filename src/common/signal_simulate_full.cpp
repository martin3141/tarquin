#include "CSignal.hpp"
#include "fast_sim.hpp"
#include "Options.hpp"
#include "pulse_sequences.hpp"
#include <fstream>
#include <vector>

using namespace std;

namespace tarquin {

	int signal_simulate(std::vector<std::vector<double> >& doubmat, const CFID& fidMatch, cvm::cmatrix& matG, const Options& options, std::string name)
	{
        //std::cout << name << std::endl;
		const treal ref = 4.65;
		const treal tau = fidMatch.GetEchoTime();
		
        const treal delay = options.GetAcqDelay();
		
		// assumed for the purposes of the PRESS sequence		
        const treal TE1 = options.GetPRESS_TE1();
		const treal TE2 = tau - TE1;
        const treal TE = tau;
		const treal TM = options.GetSTEAM_TM();

        // TODO add these as options...
		const int cpmg_pulses = 10;
        const treal t1 = TE/4.0;
        const treal t2 = TE/4.0;
        const treal t3 = TE/4.0;
        const treal t4 = TE/4.0;
		
		const treal fs = fidMatch.GetSamplingFrequency();
		const treal B0 = fidMatch.GetTransmitterFrequency();	
		const integer N = fidMatch.GetNumberOfPoints();
		const treal lambda = 0.0;
		
		// find number of spins
		int tot_spin_num = doubmat.size();

		// populate doubvec
		std::vector<double> doubvec;
		for(int n = 0; n < tot_spin_num; n++)
			doubvec.push_back(doubmat[n][0]);
		
		// find the minium of doubvec and put into spinID_start
		std::vector<double>::const_iterator iter = min_element(doubvec.begin(), doubvec.end());
		double spinID_start = *iter;

		// spin ID vector
		drv group_vec(tot_spin_num);
		for(int n = 0; n < tot_spin_num; n++)
			group_vec(n) = doubmat[n][0]-spinID_start+1;

		// read spin_num vector
		dcv spin_num_vec(tot_spin_num);
		for(int n = 0; n < tot_spin_num; n++)
			spin_num_vec(n) = doubmat[n][1];
		
		// read quantum_num vector
		drv spin_vec(tot_spin_num);
		for(int n = 0; n < tot_spin_num; n++)
			spin_vec(n) = doubmat[n][2];

		// read chem shift vector
		drv chem_shift_vec(tot_spin_num);
		for(int n = 0; n < tot_spin_num; n++)
			chem_shift_vec(n) = doubmat[n][3];
		
		// read the J-coupling matrix	
		drm j_coupling_mat(tot_spin_num,tot_spin_num);
		j_coupling_mat.setZero(tot_spin_num,tot_spin_num);
		
		for(int n = 0; n < tot_spin_num; n++)
			for(int m = 0; m < n; m++)
				j_coupling_mat(m,n) = doubmat[n][4+m];

        //std::cout << j_coupling_mat << std::endl;

		dcm time_sig_mat;
        
        pul_seq_e pul_seq = options.GetPulSeq();
		if ( pul_seq == PRESS )
			press(spin_vec, chem_shift_vec, j_coupling_mat, group_vec, spin_num_vec, B0, fs, N, ref, lambda, time_sig_mat, TE1, TE2);
        
        else if ( pul_seq == STEAM )
			steam(spin_vec, chem_shift_vec, j_coupling_mat, group_vec, spin_num_vec, B0, fs, N, ref, lambda, time_sig_mat, TE, TM);

		else if ( pul_seq == LASER )
			laser(spin_vec, chem_shift_vec, j_coupling_mat, group_vec, spin_num_vec, B0, fs, N, ref, lambda, time_sig_mat, tau);
    
        else if ( pul_seq == SEMI_LASER )
			semi_laser(spin_vec, chem_shift_vec, j_coupling_mat, group_vec, spin_num_vec, B0, fs, N, ref, lambda, time_sig_mat, t1, t2, t3, t4);
        
        else if ( pul_seq == PULSE_ACQUIRE )
			pulse_acquire(spin_vec, chem_shift_vec, j_coupling_mat, group_vec, spin_num_vec, B0, fs, N, ref, lambda, time_sig_mat, delay);
        
        else if ( pul_seq == CPMG )
			cpmg(spin_vec, chem_shift_vec, j_coupling_mat, group_vec, spin_num_vec, B0, fs, N, ref, lambda, time_sig_mat, tau, cpmg_pulses);

		else if ( pul_seq == SPIN_ECHO )
			spin_echo(spin_vec, chem_shift_vec, j_coupling_mat, group_vec, spin_num_vec, B0, fs, N, ref, lambda, time_sig_mat, tau);

        else if ( pul_seq == MEGA_PRESS )
			press(spin_vec, chem_shift_vec, j_coupling_mat, group_vec, spin_num_vec, B0, fs, N, ref, lambda, time_sig_mat, TE1, TE2);

        else if ( pul_seq == SHAPED_PRESS )
			shaped_press(spin_vec, chem_shift_vec, j_coupling_mat, group_vec, spin_num_vec, B0, fs, N, ref, lambda, time_sig_mat, TE1, TE2);

        else if ( pul_seq == PROFILE_PRESS )
			profile_press(spin_vec, chem_shift_vec, j_coupling_mat, group_vec, spin_num_vec, B0, fs, N, ref, lambda, time_sig_mat, TE1, TE2);

		else 
		{
			std::cout << "Error pulse sequence not recognised" << std::endl;
			return 1;
		}
        
        //std::cout << name << std::endl << std::flush;
        //std::cout << time_sig_mat(0,0) << std::endl << std::flush;

		// convert time_sig_mat to matG
		matG.resize(time_sig_mat.rows(), time_sig_mat.cols());
		for (int row = 0; row < time_sig_mat.rows(); row++)
			for (int col = 0; col < time_sig_mat.cols(); col++)
				matG(row+1, col+1) = time_sig_mat(row, col);

		return 0;
	}

    int signal_simulate_full(std::vector<std::vector<double> >& full_doubmat, CSignal& signal, const CFID& fidMatch, const Options& options, std::string name)
    {
    //std::cout << name << " simlation started." << std::endl;

	// begin to find chunks of doubmat_full to append to doubmat_vec
	std::vector<std::vector<std::vector<double> > > doubmat_vec;

	// find dimensions of csv sheet
	int rows = full_doubmat.size();

	// int cols = full_doubmat[0].size();
	std::vector<int> chunk_gap;
	for(int r = 0; (r < rows); r++)
	    // make a list of non numerical rows in the csv sheet
	    if (isnan(full_doubmat[r][0]))
				chunk_gap.push_back(r);

	// add one for good luck
	chunk_gap.push_back(rows);

	// populate doubmat_vec
	for(size_t r = 0; (r<(chunk_gap.size()-1)); r++) {
	    std::vector<std::vector<double> > chunk_doubmat;
	    for(int rchunk = chunk_gap[r]+1; (rchunk < chunk_gap[r+1]); rchunk++) { 
		chunk_doubmat.push_back(full_doubmat[rchunk]);
	    }
	    doubmat_vec.push_back(chunk_doubmat);
	}

	// outputs the chunks, keep code for debugging...
	/*
	   for(int z = 0; z<doubmat_vec.size(); z++) {
	   for(int r = 0; r<doubmat_vec[z].size(); r++) {
	   for(int c = 0; c<doubmat_vec[z][r].size(); c++) {
	   cout << doubmat_vec[z][r][c] << " ";
	   }
	   cout << endl;
	   }
	   cout << endl;
	   }
	   */

	int groups = doubmat_vec[0].size();
	const integer N = fidMatch.GetNumberOfPoints();

	// number of groups per chunk
	std::vector<int> chunk_groups;

	// find number of groups per chunk...
	for(size_t chunk_no = 1; chunk_no < doubmat_vec.size(); chunk_no++) {

	    // populate doubvec
	    std::vector<double> doubvec;
	    for(size_t n = 0; n < doubmat_vec[chunk_no].size(); n++)
		doubvec.push_back(doubmat_vec[chunk_no][n][0]);

	    // find the minium of doubvec and put into spinID_start
	    std::vector<double>::const_iterator itermin = min_element(doubvec.begin(), doubvec.end());
	    int spinID_start = (int) *itermin;

	    // find the maximum of doubvec and put into spinID_start
	    std::vector<double>::const_iterator itermax = max_element(doubvec.begin(), doubvec.end());
	    int spinID_end = (int) *itermax;
	    chunk_groups.push_back(spinID_end-spinID_start+1);
	}

	// cumlative group chunk number
	int cum_group_chunk_no = 1;

	// all the groups are stored here until summation
	cvm::cmatrix matG(N, groups);

	// simulate each group
	for(size_t chunk_no = 1; chunk_no < doubmat_vec.size(); chunk_no++) {

	    cvm::cmatrix chunkG;	
	    signal_simulate(doubmat_vec[chunk_no], fidMatch, chunkG, options, name); 

        //std::cout << chunkG;
        //std::cout << "Simulating spin system " << chunk_no << " of " << doubmat_vec.size() - 1 << std::endl;

	    matG.assign(1, cum_group_chunk_no,chunkG);
	    cum_group_chunk_no = cum_group_chunk_no + chunk_groups[chunk_no-1];
	}

	// read in damping values
	vector<treal> alpha;
	vector<treal> beta;
    treal max_alpha = -std::numeric_limits<treal>::infinity();
    treal max_beta = -std::numeric_limits<treal>::infinity();

	for( size_t y = 0; y < doubmat_vec[0].size(); y++ ) {
        treal tmp_alpha = doubmat_vec[0][y][4];
        treal tmp_beta = doubmat_vec[0][y][5];
        alpha.push_back(tmp_alpha);
        beta.push_back(tmp_beta);
        if ( tmp_alpha > max_alpha )
            max_alpha = tmp_alpha;
        if ( tmp_beta > max_beta )
            max_beta = tmp_beta;
	}
    
    /*if ( max_alpha > 10 || max_beta > 10 ) 
    {
        //std::cout << "True" << std::endl;
        broad_sig.push_back(true);
    }
    else
    {
        //std::cout << "False" << std::endl;
        broad_sig.push_back(false);
    }*/

	// apply some damping
	const treal fs = fidMatch.GetSamplingFrequency();
	cvm::rvector time_scale(N);
	for(int n = 0; n < N; n++) 
	    time_scale(n+1) = n*(1.0/fs); 



	for(int r = 0; r < matG.msize(); r++) {

	    //tcomplex z = 0;

	    for(int c = 0; c < matG.nsize(); c++) {
		matG(r+1,c+1) = matG(r+1,c+1)*exp(-alpha[c]*time_scale(r+1)-beta[c]*pow(time_scale(r+1),2));

		//z = z + matG(r+1, c+1);
	    }
	    //signal.m_fids[0].appendSample(z);	
	}

	// read in group values
	vector<integer> group;
	vector<integer> independance;
	for( size_t y = 0; y < doubmat_vec[0].size(); y++ ) {
	    group.push_back( (int) doubmat_vec[0][y][0]);
	    independance.push_back( (int) doubmat_vec[0][y][1]);
	}

	// find number of independant groups	
	size_t max = 0;
	for( size_t y = 0; y < independance.size(); y++ ) {
	    if (independance[y] > (int) max)
		max = independance[y];
	}

	// resize m_fids 
	signal.m_fids.resize(max);

	// define the matrix of signals
	cvm::cmatrix matG_indep(N, max);

	matG_indep.set(0);
	for( size_t y = 1; y <= independance.size(); y++ ) {
	    for( size_t m = 1; m <= max; m++ ) {
		if (independance[y-1] == (int) m) {
		    for( size_t n = 1; (int) n <= N ; n++ ) {
			matG_indep(n,m) += matG(n,y);
		    }
		}
	    }
	}

	// transfer simulation results to FID object
	for( int c = 0; c < matG_indep.nsize(); c++ ) {
	    signal.m_fids[c].AppendFromVector(matG_indep(c+1));

		//const integer N = fidMatch.GetNumberOfPoints();
	    signal.m_fids[c].SetSamplingFrequency(fidMatch.GetSamplingFrequency());
	    signal.m_fids[c].SetTransmitterFrequency(fidMatch.GetTransmitterFrequency());
	    // signal.m_fids[c].SetPPMRef(fidMatch.GetPPMRef());
	    // hard coded
		pair_vec default_ref;
		default_ref.push_back(std::make_pair(4.65, true));
	    signal.m_fids[c].SetPPMRef(default_ref);
	}

	// do fft
	//cvm::cmatrix freq_sig(fft(matG));
	//cvm::cmatrix freq_sig(fft(matG_indep));
	// fftshift	
	//freq_sig = fftshift(freq_sig);
	// get PPM scale
	// cvm::rvector freq_scale = fidMatch.getPPMScale();
	// and plot
	//plot(freq_sig);

    //std::cout << name << " simlation completed." << std::endl;
	return 0;
    }
}
