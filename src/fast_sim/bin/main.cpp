#include <iostream>
#include <fstream>
#include "fast_sim.hpp"
#include "pulse_sequences.hpp"

int main()
{
	// including the following 2 lines makes the eigenval decomp faster
	// no idea why, perhaps a bug in eigen? TODO find out why
    dcm Ham(2,2);
    Eigen::SelfAdjointEigenSolver<dcm> Hsolver(Ham);

	// set the number of spins
	size_t spin_no = 5;

	// set all elements to 0.5
	drv spin_vec = drv::Constant(spin_no, 0.5);
	
    // set up the spin system parameters
    double B0 = 63e6*2;
	drv chem_shift_vec(spin_no);

    // GABA
    /*
	chem_shift_vec(0) = 3.0128;
	chem_shift_vec(1) = 3.0128;
	chem_shift_vec(2) = 1.889;
	chem_shift_vec(3) = 1.889;
	chem_shift_vec(4) = 2.284;
	chem_shift_vec(5) = 2.284;
    */

	// gln
	chem_shift_vec(0) = 3.753;
	chem_shift_vec(1) = 2.129;
	chem_shift_vec(2) = 2.109;
	chem_shift_vec(3) = 2.432;
	chem_shift_vec(4) = 2.454;

	// lac
    /*
	chem_shift_vec(0) = 4.0974;
	chem_shift_vec(1) = 1.3142;
	chem_shift_vec(2) = 1.3142;
	chem_shift_vec(3) = 1.3142;
    */
    
	drm j_coupling_mat(spin_no, spin_no);

    //GABA
    /*
	j_coupling_mat(0,1) = -12.021;
	j_coupling_mat(0,2) = 5.372;
	j_coupling_mat(0,3) = 7.127;
	j_coupling_mat(1,2) = 10.578;
	j_coupling_mat(1,3) = 6.982;
	j_coupling_mat(2,3) = -13.121;
	j_coupling_mat(2,4) = 7.755;
	j_coupling_mat(2,5) = 7.432;
	j_coupling_mat(3,4) = 6.173;
	j_coupling_mat(3,5) = 7.933;
	j_coupling_mat(4,5) = -10.744;
    */


	// lac
    /*
	j_coupling_mat(0,1) = 6.933;
	j_coupling_mat(0,2) = 6.933;
	j_coupling_mat(0,3) = 6.933;
    */
	
    // gln
	j_coupling_mat(0,1) = 5.847;
	j_coupling_mat(0,2) = 6.5;
	j_coupling_mat(1,2) = -14.504;
	j_coupling_mat(1,3) = 9.165;
	j_coupling_mat(1,4) = 6.347;
	j_coupling_mat(2,3) = 6.324;
	j_coupling_mat(2,4) = 9.209;
	j_coupling_mat(3,4) = -15.371;
	
	//drv group_vec = drv::Constant(spin_no, 2);
	drv group_vec = drv::Constant(spin_no, 1);
	drv spin_num_vec = drv::Constant(spin_no, 1);
	//group_vec(0) = 1;

    dcm time_sig_mat;
	double fs = 2000;
	size_t N = 2048;
	double ref = 4.7;
	double lambda = 8;
	//double tau = 0.03;
	//double tau = 0.03;
	double TE1 = 0.010;
	double TE2 = 0.026;
	//double TE = tau;
	//double TM = 0.012;
	//int cpmg_pulses = 6;

	//pulse_acquire(spin_vec, chem_shift_vec, j_coupling_mat, group_vec, spin_num_vec, B0, fs, N, ref, lambda, time_sig_mat);
	//spin_echo(spin_vec, chem_shift_vec, j_coupling_mat, group_vec, spin_num_vec, B0, fs, N, ref, lambda, time_sig_mat, tau);
	press(spin_vec, chem_shift_vec, j_coupling_mat, group_vec, spin_num_vec, B0, fs, N, ref, lambda, time_sig_mat, TE1, TE2);
	//mega_press(spin_vec, chem_shift_vec, j_coupling_mat, group_vec, spin_num_vec, B0, fs, N, ref, lambda, time_sig_mat, TE1, TE2);
	//steam(spin_vec, chem_shift_vec, j_coupling_mat, group_vec, spin_num_vec, B0, fs, N, ref, lambda, time_sig_mat, TE, TM);
	//laser(spin_vec, chem_shift_vec, j_coupling_mat, group_vec, spin_num_vec, B0, fs, N, ref, lambda, time_sig_mat, tau);
	//cpmg(spin_vec, chem_shift_vec, j_coupling_mat, group_vec, spin_num_vec, B0, fs, N, ref, lambda, time_sig_mat, tau, cpmg_pulses);

	// write FID to file	
	//std::ofstream FID("file_out.csv", std::ios::out);
    //FID << time_sig_mat.col(0) << std::endl;

    std::string strFilename = "3T_asym.dpt";
    std::ofstream fout(strFilename.c_str(), std::ios::out);
	// setup the stream so that things are formatted properly
	fout.setf(std::ios::scientific | std::ios::left);
	fout.precision(8);

    fout << "Dangerplot_version\t" << "2.0" << std::endl;
	fout << "Number_of_points\t" << N << std::endl;
	fout << "Sampling_frequency\t" << fs << std::endl;
	fout << "Transmitter_frequency\t" << B0 << std::endl;
	fout << "Phi0\t" << 0 << std::endl;
	fout << "Phi1\t" << 0 << std::endl;
	fout << "PPM_reference\t" << ref << std::endl;
	fout << "Echo_time\t" << 0 << std::endl;
	fout << "Real_FID\t" << "Imag_FID\t" << std::endl;

	for(int n = 0; n < N; n++)
		fout << real(time_sig_mat.col(0)(n)) << "\t" << imag(time_sig_mat.col(0)(n)) << std::endl;

	return 0;
}
