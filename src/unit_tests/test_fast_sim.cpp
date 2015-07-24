#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>
#include <stdio.h>
#include <fstream>
#include <boost/filesystem.hpp>

#include "test_util.hpp"
#include "fast_sim.hpp"
#include "pulse_sequences.hpp"

namespace fs = boost::filesystem;

BOOST_AUTO_TEST_SUITE( fast_sim )

BOOST_AUTO_TEST_CASE( lacate_press )
{
    // set the number of spins
	size_t spin_no = 4;

	// set all elements to 0.5
	drv spin_vec = drv::Constant(spin_no, 0.5);
	
    // set up the spin system parameters
    double B0 = 1.27786e8;
	drv chem_shift_vec(spin_no);

	// lac
	chem_shift_vec(0) = 4.0974;
	chem_shift_vec(1) = 1.3142;
	chem_shift_vec(2) = 1.3142;
	chem_shift_vec(3) = 1.3142;

	drm j_coupling_mat(spin_no, spin_no);
	j_coupling_mat.setZero(spin_no, spin_no);
	
	// lac
	j_coupling_mat(0,1) = 6.933;
	j_coupling_mat(0,2) = 6.933;
	j_coupling_mat(0,3) = 6.933;
	
	drv group_vec = drv::Constant(spin_no, 1);
	dcv spin_num_vec = dcv::Constant(spin_no, 1);

	double fs = 2000;
	size_t N = 2048;
	double ref = 4.7;
	double lambda = 5;
	double TE1 = 0.01;
	double TE2 = 0.02;
    
    // initialise the spin system
	spin_sys sys(spin_vec, chem_shift_vec, j_coupling_mat, B0);
    
    // check H
    dcm H;
    sys.get_H(H);
    H.resize(H.rows()*H.cols(),1);

	//write_ref_data_eigen(H.col(0), "../../unit_tests/data/fast_sim_H.ref");

    dcv H_ref;
    H_ref = read_ref_data_eigen("data/fast_sim_H.ref");
    for( int i = 0; i < H.size(); ++i )
	{
		double x = H.col(0)[i].real();
		double y = H.col(0)[i].imag();
		
		BOOST_CHECK_LT( std::fabs(H_ref[i].real() - x), 1e-4 );
		BOOST_CHECK_LT( std::fabs(H_ref[i].imag() - y), 1e-4 );
	}

    // check eigen vectors
    dcm eig_vecs;
    sys.get_eig_vecs(eig_vecs);

    //std::cout << eig_vecs << std::endl;

    eig_vecs.resize(eig_vecs.rows()*eig_vecs.cols(),1);

	//write_ref_data_eigen(eig_vecs.col(0), "../../unit_tests/data/fast_sim_eig_vecs.ref");
    /* 
    dcv eig_vecs_ref;
    eig_vecs_ref = read_ref_data_eigen("data/fast_sim_eig_vecs.ref");
    for( int i = 0; i < eig_vecs_ref.size(); ++i )
	{
		float x = eig_vecs.col(0)[i].real();
		float y = eig_vecs.col(0)[i].imag();
		
		BOOST_REQUIRE_EQUAL( eig_vecs_ref[i].real(), x );
		BOOST_REQUIRE_EQUAL( eig_vecs_ref[i].imag(), y );
	}
    */

    // set the state of the system to be -Fy
    scm rho;
    sys.Fy(rho, spin_num_vec);
    rho *= -1;
	sys.set_state(rho);
    
    // check initial state of system
    dcm init_state;
    sys.get_state(init_state);
    init_state.resize(init_state.rows()*init_state.cols(),1);
	//write_ref_data_eigen(init_state.col(0), "../../unit_tests/data/fast_sim_init_state.ref");
    
    dcv init_state_ref;
    init_state_ref = read_ref_data_eigen("data/fast_sim_init_state.ref");
    for( int i = 0; i < init_state_ref.size(); ++i )
	{
		float x = init_state.col(0)[i].real();
		float y = init_state.col(0)[i].imag();
		
		BOOST_REQUIRE_EQUAL( init_state_ref[i].real(), x );
		BOOST_REQUIRE_EQUAL( init_state_ref[i].imag(), y );
	}

    dcm delay_a;
	dcm delay_b;
	sys.get_delay_ops(TE1/2.0, delay_a, delay_b);
    
	//write_ref_data_eigen(delay_a.col(0), "../../unit_tests/data/fast_sim_delay_a_state.ref");

    dcv delay_a_state_ref;
    delay_a_state_ref = read_ref_data_eigen("data/fast_sim_delay_a_state.ref");
    for( int i = 0; i < delay_a_state_ref.size(); ++i )
	{
		float x = delay_a.col(0)[i].real();
		float y = delay_a.col(0)[i].imag();
		
		BOOST_REQUIRE_EQUAL( delay_a_state_ref[i].real(), x );
		BOOST_REQUIRE_EQUAL( delay_a_state_ref[i].imag(), y );
	}
    
    //write_ref_data_eigen(delay_b.col(0), "../../unit_tests/data/fast_sim_delay_b_state.ref");

    dcv delay_b_state_ref;
    delay_b_state_ref = read_ref_data_eigen("data/fast_sim_delay_b_state.ref");
    for( int i = 0; i < delay_b_state_ref.size(); ++i )
	{
		float x = delay_b.col(0)[i].real();
		float y = delay_b.col(0)[i].imag();
		
		BOOST_REQUIRE_EQUAL( delay_b_state_ref[i].real(), x );
		BOOST_REQUIRE_EQUAL( delay_b_state_ref[i].imag(), y );
	}

	sys.delay(TE1/2.0);

    dcm post_d1_state;
    sys.get_state(post_d1_state);
    post_d1_state.resize(post_d1_state.rows()*post_d1_state.cols(),1);
	//write_ref_data_eigen(post_d1_state.col(0), "../../unit_tests/data/fast_sim_post_d1_state.ref");
    
    /*
    dcv post_d1_state_ref;
    post_d1_state_ref = read_ref_data_eigen("data/fast_sim_post_d1_state.ref");
    for( int i = 0; i < post_d1_state_ref.size(); ++i )
	{
		float x = post_d1_state.col(0)[i].real();
		float y = post_d1_state.col(0)[i].imag();
		
		BOOST_REQUIRE_EQUAL( post_d1_state_ref[i].real(), x );
		BOOST_REQUIRE_EQUAL( post_d1_state_ref[i].imag(), y );
	}
    */

	dcm a;
	dcm b;
    sys.get_pulse_ops(180, "y", a, b);
	sys.apply_op(a, b);
	sys.delay((TE1+TE2)/2.0);
	sys.apply_op(a, b);
	sys.delay(TE2/2.0);

    dcm pre_acq_state;
    sys.get_state(pre_acq_state);
    pre_acq_state.resize(pre_acq_state.rows()*pre_acq_state.cols(),1);
	//write_ref_data_eigen(pre_acq_state.col(0), "../../unit_tests/data/fast_sim_pre_acq_state.ref");
    
    /*
    dcv pre_acq_state_ref;
    pre_acq_state_ref = read_ref_data_eigen("data/fast_sim_pre_acq_state.ref");
    for( int i = 0; i < pre_acq_state_ref.size(); ++i )
	{
		float x = pre_acq_state.col(0)[i].real();
		float y = pre_acq_state.col(0)[i].imag();
		
		BOOST_REQUIRE_EQUAL( pre_acq_state_ref[i].real(), x );
		BOOST_REQUIRE_EQUAL( pre_acq_state_ref[i].imag(), y );
	}
    */

	// acquire the result
    dcm time_sig_mat;
	sys.acquire(time_sig_mat, fs, N, ref, lambda, group_vec);
    
    //write_ref_data_eigen(time_sig_mat.col(0), "../../unit_tests/data/fast_sim_time_sig.ref");
    
    dcv time_sig_ref;
    time_sig_ref = read_ref_data_eigen("data/fast_sim_time_sig.ref");
    for( int i = 0; i < time_sig_ref.size(); ++i )
	{
		float x = time_sig_mat.col(0)[i].real();
		float y = time_sig_mat.col(0)[i].imag();
		
		BOOST_REQUIRE_EQUAL( time_sig_ref[i].real(), x );
		BOOST_REQUIRE_EQUAL( time_sig_ref[i].imag(), y );
	}
}

/*
BOOST_AUTO_TEST_CASE( lacate_press_full ) // In progress...
{
    // set the number of spins
	size_t spin_no = 4;

	// set all elements to 0.5
	drv spin_vec = drv::Constant(spin_no, 0.5);
	
    // set up the spin system parameters
    double B0 = 1.27786e8;
	drv chem_shift_vec(spin_no);

	// lac
	chem_shift_vec(0) = 4.0974;
	chem_shift_vec(1) = 1.3142;
	chem_shift_vec(2) = 1.3142;
	chem_shift_vec(3) = 1.3142;

	drm j_coupling_mat(spin_no, spin_no);
	
	// lac
	j_coupling_mat(0,1) = 6.933;
	j_coupling_mat(0,2) = 6.933;
	j_coupling_mat(0,3) = 6.933;
	
	drv group_vec = drv::Constant(spin_no, 1);
	drv spin_num_vec = drv::Constant(spin_no, 1);

    dcm time_sig_mat;
	double fs = 2000;
	size_t N = 2048;
	double ref = 4.7;
	double lambda = 5;
	double TE1 = 0.01;
	double TE2 = 0.02;
    
    press(spin_vec, chem_shift_vec, j_coupling_mat, group_vec, spin_num_vec, B0, fs, N, ref, lambda, time_sig_mat, TE1, TE2);
    
    //std::cout << time_sig_mat.col(0);

}
*/

BOOST_AUTO_TEST_SUITE_END()
