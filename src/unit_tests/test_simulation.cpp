#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>
#include <stdio.h>
#include <fstream>

#include <boost/filesystem.hpp>

#include "common.hpp"
#include "Options.hpp"
#include "CFID.hpp"
#include "CBasis.hpp"
#include "Workspace.hpp"
#include "preprocess.hpp"
#include "proto_buf.hpp"
#include "test_util.hpp"

namespace fs = boost::filesystem;

BOOST_AUTO_TEST_SUITE( simulator )

BOOST_AUTO_TEST_CASE( simulation_test )
{
	// this contains all final and intermediate results
	tarquin::Workspace workspace;

	// this is where all messages go
	tarquin::CBoswell log(tarquin::LOG_NOWHERE);

	// this is the FID we will be preprocessing, then fitting to
	tarquin::CFID& fidraw = workspace.GetFIDRaw();

	// options
	tarquin::Options& options = workspace.GetOptions();
	options.SetFormat(tarquin::PHILIPS);

	// load the water suppressed FID file
	fidraw.Load("data/philips_spar_sdat_WS.SDAT", options, workspace, log);
	fidraw.SetPPMRef(options.GetRef());

	tarquin::CBasis& basis = workspace.GetBasis();
	BOOST_REQUIRE_EQUAL(basis.Simulate("data", fidraw, options, log), true);

    // bit of a hack as the first simulated file is a dpt file
	cvm::cvector s = basis.GetBasisMatrix()(2);
    
    
    // might want to plot the results for testing purposes
	//tarquin::coord raw_spec = {1, 1, 1};
    /*
	tarquin::coord raw_spec(1, 1, 1);
	cvm::cvector S = tarquin::fft(s);
	S = tarquin::fftshift(S);
	cvm::rvector freq_scale = fidraw.GetPPMScale(raw_spec);
	tarquin::plot(freq_scale,S);
    */

	//write_ref_data(s, "../../unit_tests/data/simulation.ref");
	
	// compare with ref data
	// load in the reference data
	cvm::cvector ref = read_ref_data("data/simulation.ref");
	for( int i = 0; i < ref.size(); ++i )
	{
		float x = s[i+1].real();
		float y = s[i+1].imag();
		
		BOOST_REQUIRE_EQUAL( ref[i+1].real(), x );
		BOOST_REQUIRE_EQUAL( ref[i+1].imag(), y );
	}

}

BOOST_AUTO_TEST_CASE( simulation_test_internal )
{
	// this contains all final and intermediate results
	tarquin::Workspace workspace;

	// this is where all messages go
	tarquin::CBoswell log(tarquin::LOG_NOWHERE);

	// this is the FID we will be preprocessing, then fitting to
	tarquin::CFID& fidraw = workspace.GetFIDRaw();

	// options
	tarquin::Options& options = workspace.GetOptions();
	options.SetFormat(tarquin::PHILIPS);
    options.SetIntBasisSet(tarquin::PROTON_BRAIN_NO_PCR);

	// load the water suppressed FID file
	fidraw.Load("data/philips_spar_sdat_WS.SDAT", options, workspace, log);
	fidraw.SetPPMRef(options.GetRef());

	tarquin::CBasis& basis = workspace.GetBasis();
	const tarquin::Options& options_const = workspace.GetOptions();
	BOOST_REQUIRE_EQUAL(basis.Simulate(fidraw, options_const, log), true);

	cvm::cmatrix s = basis.GetBasisMatrix();
    
    const std::vector<std::string> names = basis.GetSignalNames();
    for ( size_t n = 0; n < names.size(); n++)
    {
		BOOST_TEST_MESSAGE("Checking : " << names[n]);

        cvm::cvector simulated = s(n+1);
        //write_ref_data(simulated, "../../unit_tests/data/"+names[n]+".ref");

        // load in the reference data
        cvm::cvector ref = read_ref_data("data/"+names[n]+".ref");
    
        // might want to plot the results for testing purposes
        /*
        cvm::cmatrix plot_mat(simulated.size(),2);
        plot_mat(1) = simulated;
        plot_mat(2) = ref;
        tarquin::coord raw_spec = {1, 1, 1};
        cvm::cmatrix S = tarquin::fft(plot_mat);
        S = tarquin::fftshift(S);
        cvm::rvector freq_scale = fidraw.GetPPMScale(raw_spec);
        tarquin::plot(freq_scale,S);
        */

        // compare with ref data
		double norm_diff = (ref-simulated).norm();

		BOOST_TEST_MESSAGE("norm_diff = " << norm_diff);
		BOOST_CHECK_LT( norm_diff, 0.004 );
    }
}

BOOST_AUTO_TEST_SUITE_END()
