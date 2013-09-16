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

BOOST_AUTO_TEST_SUITE( fast )

BOOST_AUTO_TEST_CASE( preproc )
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

	tarquin::CFID& fidWater = workspace.GetFIDWater();

	// initialise the water FID to have the same parameters as the WS FID
	fidWater.SetParametersFromFID(fidraw);

	// load the water unsuppressed FID file
	fidWater.Load("data/philips_spar_sdat_W.SDAT", options, workspace, log);

	// this is the FID after some preprocessing, the FID we will actually fit to
	tarquin::CFID& fidproc = workspace.GetFIDProc();

	tarquin::Preprocessor preprocess(workspace, log);

	// preprocess the raw fid
	preprocess(); 

    cvm::cvector test = fidproc.GetVectorFID(0);
    //write_ref_data(test, "../../unit_tests/data/proc.ref");
    //tarquin::plot(test);
    
    // next line doesn't write a good file for some reason...
	//fidproc.SaveToBinFile("../../unit_tests/data/proc.ref");

	
	// compare with ref data
	// load in the reference data
	cvm::cvector ref = read_ref_data("data/proc.ref");
    //tarquin::plot(test);

    //cvm::cvector diff = ref - test;
    //tarquin::plot(diff);

	for( int i = 0; i < ref.size(); ++i )
	{
		float x = fidproc.GetVectorFID(0)[i+1].real();
		float y = fidproc.GetVectorFID(0)[i+1].imag();
        
		BOOST_REQUIRE_EQUAL( ref[i+1].real(), x );
		BOOST_REQUIRE_EQUAL( ref[i+1].imag(), y );
	}

}

BOOST_AUTO_TEST_CASE( auto_phase )
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
	
	tarquin::coord proc_spec(1, 1, 1); 

	tarquin::AutoPhaseNew(proc_spec, fidraw, options, log);
    
    cvm::cvector test = fidraw.GetVectorFID(0);
    //write_ref_data(test, "../../unit_tests/data/phase.ref");

    // next line doesn't write a good file for some reason...
	//fidraw.SaveToBinFile("../../unit_tests/data/phase.ref");

	// compare with ref data
	// load in the reference data
	cvm::cvector ref = read_ref_data("data/phase.ref");
	for( int i = 0; i < ref.size(); ++i )
	{
		float x = fidraw.GetVectorFID(0)[i+1].real();
		float y = fidraw.GetVectorFID(0)[i+1].imag();
		
		BOOST_REQUIRE_EQUAL( ref[i+1].real(), x );
		BOOST_REQUIRE_EQUAL( ref[i+1].imag(), y );
	}

	// check new phi0, phi1
	//BOOST_REQUIRE_CLOSE(fidraw.GetPhi0(proc_spec), 2.07345, 0.001);
	BOOST_REQUIRE_CLOSE(fidraw.GetPhi0(proc_spec), 0.251327, 0.001);
	BOOST_REQUIRE_CLOSE(fidraw.GetPhi1(proc_spec), 0.0000, 0.001);
	
	//std::cout << fidraw.GetPhi0() << std::endl;
	//std::cout << fidraw.GetPhi1() << std::endl;

}

BOOST_AUTO_TEST_CASE( auto_ref )
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
	
	// check starting ref is ok	
	BOOST_REQUIRE_CLOSE(fidraw.GetPPMRef(0), 4.65, 0.001);

	tarquin::coord proc_spec(1, 1, 1); 

	// auto reference using the convolution method
	tarquin::AutoReferenceCorr(proc_spec, options, fidraw, true, false, false, log);

    //std::cout << fidraw.GetPPMRef(0) << std::endl;
	// check new ref is ok	
	BOOST_REQUIRE_CLOSE(fidraw.GetPPMRef(0), 4.66719, 0.001);
}

BOOST_AUTO_TEST_CASE( water_ref )
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

	tarquin::CFID& fidWater = workspace.GetFIDWater();

	// initialise the water FID to have the same parameters as the WS FID
	fidWater.SetParametersFromFID(fidraw);

	// load the water unsuppressed FID file
	fidWater.Load("data/philips_spar_sdat_W.SDAT", options, workspace, log);
	
	tarquin::coord proc_spec(1, 1, 1); 

	tarquin::treal amp = ComputeWaterNormalisation(proc_spec, fidWater, log);

	//std::cout << amp << std::endl;

	BOOST_REQUIRE_CLOSE(amp, 0.2817374618, 0.001);
}

BOOST_AUTO_TEST_SUITE_END()
