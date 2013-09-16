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
#include "export_data.hpp"

namespace fs = boost::filesystem;

BOOST_AUTO_TEST_CASE( nb_fit_test_no_iters )
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
	options.SetMaxIters(0);

	// load the water suppressed FID file
	fidraw.Load("data/philips_spar_sdat_WS.SDAT", options, workspace, log);
	fidraw.SetPPMRef(options.GetRef());

	tarquin::CFID& fidWater = workspace.GetFIDWater();
	tarquin::CFID& fidproc = workspace.GetFIDProc();

	// initialise the water FID to have the same parameters as the WS FID
	fidWater.SetParametersFromFID(fidraw);

	// load the water unsuppressed FID file
	fidWater.Load("data/philips_spar_sdat_W.SDAT", options, workspace, log);

    // give it a name so that water scaling is performed
    options.SetFilenameWater("dummy");

	tarquin::Preprocessor preprocess(workspace, log);

	// preprocess the raw fid
	preprocess(); 
	tarquin::CBasis& basis = workspace.GetBasis();

    // simuate basis "on the fly"
	const tarquin::Options& options_const = workspace.GetOptions();
    basis.Simulate(fidproc, options_const, log);
    
    // pre-compiled basis
	//BOOST_REQUIRE_EQUAL(load_basis("data/test_basis.xml", basis), true);
	//BOOST_REQUIRE_EQUAL(basis.check(fidraw, log), true);

	// do the analysis
	BOOST_REQUIRE_EQUAL(RunTARQUIN(workspace, log), true);

	const tarquin::rvec_stdvec& ahat_vec = workspace.GetAmplitudesNormalised();
	cvm::rvector ahat = ahat_vec[0];
	
    
    // version 4.2.7
     //double ahat_array[] = {0.123933, 0, 4.31036, 0, 1.65386, 1.49086, 0.508613, 4.07789, 1.00659, 0.0746496, 3.63818, 0.170727, 1.32006, 1.71451, 0, 0, 0.575216, 0.195786, 2.70968, 0, 5.66233, 6.3855, 0.25126, 0.0209098, 0.117753, 0};

     // version 4.2.8
    //double ahat_array[] = {0.123896, 0, 4.31036, 0, 1.65386, 1.49085, 0.508613, 4.07789, 1.00658, 0.0746535, 3.63818, 0.170318, 1.32007, 1.71495, 0, 0, 0.575213, 0.195716, 2.70982, 0, 5.66232, 6.3855, 0.25126, 0.0209142, 0.117753, 0};
    
    //std::cout << std::endl << ahat << std::endl;

    // new version 
    //double ahat_array[] = {0.187897, 0.991977, 4.47883, 0, 0.931229, 0, 1.01458, 4.29613, 1.04516, 1.12278, 3.37964, 0.186356, 1.01185, 0.79202, 0, 0, 1.27544, 0, 2.11104, 0.246311, 6.80434, 5.70221, 1.17356, 0, 0.120588, 0};

    // 4.2.11
    double ahat_array[] = {0.298991, 1.28129, 4.52445, 0, 1.13412, 0, 1.24417, 4.32446, 1.06926, 1.23649, 3.58885, 0.255948, 0.996417, 0.882861, 0, 0, 1.37716, 0.0562968, 2.15641, 0.368497, 6.65246, 5.77977, 1.21318, 0, 0.162392, 0};

	// ensure the two are close
	for( int i = 0; i < ahat.size(); ++i )
		BOOST_CHECK_CLOSE( ahat_array[i], ahat(i+1), 0.1 );

}


BOOST_AUTO_TEST_CASE( nb_fit_test )
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
	tarquin::CFID& fidproc = workspace.GetFIDProc();

	// initialise the water FID to have the same parameters as the WS FID
	fidWater.SetParametersFromFID(fidraw);

	// load the water unsuppressed FID file
	fidWater.Load("data/philips_spar_sdat_W.SDAT", options, workspace, log);

    // give it a name so that water scaling is performed
    options.SetFilenameWater("dummy");

	tarquin::Preprocessor preprocess(workspace, log);

	// preprocess the raw fid
	preprocess(); 
	tarquin::CBasis& basis = workspace.GetBasis();

    // simuate basis "on the fly"
    basis.Simulate(fidproc, options, log);
	
	// do the analysis
	BOOST_REQUIRE_EQUAL(RunTARQUIN(workspace, log), true);

	const tarquin::rvec_stdvec& ahat_vec = workspace.GetAmplitudesNormalised();
	cvm::rvector ahat = ahat_vec[0];
	
    //std::cout << std::endl << ahat << std::endl;
    
    // version 4.2.7
    //double ahat_array[] = {0.185446, 2.26678, 5.57108, 0, 1.25154, 0.900754, 0.998319, 4.33083, 0.941118, 1.0519, 3.46246, 0, 1.99855, 1.9584, 0.244112, 1.38244, 1.16133, 0.129494, 1.98681, 0.818801, 7.67897, 6.62238, 0, 0.224791, 0.232873, 0};
    
    // version 4.2.8
    //double ahat_array[] = {0.186736, 2.26485, 5.57102, 0, 1.25259, 0.901531, 0.99713, 4.32977, 0.941039, 1.052, 3.46265, 0, 1.99651, 1.95962, 0.243793, 1.38187, 1.16648, 0.13018, 1.98199, 0.819814, 7.67822, 6.62255, 0, 0.224926, 0.232789, 0};

    // version 4.2.8 with fp correction
    //double ahat_array[] = {0.131951, 1.75534, 4.77611, 0, 2.08089, 1.00926, 1.06783, 4.14403, 0.954601, 0.765174, 3.497, 0.0340521, 0.783153, 2.98286, 0, 0, 1.62311, 0, 3.63711, 1.00182, 8.58571, 6.48766, 0.136921, 0.214334, 0.217733, 0};

    // new version 
    //double ahat_array[] = {0.283047, 2.22998, 4.72108, 0, 0.837213, 0, 1.55278, 5.3586, 1.07697, 1.21919, 3.71015, 0.198904, 1.44094, 1.50791, 0, 0, 1.60879, 0.230414, 2.26488, 1.11687, 7.38569, 6.48249, 0.602102, 0, 0.216289, 0};

    // 4.2.11
    double ahat_array[] = {0.188134, 1.22779, 6.16942, 0, 1.40724, 0.629777, 1.15256, 5.34849, 1.50427, 1.55262, 3.55217, 0.12896, 0.349893, 0, 0, 0, 2.94453, 0.428997, 1.53264, 1.48881, 11.0323, 5.35602, 1.68976, 0, 0.442839, 0};
    

	// ensure the two are close
	for( int i = 0; i < ahat.size(); ++i )
		BOOST_CHECK_CLOSE( ahat_array[i], ahat(i+1), 0.1 );

}
