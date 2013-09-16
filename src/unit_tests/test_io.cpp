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

/*!
 * Attempt to load a DPT file.
 */

BOOST_AUTO_TEST_SUITE( fast )

BOOST_AUTO_TEST_CASE( dpt_spectroscopy_read )
{
	// the testing script will have copied this file to the correct place
	std::string file = "data/dpt_test_input.dpt";

	// set the options
	tarquin::Options options;
	tarquin::Workspace workspace;
	options.SetFormat(tarquin::DANGER);

	tarquin::CBoswell log(tarquin::LOG_NOWHERE);

	tarquin::CFID fid;
	fid.Load(file, options, workspace, log);

	// test that we have read all the parameters correctly
	BOOST_REQUIRE_EQUAL( fid.GetNumberOfPoints(),      1024          );
	BOOST_REQUIRE_CLOSE( fid.GetSamplingFrequency(),   1000,      0.1);
	BOOST_REQUIRE_CLOSE( fid.GetTransmitterFrequency(),6.3684e+07,0.1);
	// TODO: check more parameters here

	// TODO: compare each element of the FID against a binary file that we should construct with the correct
	// data
}

BOOST_AUTO_TEST_CASE( philips_spar_sdat_spectroscopy_read )
{
	std::vector<std::string> test_files;
	test_files.push_back( "data/philips_spar_sdat_WS.SDAT" );
	test_files.push_back( "data/philips_spar_sdat_W.SDAT" );

	tarquin::CBoswell log(tarquin::LOG_NOWHERE);

	for( size_t i = 0; i < test_files.size(); ++i )
	{
		// the file that we are going to work with
		fs::path file = test_files[i];

		// what is the name of the reference data associated with this file?
		fs::path ref_file_base = file.stem();

		// the testing script will have copied this file to the correct place
		std::string ref_file = (file.parent_path()/ref_file_base).string() + ".ref";

		// set the options
		tarquin::Options options;
	    tarquin::Workspace workspace;
		options.SetFormat(tarquin::PHILIPS);

		// attempt to load the FID
		tarquin::CFID fid;
		fid.Load(file.string(), options, workspace, log);

		// test that we have read all the parameters correctly
		BOOST_REQUIRE_EQUAL( fid.GetNumberOfPoints(),      1024          );
		BOOST_REQUIRE_CLOSE( fid.GetSamplingFrequency(),   2000,      0.1);
		BOOST_REQUIRE_CLOSE( fid.GetTransmitterFrequency(),127.786142e6, 0.1);

		// load in the reference data
		cvm::cvector ref = read_ref_data(ref_file);
		for( int i = 0; i < ref.size(); ++i )
		{
			BOOST_REQUIRE_EQUAL( ref[i+1], fid.GetVectorFID(0)[i+1] );
		}
	}
}

BOOST_AUTO_TEST_CASE( siemens_ima_spectroscopy_read )
{
	std::vector<std::string> test_files;
	test_files.push_back( "data/siemens_ima_WS.ima" );

	tarquin::CBoswell log(tarquin::LOG_NOWHERE);

	for( size_t i = 0; i < test_files.size(); ++i )
	{
		// the file that we are going to work with
		fs::path file = test_files[i];

		// what is the name of the reference data associated with this file?
		fs::path ref_file_base = file.stem();

		// the testing script will have copied this file to the correct place
		std::string ref_file = (file.parent_path()/ref_file_base).string() + ".ref";

		// set the options
		tarquin::Options options;
	    tarquin::Workspace workspace;
		options.SetFormat(tarquin::SIEMENS);

		// attempt to load the FID
		tarquin::CFID fid;
		fid.Load(file.string(), options, workspace, log);

		// test that we have read all the parameters correctly
		BOOST_REQUIRE_EQUAL( fid.GetNumberOfPoints(),      2048          );
		BOOST_REQUIRE_CLOSE( fid.GetSamplingFrequency(),   2000,      0.1);
		BOOST_REQUIRE_CLOSE( fid.GetTransmitterFrequency(),63683831, 0.1);

		// load in the reference data
		cvm::cvector ref = read_ref_data(ref_file);
		for( int i = 0; i < ref.size(); ++i )
		{
			BOOST_REQUIRE_EQUAL( ref[i+1], fid.GetVectorFID(0)[i+1] );
		}
	}
}

BOOST_AUTO_TEST_CASE( siemens_rda_spectroscopy_read )
{
	std::vector<std::string> test_files;
	test_files.push_back( "data/siemens_rda_WS.rda" );

	tarquin::CBoswell log(tarquin::LOG_NOWHERE);

	for( size_t i = 0; i < test_files.size(); ++i )
	{
		// the file that we are going to work with
		fs::path file = test_files[i];

		// what is the name of the reference data associated with this file?
		fs::path ref_file_base = file.stem();

		// the testing script will have copied this file to the correct place
		std::string ref_file = (file.parent_path()/ref_file_base).string() + ".ref";

		// set the options
		tarquin::Options options;
	    tarquin::Workspace workspace;
		options.SetFormat(tarquin::RDA);

		// attempt to load the FID
		tarquin::CFID fid;
		fid.Load(file.string(), options, workspace, log);

		// test that we have read all the parameters correctly
		BOOST_REQUIRE_EQUAL( fid.GetNumberOfPoints(),      2048          );
		BOOST_REQUIRE_CLOSE( fid.GetSamplingFrequency(),   2000,      0.1);
		BOOST_REQUIRE_CLOSE( fid.GetTransmitterFrequency(),63683831, 0.1);

		// load in the reference data
		cvm::cvector ref = read_ref_data(ref_file);
		for( int i = 0; i < ref.size(); ++i )
		{
			BOOST_REQUIRE_EQUAL( ref[i+1], fid.GetVectorFID(0)[i+1] );
		}
	}
}

BOOST_AUTO_TEST_SUITE_END()
