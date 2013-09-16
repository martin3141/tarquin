#include <stdlib.h>
#include "../common/CBasis.hpp"
#include "../common/common.hpp"

#include <stdlib.h>
#include <algorithm>
#include <complex>
#include <sstream>
//#include "boost/filesystem.hpp"   // includes all needed Boost.Filesystem declarations

using namespace tarquin;
using namespace std;
namespace fs = boost::filesystem;

int main(int argc, char* argv[])
{
    // name of the directory to store the results
    fs::path out_path( "../../data/simulated/SNR" );
    fs::path out_path_fit( "./" );

    // name of the directory to store the raw data 
    fs::path out_path_raw_data( out_path / "RAW_data" );
    fs::path out_path_raw_data_fit( out_path_fit / "RAW_data" );

    // name of the directory to store the TARQUIN fits 
    fs::path out_path_TARQUIN( out_path / "TARQUIN" );
    fs::path out_path_TARQUIN_fit( out_path_fit / "TARQUIN" );

    // name of the directory to store the LCModel fits 
    fs::path out_path_LCM( out_path / "LCModel" );
    fs::path out_path_LCM_fit( out_path_fit / "LCModel" );
    
    // name of the directory to store the simulated amplitudes 
    fs::path out_path_amps( out_path / "amplitudes" );

    // delete old results path if it exists
    if ( fs::exists( out_path ) ) {
	std::cout << "Removing old simulations.\n";
	fs::remove_all( out_path );
    }

    // create a new directory for the results
    fs::create_directory( out_path );
    
    // create a new directory for the raw data
    fs::create_directory( out_path_raw_data );
    // create a new directory for the TARQUIN fits 
    fs::create_directory( out_path_TARQUIN );
    // create a new directory for the LCModel fits 
    fs::create_directory( out_path_LCM );
    // create a new directory for the simulated ampltudes 
    fs::create_directory( out_path_amps );

    // create the dummy FID object
    CFID fid;    
    fid.SetSamplingFrequency(1000.0);
    fid.SetTransmitterFrequency(6.3684e7);
    fid.SetEchoTime(30e-3);
    fid.SetPPMRef(4.7);
    fid.SetNumberOfPoints(1024);
    std::string strBasis = "../../basis/LCM.xml";
    std::string strBasis_fit = "../../../basis/LCM.xml";

    // create the basis object
    CBasis basis;

    // load a test basis
    //load_xml(basis, strBasis);
    //basis.GenerateNonSTLstruc();

    // define matrix of signals (no groups) S 
    const cvm::cmatrix& Q = basis.GetBasisMatrix();
    cvm::cmatrix S = Q;

    // generate a water signal such that the metabolite amplitudes are scaled to unity
    treal dt = 1.0 / fid.GetSamplingFrequency();
    treal t;
    size_t N = fid.GetNumberOfPoints();
    cvm::cvector waterFID(N);
    for(size_t n = 0; n < N; n++) {
	t = n*dt;
	waterFID(n+1) = 0.7*35880.0*exp(-5*t*t);
    }

    // initialise random number generator
    srand(1);

    // a represents signal amplitudes
    cvm::cvector a(S.nsize());

    // randgen - number of random spectra to generate per paramter permutation
    int randgen = 5;
    
    // randgen - number of variable permutations to run
    int para = 5;

    cvm::cvector y(S.msize());

    for(int p = 1; p <= para ;p++) {
    for(int m = 1; m <= randgen ;m++) {
	for(int n = 1; n < S.nsize()+1; n++) {	
	    // mean 10 stdev 5 (amplitude of basis column)
	    //a(n) = GeneratePosGauss(10,15);
	    a(n) = 1;
	}
	//a(S.nsize()) = 0;

	// convert m to string
	stringstream ss;
	string FileNum;

	ss.fill('0');
	ss.width(3);
	ss << p;
	ss << "-";
	ss.width(3);
	ss << m;
	ss >> FileNum;

	// add components together
	y = (S*a);

	// add noise
	//y += GenerateNoise(y.size(),0.2*(p-1));

	for(int n = 1; n <= y.size(); ++n ) {
	    t = (n-1)*dt;
	    y(n) = y(n) * exp(-5*t*t);
	}

	// distort the first few points
	complex <double> dista(8,8);
	complex <double> distb(4,4);
	//y(1) = y(1)*dista;
	//y(2) = y(2)*distb;
	//y(3) = y(3)*(complex <double>)1.5;

	// apply additional damping
	/*for(size_t n = 0; n < N; n++) {
	    t = n*dt;
	    y(n+1) = y(n+1)*exp(-50*t*t);
	}*/

	// path to simulated water file
	string W_file = FileNum + "_W_RAW";
	fs::path W_path( out_path_raw_data / W_file );
	fs::path W_path_fit( out_path_raw_data_fit / W_file );

	// path to simulated water supressed file
	string WS_file = FileNum + "_WS_RAW";
	fs::path WS_path( out_path_raw_data / WS_file );
	fs::path WS_path_fit( out_path_raw_data_fit / WS_file );

	// path to TARQUIN output image
	string eps_file = FileNum + "_plot.eps";
	fs::path image_path_eps( out_path_TARQUIN_fit / eps_file );
	//fs::path image_path_fit( out_path_TARQUIN_fit / eps_file );

	// path to LCModel output image
	string ps_file = FileNum + "_plot.ps";
	fs::path image_path_ps( out_path_LCM_fit / ps_file );
	//fs::path image_path_fit_LCM( out_path_TARQUIN_fit / eps_file );

	string csv_file = FileNum + "_results.csv";
	// path to output csv file 
	// fs::path csv_path( out_path_TARQUIN / csv_file );
	fs::path csv_path( out_path_TARQUIN_fit / csv_file );
	fs::path csv_path_LCM( out_path_LCM_fit / csv_file );

	// now save y as LCM format
	CFID WS_out = fid;
	WS_out.InitialiseFromVector(y);
	WS_out.SaveToFileLCM(WS_path.native_file_string());

	// now save y as LCM format
	CFID W_out = fid;
	W_out.InitialiseFromVector(waterFID);
	W_out.SaveToFileLCM(W_path.native_file_string());

	// generate the command to run tarquin on the simulated data
	// std::string run_tarquin = "./tarquin --basis_xml "+strBasis+" --input "+WS_path.native_file_string()+" --input_w "+W_path.native_file_string()+" --format lcm --echo 0.03 --fs 1000 --ft 6.3684e7 --auto_phase false --auto_ref false --output_image "+image_path.native_file_string()+" --output_csv "+csv_path.native_file_string();

	//std::string run_tarquin = "../../../src/bin/tarquin --basis_xml "+strBasis+" --input "+WS_path.native_file_string()+" --input_w "+W_path.native_file_string()+" --format lcm --echo 0.03 --fs 1000 --ft 6.3684e7 --output_image "+image_path.native_file_string()+" --output_csv "+csv_path.native_file_string()+" --water_eddy true";
	

	std::string run_tarquin = "../../../src/bin/tarquin --basis_xml "+strBasis_fit+" --input "+WS_path_fit.native_file_string()+" --input_w "+W_path_fit.native_file_string()+" --format lcm --echo 0.03 --fs 1000 --ft 6.3684e7 --output_image "+image_path_eps.native_file_string()+" --output_csv "+csv_path.native_file_string()+" --water_eddy true";

	// generate tarquin batch job
	fs::path tarquin_batch( out_path / "goTARQUIN" );
	std::string file_name = tarquin_batch.native_file_string();
	ofstream file( file_name.c_str(), ios::app );
	file << run_tarquin << std::endl;
	file.close();

	std::string LCMpath = "~/.lcmodel/bin/lcmodel";
	// generate LCModel batch job
	fs::path LCM_batch( out_path / "goLCModel" );
	file_name = LCM_batch.native_file_string();
	ofstream LCM_file( file_name.c_str(), ios::app );
	LCM_file << "echo \"";
	LCM_file << " \\$LCMODL";
	LCM_file << " OWNER='Institute of Child Health, University of Birmingham'";
	LCM_file << " KEY(1)=599209454";
	LCM_file << " FILBAS='/export/spare/acp/.lcmodel/basis-sets/Simulated_1p5T.basis'";
	LCM_file << " FILRAW='" << WS_path_fit.native_file_string() << "'";
	LCM_file << " FILH2O='" << W_path_fit.native_file_string() << "'";
	LCM_file << " DOWS=T";
	LCM_file << " DOECC=T";
	LCM_file << " FILPS='" << image_path_ps.native_file_string() << "'";
	LCM_file << " LCSV=11";
	LCM_file << " FILCSV='" << csv_path_LCM.native_file_string() << "'";
	LCM_file << " SPTYPE='tumor'";
	LCM_file << " HZPPPM=6.3684e+01";
	LCM_file << " DELTAT=1.000e-03";
	LCM_file << " NUNFIL=1024";
	LCM_file << " NRATIO=0";
	LCM_file << " NSIMUL=11";
	LCM_file << " \\$END\" | ";
	LCM_file << LCMpath << std::endl;
	LCM_file.close();
	
	// run the command
	//system(run_tarquin.c_str());


	// output the simulated amplitudes to a csv file
	// path to simulated amplitudes csv file 
	string sim_amp_file = FileNum + "_amp.csv";
	fs::path sim_amp_path( out_path_amps / sim_amp_file );
	std::string strFilename = sim_amp_path.native_file_string();
	std::ofstream fout(strFilename.c_str(), std::ios::out);
	fout << std::showpoint;
	fout << std::scientific;

	std::vector < std::string > signal_names = basis.GetSignalNames();	
	std::vector < std::string > signal_names_sorted = signal_names;

	std::sort(signal_names_sorted.begin(), signal_names_sorted.end()); 

	double tempamp;
	for(integer n = 1; n < a.size(); n++) {
	    fout << signal_names_sorted[n-1] << ",";
	}
	// last one without the comma	
	fout << signal_names_sorted[a.size()-1] << std::endl;


	for(integer n = 1; n < a.size(); n++) {
	    // find index of a now the list is sorted
	    for(integer m = 1; m < a.size()+1; m++) {
		if ( signal_names_sorted[n-1] == signal_names[m-1] ){
		    tempamp = a(m).real();
		    break;
		}
	    }
	    fout << tempamp << ",";	
	}

	// last one without the comma	
	integer n = a.size();
	for(integer m = 1; m < a.size()+1; m++) {
	    if ( signal_names_sorted[n-1] == signal_names[m-1] ){
		tempamp = a(m).real();
		break;
	    }
	}
	fout << tempamp;
    }
    }

    return 0;
}
