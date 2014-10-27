#include "CBasis.hpp"
#include "CCSVFile.hpp"
#include "signal_simulate_full.hpp"
#include "basis_set.hpp"
#include "Options.hpp"
#include "fidio/CFIDReaderDPT.hpp"

#include "threadpool/boost/threadpool.hpp"
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>


using namespace std;

namespace fs = boost::filesystem;


namespace 
{

/*!
 * Determines if a signal has a name that means it is a "broad signal".
 */
bool has_broad_signal_name( const std::string& fname )
{
	if ( fname.substr(0,2) == "MM" || fname.substr(0,3) == "Lip" ) 
	{
		return true;
	}
	else
	{
		return false;
	}
}

} // namespace


/*!
 * \brief Initialise the internal variables of a Basis object to be ready for analysis.
 */
void tarquin::CBasis::Initialise(float ppm_ref)
{
	applyRef(ppm_ref);	
	makeBasisMatrix();
	makeGroupMatrix();
}


namespace tarquin
{

struct SimulationJob
{
	std::vector< std::vector<double> > metab;    // simulation parameters
	const CFID*                        fidMatch; // the FID we simulating to fit to
	CSignal*                           signal;   // the signal we are producing
	const Options*                     options;   // the signal we are producing
	std::string                        name;   // the signal name we are producing

	SimulationJob() :
		fidMatch(NULL),
		signal(NULL) { }

	void run()
	{
		signal_simulate_full(metab, *signal, *fidMatch, *options, name);
	}
};

} // namespace tarquin

// simulate using a standard internal basis set
bool tarquin::CBasis::Simulate(const CFID& fidMatch, const Options& options, CBoswell& log)
{
	m_strBasisPath = "internal";
	
    // the metabolites in the basis set
	std::vector< basis_vector_e > metabolites;

    if ( options.GetIntBasisSet() == PROTON_BRAIN )
    {
        metabolites.push_back( BV_ALA );
        metabolites.push_back( BV_ASP );
        metabolites.push_back( BV_CR );
        metabolites.push_back( BV_CRCH2 );
        metabolites.push_back( BV_GABA );
        metabolites.push_back( BV_GLC );
        metabolites.push_back( BV_GLN );
        metabolites.push_back( BV_GLU );
        metabolites.push_back( BV_GPC );
        metabolites.push_back( BV_GUA );
        metabolites.push_back( BV_INS );
        metabolites.push_back( BV_LAC );
        metabolites.push_back( BV_LIP09 );
        metabolites.push_back( BV_LIP13A );
        metabolites.push_back( BV_LIP13B );
        metabolites.push_back( BV_LIP20 );
        metabolites.push_back( BV_MM09 );
        metabolites.push_back( BV_MM12 );
        metabolites.push_back( BV_MM14 );
        metabolites.push_back( BV_MM17 );
        metabolites.push_back( BV_MM20 );
        metabolites.push_back( BV_NAA );
        metabolites.push_back( BV_NAAG );
        metabolites.push_back( BV_PCH );
        metabolites.push_back( BV_PCR );
        metabolites.push_back( BV_SCYLLO );
        metabolites.push_back( BV_TAU );
    }
    else if ( options.GetIntBasisSet() == PROTON_BRAIN_NO_PCR )
    {
        metabolites.push_back( BV_ALA );
        metabolites.push_back( BV_ASP );
        metabolites.push_back( BV_CR );
        metabolites.push_back( BV_CRCH2 );
        metabolites.push_back( BV_GABA );
        metabolites.push_back( BV_GLC );
        metabolites.push_back( BV_GLN );
        metabolites.push_back( BV_GLU );
        metabolites.push_back( BV_GPC );
        metabolites.push_back( BV_GUA );
        metabolites.push_back( BV_INS );
        metabolites.push_back( BV_LAC );
        metabolites.push_back( BV_LIP09 );
        metabolites.push_back( BV_LIP13A );
        metabolites.push_back( BV_LIP13B );
        metabolites.push_back( BV_LIP20 );
        metabolites.push_back( BV_MM09 );
        metabolites.push_back( BV_MM12 );
        metabolites.push_back( BV_MM14 );
        metabolites.push_back( BV_MM17 );
        metabolites.push_back( BV_MM20 );
        metabolites.push_back( BV_NAA );
        metabolites.push_back( BV_NAAG );
        metabolites.push_back( BV_PCH );
        metabolites.push_back( BV_SCYLLO );
        metabolites.push_back( BV_TAU );
    }
    else if ( options.GetIntBasisSet() == PROTON_BRAIN_GLTH )
    {
        metabolites.push_back( BV_ALA );
        metabolites.push_back( BV_ASP );
        metabolites.push_back( BV_CR );
        metabolites.push_back( BV_CRCH2 );
        metabolites.push_back( BV_GABA );
        metabolites.push_back( BV_GLC );
        metabolites.push_back( BV_GLN );
        metabolites.push_back( BV_GLTH );
        metabolites.push_back( BV_GLU );
        metabolites.push_back( BV_GPC );
        metabolites.push_back( BV_INS );
        metabolites.push_back( BV_LAC );
        metabolites.push_back( BV_LIP09 );
        metabolites.push_back( BV_LIP13A );
        metabolites.push_back( BV_LIP13B );
        metabolites.push_back( BV_LIP20 );
        metabolites.push_back( BV_MM09 );
        metabolites.push_back( BV_MM12 );
        metabolites.push_back( BV_MM14 );
        metabolites.push_back( BV_MM17 );
        metabolites.push_back( BV_MM20 );
        metabolites.push_back( BV_NAA );
        metabolites.push_back( BV_NAAG );
        metabolites.push_back( BV_PCH );
        metabolites.push_back( BV_PCR );
        metabolites.push_back( BV_SCYLLO );
        metabolites.push_back( BV_TAU );
    }
    else if ( options.GetIntBasisSet() == PROTON_BRAIN_GLY_GLTH )
    {
        metabolites.push_back( BV_ALA );
        metabolites.push_back( BV_ASP );
        metabolites.push_back( BV_CR );
        metabolites.push_back( BV_CRCH2 );
        metabolites.push_back( BV_GABA );
        metabolites.push_back( BV_GLC );
        metabolites.push_back( BV_GLN );
        metabolites.push_back( BV_GLTH );
        metabolites.push_back( BV_GLU );
        metabolites.push_back( BV_GLY );
        metabolites.push_back( BV_GPC );
        metabolites.push_back( BV_INS );
        metabolites.push_back( BV_LAC );
        metabolites.push_back( BV_LIP09 );
        metabolites.push_back( BV_LIP13A );
        metabolites.push_back( BV_LIP13B );
        metabolites.push_back( BV_LIP20 );
        metabolites.push_back( BV_MM09 );
        metabolites.push_back( BV_MM12 );
        metabolites.push_back( BV_MM14 );
        metabolites.push_back( BV_MM17 );
        metabolites.push_back( BV_MM20 );
        metabolites.push_back( BV_NAA );
        metabolites.push_back( BV_NAAG );
        metabolites.push_back( BV_PCH );
        metabolites.push_back( BV_PCR );
        metabolites.push_back( BV_SCYLLO );
        metabolites.push_back( BV_TAU );
    }
    else if ( options.GetIntBasisSet() == PROTON_BRAIN_GLY_CIT_GLTH )
    {
        metabolites.push_back( BV_ALA );
        metabolites.push_back( BV_ASP );
        metabolites.push_back( BV_CIT );
        metabolites.push_back( BV_CR );
        metabolites.push_back( BV_CRCH2 );
        metabolites.push_back( BV_GABA );
        metabolites.push_back( BV_GLC );
        metabolites.push_back( BV_GLN );
        metabolites.push_back( BV_GLTH );
        metabolites.push_back( BV_GLU );
        metabolites.push_back( BV_GLY );
        metabolites.push_back( BV_GPC );
        metabolites.push_back( BV_INS );
        metabolites.push_back( BV_LAC );
        metabolites.push_back( BV_LIP09 );
        metabolites.push_back( BV_LIP13A );
        metabolites.push_back( BV_LIP13B );
        metabolites.push_back( BV_LIP20 );
        metabolites.push_back( BV_MM09 );
        metabolites.push_back( BV_MM12 );
        metabolites.push_back( BV_MM14 );
        metabolites.push_back( BV_MM17 );
        metabolites.push_back( BV_MM20 );
        metabolites.push_back( BV_NAA );
        metabolites.push_back( BV_NAAG );
        metabolites.push_back( BV_PCH );
        metabolites.push_back( BV_PCR );
        metabolites.push_back( BV_SCYLLO );
        metabolites.push_back( BV_TAU );
    }
    else if ( options.GetIntBasisSet() == PROTON_BRAIN_FULL )
    {
        metabolites.push_back( BV_ALA );
        metabolites.push_back( BV_ASP );
        metabolites.push_back( BV_CIT );
        metabolites.push_back( BV_CR );
        metabolites.push_back( BV_CRCH2 );
        metabolites.push_back( BV_GABA );
        metabolites.push_back( BV_GLC );
        metabolites.push_back( BV_GLN );
        metabolites.push_back( BV_GLTH );
        metabolites.push_back( BV_GLU );
        metabolites.push_back( BV_GLY );
        metabolites.push_back( BV_GPC );
        metabolites.push_back( BV_INS );
        metabolites.push_back( BV_LAC );
        metabolites.push_back( BV_LIP09 );
        metabolites.push_back( BV_LIP13A );
        metabolites.push_back( BV_LIP13B );
        metabolites.push_back( BV_LIP20 );
        metabolites.push_back( BV_MM09 );
        metabolites.push_back( BV_MM12 );
        metabolites.push_back( BV_MM14 );
        metabolites.push_back( BV_MM17 );
        metabolites.push_back( BV_MM20 );
        metabolites.push_back( BV_MM38 );
        metabolites.push_back( BV_NAA );
        metabolites.push_back( BV_NAAG );
        metabolites.push_back( BV_PCH );
        metabolites.push_back( BV_PCR );
        metabolites.push_back( BV_PETH );
        metabolites.push_back( BV_SCYLLO );
        metabolites.push_back( BV_TAU );
    }
    else if ( options.GetIntBasisSet() == PROTON_BRAIN_LE )
    {
        metabolites.push_back( BV_CR );
        metabolites.push_back( BV_CRCH2 );
        metabolites.push_back( BV_LAC );
        metabolites.push_back( BV_LIP09 );
        metabolites.push_back( BV_LIP13A );
        metabolites.push_back( BV_LIP13B );
        metabolites.push_back( BV_PCR );
        metabolites.push_back( BV_TCHO );
        metabolites.push_back( BV_TNAA );
    }
    else if ( options.GetIntBasisSet() == PROTON_MEGAPRESS_GABA )
    {
        metabolites.push_back( BV_MM09_MEGA );
        metabolites.push_back( BV_GABA_A );
        metabolites.push_back( BV_GABA_B );
        metabolites.push_back( BV_NAA_MEGA );
        metabolites.push_back( BV_GLX_A );
        metabolites.push_back( BV_GLX_B );
        metabolites.push_back( BV_GLX_C );
        metabolites.push_back( BV_GLX_D );
    }
    else if ( options.GetIntBasisSet() == PROTON_BRAINO )
    {
        metabolites.push_back( BV_CHO_RT );
        metabolites.push_back( BV_CR_RT );
        metabolites.push_back( BV_CRCH2_RT );
        metabolites.push_back( BV_GLU_RT );
        metabolites.push_back( BV_INS_RT );
        metabolites.push_back( BV_LAC_RT );
        metabolites.push_back( BV_NAA_RT );
    }
    else if ( options.GetIntBasisSet() == PHOSPH_BRAIN_DECOUP )
    {
        metabolites.push_back( BV_31P_ATP );
		metabolites.push_back( BV_31P_GPC );
		metabolites.push_back( BV_31P_GPE );
		metabolites.push_back( BV_31P_NADH );
		metabolites.push_back( BV_31P_PCH );
		metabolites.push_back( BV_31P_PCR );
		metabolites.push_back( BV_31P_PE );
		metabolites.push_back( BV_31P_PI );
    }
    else
    {
		log.Out(LOG_ERROR) << "Requested internal basis set not found";
        return false;
    }
    
	// these contain details for each metabolite
	m_signals.resize( metabolites.size() );
	m_vecSignalFiles.resize( metabolites.size() );
	m_broad_sig.resize( metabolites.size() );

	std::vector<SimulationJob> jobs( metabolites.size() );

	// simulate each metabolite
	for( size_t i = 0; i < metabolites.size(); ++i )
	{
		// get the description/name
		std::string fname;
	  	getMetaboliteDescription(metabolites[i], fname);
		
		// get the hardcoded parameters
		getMetaboliteMatrix(metabolites[i], jobs[i].metab, fidMatch.GetTransmitterFrequency());

		// store the filename 
		m_vecSignalFiles[i] = fs::path(fname);

		// is it broad?
		m_broad_sig[i] = has_broad_signal_name(fname);
		
		// parameters for the job
		jobs[i].fidMatch = &fidMatch;
		jobs[i].signal   = &m_signals[i];
		jobs[i].name     = fname;
		jobs[i].options  = &options;
	}

	int num_threads;
    if ( options.GetThreads() == 0 )
        num_threads = boost::thread::hardware_concurrency();
    else
        num_threads = options.GetThreads();

	boost::threadpool::pool thread_pool(num_threads);

	log.LogMessage(LOG_INFO, "Simulating using %d threads.", num_threads);

	//boost::posix_time::ptime time_start = boost::posix_time::microsec_clock::local_time();

	for( size_t i = 0; i < jobs.size(); ++i )
	{
		boost::threadpool::schedule(thread_pool, boost::bind(&SimulationJob::run, &jobs[i]));
	}

	thread_pool.wait();
	log.LogMessage(LOG_INFO, "All simulation threads stopped.");

	//boost::posix_time::ptime time_end = boost::posix_time::microsec_clock::local_time();
	//boost::posix_time::time_duration dur = time_end - time_start; 
	//std::cout << "\nSimulation took: " << dur.total_milliseconds() << std::flush;

	makeBasisMatrix();
	makeGroupMatrix();

	return true;
}

// simulate from a csv dir
bool tarquin::CBasis::Simulate(std::string strBasisPath, const CFID& fidMatch, const Options& options, CBoswell& log)
{
	m_strBasisPath = strBasisPath;
	//
	// First step, find a list of files from which to simulate.
	//
	// boost will enumerate this directory for us
	fs::path full_path(m_strBasisPath);

	// points passed the last file in the directory
	fs::directory_iterator end_iter;

	try {
		// for each file in the directory get a path object for it
		for( fs::directory_iterator dir_itr(full_path); dir_itr != end_iter; dir_itr++) {
			fs::path current = *dir_itr;

			std::string current_str = current.string();

			// csv file and not a directory so store this file
			if( false == fs::is_directory(*dir_itr) && ( current_str.substr(current_str.size()-4,4) == ".csv" || current_str.substr(current_str.size()-4,4) == ".dpt" ) ) 
				m_vecSignalFiles.push_back(*dir_itr);
		}
	}
	catch ( const exception& ex ) {

		log.Out(LOG_ERROR) << "basis directory does not exist or is empty: " << ex.what();
		return false;
	}

	// allocate space for the number of signals
	m_signals.resize( m_vecSignalFiles.size() );

	//
	// Third step, call the simulator for each file.
	//
	integer nMetabolite = 1;
	for( vector<fs::path>::iterator itFile = m_vecSignalFiles.begin(); 
			itFile != m_vecSignalFiles.end(); itFile++ ) {

		std::string strMessage = "Simulating '" + itFile->string() + "'";
		log.BeginTask(strMessage);
        
        // find the filename	
		int last_sep = itFile->string().find_last_of(filesep);
		std::string fname = itFile->string().substr(last_sep+1,itFile->string().size()); 

		boost::filesystem::path full_path = fname;

        if ( itFile->string().substr(itFile->string().size()-4,4) == ".csv" )
        {

            CCSVFile file;
            if( false == file.load(itFile->string()) )
                return false;

            log.UpdateTask(".");

            // read in CSV file
            std::vector<std::vector<double> >& file_doubmat = file.getDoubleMatrix();

            // print the doubmat
            /*
               for (vector< vector<double> >::size_type u = 0; u < file_doubmat.size(); u++) {
               for (vector<double>::size_type v = 0; v < file_doubmat[u].size(); v++) {
               cout << file_doubmat[u][v] << " ";
               }
               cout << endl;
               }
             */


            signal_simulate_full(file_doubmat, m_signals[nMetabolite-1], fidMatch, options, fname);
        }
        else // looks like a dpt file
        {
            // read in the dpt fid
            CFID dpt_fid;
            Options dummy_opts;

            CFIDReaderDPT reader(dpt_fid, log);
		    reader.Load(itFile->string(), dummy_opts, log);

            m_signals[nMetabolite-1].m_fids.resize(1);
            m_signals[nMetabolite-1].m_fids[0].AppendFromVector(dpt_fid.GetVectorFID()[0]);

            m_signals[nMetabolite-1].m_fids[0].SetSamplingFrequency(fidMatch.GetSamplingFrequency());
            m_signals[nMetabolite-1].m_fids[0].SetTransmitterFrequency(fidMatch.GetTransmitterFrequency());


            pair_vec default_ref;
            default_ref.push_back(std::make_pair(4.65, true));
            m_signals[nMetabolite-1].m_fids[0].SetPPMRef(default_ref);

        }

		m_broad_sig.push_back(has_broad_signal_name(full_path.filename().string()));

		log.UpdateTask(".");
		//std::cout << fidMatch.GetPPMRef() << std::endl;

		// the next file will go into this column of m_matSignals
		nMetabolite++;
		log.EndTask("done.");
	}

    if ( !makeBasisMatrix() )
        return false;

	makeGroupMatrix();

	return true;
}

bool tarquin::CBasis::SimulateCSV(std::string strBasisPath, const Options& options, const CFID& fidMatch)
{
	m_strBasisPath = strBasisPath;

	fs::path full_path(strBasisPath);

	m_signals.resize( 1 );

	CCSVFile file;
	if( false == file.load(full_path.string()) )
		return false;

	// read in CSV file
	std::vector<std::vector<double> >& file_doubmat = file.getDoubleMatrix();

	// find the filename	
	int last_sep = file.getFileName().find_last_of(filesep);
	std::string fname = file.getFileName().substr(last_sep+1,file.getFileName().size()); 

	signal_simulate_full(file_doubmat,m_signals[0], fidMatch, options, fname);
	m_broad_sig.push_back( has_broad_signal_name(fname) );

	if ( !makeBasisMatrix() )
        return false;

	makeGroupMatrix();

	return true;
}

bool tarquin::CBasis::makeBasisMatrix()
{
    if ( m_signals.size() < 1)
        return false;

    if ( m_signals[0].size() < 1)
        return false;

	assert( m_signals.size() > 0 );
	assert( m_signals[0].size() > 0 );

	// get the first fid in the collection of signals
	CFID& fid = *(m_signals[0].begin()); 

	// the number of points in the signal
	integer N = fid.GetNumberOfPoints();

	// the number of basis vectors
	integer M = m_signals.size();

	integer m = 1;

	m_matBasis.resize(N, M);
	m_matBasis.set(0);

	// for each signal (basis vector) 
	for( basis_iterator itB = m_signals.begin(); itB != m_signals.end(); itB++ ) {

		// for each fid in the current basis vector
		for( CSignal::fid_iterator itF = itB->begin(); itF != itB->end(); itF++ ) {

			cvm::cvector s = itF->GetVectorFID(0);
			m_matBasis(m) += s;
			//std::cout << s;
		}

		// column (e.g. metabolite)
		m++;
	}
	//plot(fft(m_matBasis));
    	//
	//m_matBasis.set(0);
    return true;
}

// also makes summation matrix
void tarquin::CBasis::makeGroupMatrix()
{
	assert( m_signals.size() > 0 );
	assert( m_signals[0].size() > 0 );

	// get the first fid in the collection of signals
	CFID& fid = *(m_signals[0].begin()); 

	// the number of points in the signal
	integer N = fid.GetNumberOfPoints();

	// the number of groups
	integer M = 0;
	for( basis_iterator itB = m_signals.begin(); itB != m_signals.end(); itB++ ) 
		for( CSignal::fid_iterator itF = itB->begin(); itF != itB->end(); itF++ ) 
			M++;

	m_matGroups.resize(N, M);
	m_matSum.resize(M, m_matBasis.nsize() );
	integer m = 1;

	integer q = 1;

	// for each signal (basis vector) 
	for( basis_iterator itB = m_signals.begin(); itB != m_signals.end(); itB++ ) {

		// for each fid in the current basis vector
		for( CSignal::fid_iterator itF = itB->begin(); itF != itB->end(); itF++ ) {
			m_matGroups(m) = itF->GetVectorFID(0);
			m_matSum(m,q) = 1;

			// 'index_tracker' in matlab code
			m_indexTracker.push_back(q);

			// advance column of group matrix == row of summation matrix
			m++;
		}

		// column of summation matrix == column of basis matrix
		q++;
	}

	//cvm::cmatrix A = (m_matBasis - m_matGroups * m_matSum);
	//std::cout << "\ndifference is " << A.norm1() << std::endl;
    
    /*
    cvm::cmatrix plot_mat = fft(m_matGroups);
    plot_mat = fftshift(plot_mat);
    plot(plot_mat);
    */

}

void tarquin::CBasis::SaveLCM(string file_name, CFID fid)
{
    std::ofstream basisfile;
    basisfile.open(file_name.c_str());

    int zf = 1;

    basisfile << setiosflags(ios::uppercase);
    std::vector<std::string> sig_names = GetSignalNames();

    coord voxel(1, 1, 1);

    basisfile << " $SEQPAR" << std::endl;
    // assume a metab alpha of 2.5
    basisfile << " FWHMBA =  " << 2.5/M_PI/GetTransmitterFrequency()*1.0e6 << " ," << std::endl;
    basisfile << " HZPPPM =  " << setprecision(8) << fid.GetTransmitterFrequency()*1.0e-6  << "," << std::endl;
    basisfile << " ECHOT =  " << fid.GetEchoTime()*1000.0 << "," << std::endl;
    basisfile << " SEQ = 'PRESS' $END" << std::endl;
    
    basisfile << " $BASIS1" << std::endl;
    basisfile << " IDBASI = 'TARQUIN simulated'," << std::endl;
    basisfile << " FMTBAS = '(6E13.5)'," << std::endl;
    basisfile << " BADELT =  " << 1.0/fid.GetSamplingFrequency() << "," << std::endl;
    basisfile << " NDATAB = " << zf * fid.GetNumberOfPoints() << " $END" << std::endl;


    basisfile << std::scientific;
    // metabolite
    for ( int n = 0; n < m_matBasis.nsize(); n++ )
    {
        //if ( !has_broad_signal_name(sig_names[n]) && sig_names[n] != "-CrCH2" )
        if ( sig_names[n] != "-CrCH2" )
        {
            basisfile << " $NMUSED" << std::endl;
            basisfile << " FILERAW = \'" << sig_names[n] << "_simulated" << "\'," << std::endl;
            basisfile << " METABO_CONTAM = ' '," << std::endl;
            basisfile << " METABO_SINGLET = ' '," << std::endl;
            basisfile << " AUTOPH = F," << std::endl;
            basisfile << " AUTOSC = F," << std::endl;
            basisfile << " DO_CONTAM = F," << std::endl;
            basisfile << " FLATEN = F," << std::endl;
            basisfile << " CONCSC = -1.," << std::endl;
            basisfile << " DEGZER =  0.," << std::endl;
            basisfile << " DEGPAP =  0.," << std::endl;
            basisfile << " DEGPPM =  0.," << std::endl;
            basisfile << " HWDPHA =  0.100000001," << std::endl;
            basisfile << " HWDSCA =  0.100000001," << std::endl;
            basisfile << " PPMAPP =  0.5 -0.5," << std::endl;
            basisfile << " PPMAPP_CONTAM =  0.  0.," << std::endl;
            basisfile << " PPMBAS =  0.00999999978," << std::endl;
            basisfile << " PPMPK =  0.00999999978," << std::endl;
            basisfile << " PPMPK_CONTAM =  0.," << std::endl;
            basisfile << " PPMPHA =  0.," << std::endl;
            basisfile << " PPMSCA =  8.43999958," << std::endl;
            basisfile << " PPM_SPLIT = -999.," << std::endl;
            basisfile << " RINTEG =  0.," << std::endl;
            basisfile << " SDPNTS =  1.," << std::endl;
            basisfile << " XTRASH =  0. $END" << std::endl;

            basisfile << " $BASIS" << std::endl;
            basisfile << " ID = \'" << sig_names[n] << "\'," << std::endl;
            basisfile << " METABO = \'" << sig_names[n] << "\'," << std::endl;
            basisfile << " CONC = 1.0," << std::endl;
            basisfile << " TRAMP = 1.0," << std::endl;
            basisfile << " VOLUME = 1.0," << std::endl;
            basisfile << " ISHIFT = 0 $END" << std::endl;
            // data point
            size_t cnt = 0;
            cvm::cvector metab = m_matBasis(n+1);
            metab.resize(metab.size() * zf);
            cvm::cvector METAB = fft(metab);

            //std::cout << " METABO = \'" << sig_names[n] << "\'," << std::endl << std::flush;
            /*plot(metab);
            plot(METAB);*/

            for ( int m = 0; m < m_matBasis.msize(); m++ )
            {
                //basisfile << setw(13) << setprecision(5) << m_matBasis(m+1,n+1).real() << setw(13) << setprecision(5) << m_matBasis(m+1,n+1).imag();
                basisfile << setw(13) << setprecision(5) << METAB(m+1).real() << setw(13) << setprecision(5) << METAB(m+1).imag();
                if ( cnt == 2 || m == m_matBasis.msize()-1 ) 
                {
                    basisfile << std::endl;
                    cnt = 0;
                }
                else
                    cnt++;
            }
        }
    }
    basisfile.close();
}

void tarquin::CBasis::SaveTxt(string file_name, CFID fid)
{
    std::ofstream basisfile;
    basisfile.open(file_name.c_str());
  
    std::vector<std::string> sig_names = GetSignalNames();

    coord voxel(1, 1, 1);
    cvm::rvector time_scale = fid.GetTimeScale();

    //std::cout << std::endl << sig_names.size() << std::endl;
    //std::cout << m_matBasis.nsize() << std::endl;
    //std::cout << m_matBasis.msize() << std::endl;
    //cvm::cmatrix zf_m_matBasis = m_matBasis;
    //zf_m_matBasis.resize(m_matBasis.msize()*2,m_matBasis.nsize());
    //cvm::cmatrix plot_mat = fft(zf_m_matBasis);
    //plot_mat = fftshift(plot_mat);

    basisfile << "time,";
    for ( int n = 0; n < m_matBasis.nsize()-1; n++ )
        basisfile << sig_names[n] << ",";
    basisfile << sig_names[sig_names.size()-1];
    basisfile << std::endl;
    
    for ( int m = 0; m < m_matBasis.msize(); m++ )
    {
        for ( int n = 0; n < m_matBasis.nsize()-1; n++ )
        {
            if ( n == 0 )
                basisfile << time_scale(m+1) << ",";
            basisfile << m_matBasis(m+1,n+1).real() << std::showpos << m_matBasis(m+1,n+1).imag() << std::noshowpos << "j,";
        }
        basisfile << m_matBasis(m+1,m_matBasis.nsize()).real() << std::showpos << m_matBasis(m+1,m_matBasis.nsize()).imag() << std::noshowpos << "j";
        basisfile << std::endl;
    }
    basisfile.close();
}

void tarquin::CBasis::CompressGroupMatrix()
{
	// resize the group matrix and set to be identical to the basis matrix	
	m_matGroups.resize(m_matBasis.msize(), m_matBasis.nsize());
	m_matGroups = m_matBasis;

	// resize the sum matrix and set to be an identity matrix	
	m_matSum.resize(m_matGroups.nsize(), m_matGroups.nsize());
	m_matSum = cvm::eye_complex(m_matGroups.nsize());

	// update the index tracker
	for( int n = 0 ; n < m_matGroups.nsize(); n++ )
		m_indexTracker[n] = n+1;
}

void tarquin::CBasis::applyRef(treal new_ref) 
{
	// for each signal (basis vector) 
	for( basis_iterator itB = m_signals.begin(); itB != m_signals.end(); itB++ ) 
	{

		// for each fid in the current basis vector
		for( CSignal::fid_iterator itF = itB->begin(); itF != itB->end(); itF++ ) 
			itF->ShiftRef(new_ref);

	}
}

void tarquin::CBasis::GenerateSTLstruc_pb() 
{
	// convert filepaths to vector of strings
	//std::cout << m_vecSignalFiles.size();
	// get vector of names 
	std::string signal_name;
	// clear the STL vector
	m_vecSignalFilesSTL.clear();
	for( size_t n = 0 ; n < m_vecSignalFiles.size() ; n++ )
	{	
		signal_name = m_vecSignalFiles[n].string();
		m_vecSignalFilesSTL.push_back(signal_name);
	}
}

void tarquin::CBasis::GenerateNonSTLstruc_pb() 
{
	makeBasisMatrix();
	makeGroupMatrix();

	//std::string signal_name;
	boost::filesystem::path signal_name;
	//std::cout << m_vecSignalFilesSTL.size();
	for( size_t n = 0 ; n < m_vecSignalFilesSTL.size() ; n++ )
	{	
		signal_name = m_vecSignalFilesSTL[n];
		m_vecSignalFiles.push_back(signal_name);
	}

}

bool tarquin::CBasis::check(const CFID& fid, CBoswell& log)
{
	// the basis is innocent until proven guilty
	bool ret_flag = true;

	if ( abs( fid.GetNumberOfPoints() - GetNumberOfPoints()) > 0 )
	{
		log.LogMessage(LOG_ERROR, 
				"incorrect number of basis set data points (%d in FID does not match %d in basis)", 
				fid.GetNumberOfPoints(), GetNumberOfPoints());

		ret_flag = false;
	}

	if ( abs( fid.GetSamplingFrequency() - GetSamplingFrequency() ) > 0.1 )
	{
		log.LogMessage(LOG_ERROR, 
				"incorrect sampling frequency (%.2f in FID does not match %.2f in basis)",
				fid.GetSamplingFrequency(), GetSamplingFrequency());

		ret_flag = false;
	}

	// setting the pulse sequence variable is a little sketchy
	// so skip this for now, TODO

	/*if ( fid.GetPulseSequence() != GetPulseSequence() )
	  {
	  log.Out(LOG_ERROR) << "incorrect basis set pulse sequence";
	  ret_flag = false;
	  }
	  std::cout << "FID pulse seq   :" << fid.GetPulseSequence() << std::endl;
	  std::cout << "Basis pulse seq :" << GetPulseSequence() << std::endl;
	 */

	if ( abs( fid.GetEchoTime() - GetEchoTime() ) > 1 )
	{
		log.LogMessage(LOG_ERROR, 
				"incorrect echo time (%.2f in FID does not match %.2f in basis)",
				fid.GetEchoTime(),
				GetEchoTime());
		ret_flag = false;
	}

	if ( abs( fid.GetTransmitterFrequency() - GetTransmitterFrequency() ) > ( 1.0/1000.0*fid.GetTransmitterFrequency() ) )
	{
		log.LogMessage(LOG_ERROR, "incorrect basis set field strength (%.2f in FID does not match %.2f in basis",fid.GetTransmitterFrequency(),GetTransmitterFrequency());
		ret_flag = false;
	}

	return ret_flag;
}


bool tarquin::CBasis::ReadLCMBasis(std::string strBasisPath, const CFID& fidMatch, const Options& options, CBoswell& log)
{
    // attempt to open file
    std::ifstream fin(strBasisPath.c_str());

    if( false == fin.is_open() )
    {
        std::cout << "Error" << std::endl;
        //Exception e("failed to open file: %s", strBasisPath.c_str()); 
        //throw e;
    }

    size_t nLine = 0;
    
    // main header
    bool in_seqpar = false;
    bool in_basis1 = false;

    // repeated for each metab
    bool in_nmused = false;
    bool in_basis = false;
    bool in_fid = false;

    CFID fid;
    std::vector<std::complex<double> > std_spec;
    size_t num_pts = 0;
    double fs = 0;
    double ft = 0;
    double te = 0;

    // read line by line 
    for( std::string strLine; getline(fin, strLine); )
    {
        nLine++;

        // skip over empty lines
        if( 0 == strLine.length() )
            continue;

        // skip over comments
        if( '#' == strLine.at(0) )
            continue;

        // trim whitespace before and after
        boost::trim(strLine);

        if ( in_seqpar )
        {
            
            std::vector<string> tokens;
            boost::split(tokens, strLine, boost::is_any_of("="));

            std::string para_name = tokens[0];
            boost::trim(para_name);
            if ( para_name == "HZPPPM" )
            {
                std::string para = tokens[1];
                boost::trim(para);

                std::vector<string> tokens_sp;
                boost::split(tokens_sp, para, boost::is_any_of(" "));
                para = tokens_sp[0];

                boost::erase_all(para, "\'");
                boost::erase_all(para, ",");
                //std::cout << para << std::endl;

                std::istringstream strmValue;
                strmValue.clear();
                strmValue.str(para);
                strmValue >> ft;
                ft = ft * 1e6;
            }
            else if ( para_name == "ECHOT" )
            {
                std::string para = tokens[1];
                boost::trim(para);
                
                std::vector<string> tokens_sp;
                boost::split(tokens_sp, para, boost::is_any_of(" "));
                para = tokens_sp[0];

                boost::erase_all(para, "\'");
                boost::erase_all(para, ",");
                //std::cout << para << std::endl;

                std::istringstream strmValue;
                strmValue.clear();
                strmValue.str(para);
                strmValue >> te;
                te = te/1000.0;
            }

            std::vector<string> end_tokens;
            boost::split(end_tokens, strLine, boost::is_any_of(" "));
            if ( std::find(end_tokens.begin(), end_tokens.end(), "$END") != end_tokens.end() )
                in_seqpar = false;
        }
        else if ( in_basis1 )
        {
            std::vector<string> tokens;
            boost::split(tokens, strLine, boost::is_any_of("="));

            std::string para_name = tokens[0];
            boost::trim(para_name);

            if ( para_name == "BADELT" )
            {
                std::string para = tokens[1];
                boost::trim(para);
                
                std::vector<string> tokens_sp;
                boost::split(tokens_sp, para, boost::is_any_of(" "));
                para = tokens_sp[0];

                boost::erase_all(para, "\'");
                boost::erase_all(para, ",");
                //std::cout << para << std::endl;
                
                std::istringstream strmValue;
                strmValue.clear();
                strmValue.str(para);
                strmValue >> fs;
                fs = 1.0/fs;
            }
            else if ( para_name == "NDATAB" )
            {
                std::string para = tokens[1];
                boost::trim(para);

                std::vector<string> tokens_sp;
                boost::split(tokens_sp, para, boost::is_any_of(" "));
                para = tokens_sp[0];

                boost::erase_all(para, "\'");
                boost::erase_all(para, ",");
                //std::cout << para << std::endl;
                
                std::istringstream strmValue;
                strmValue.clear();
                strmValue.str(para);
                strmValue >> num_pts;
            }

            std::vector<string> end_tokens;
            boost::split(end_tokens, strLine, boost::is_any_of(" "));
            if ( std::find(end_tokens.begin(), end_tokens.end(), "$END") != end_tokens.end() )
                in_basis1 = false;
        }
        else if ( in_nmused )
        {
            // clear the FID
            fid = CFID();
            std_spec.clear();

            std::vector<string> end_tokens;
            boost::split(end_tokens, strLine, boost::is_any_of(" "));
            if ( std::find(end_tokens.begin(), end_tokens.end(), "$END") != end_tokens.end() )
                in_nmused = false;
        }
        else if ( in_basis )
        {

            std::vector<string> tokens;
            boost::split(tokens, strLine, boost::is_any_of("="));
            //std::cout << tokens[0] << std::endl;
            //std::cout << tokens[1] << std::endl;

            std::string para_name = tokens[0];
            boost::trim(para_name);
            if ( para_name == "METABO" )
            {
                std::string para = tokens[1];
                boost::trim(para);
                
                std::vector<string> tokens_sp;
                boost::split(tokens_sp, para, boost::is_any_of(" "));
                para = tokens_sp[0];

                boost::erase_all(para, "\'");
                boost::erase_all(para, ",");
                
                fid.SetFilename(para); // TODO add .csv?
                m_vecSignalFiles.push_back(fs::path(para));
                m_broad_sig.push_back(has_broad_signal_name(para));
            }
            else if ( para_name == "ISHIFT" )
            {
                std::string para = tokens[1];
                boost::trim(para);
                
                std::vector<string> tokens_sp;
                boost::split(tokens_sp, para, boost::is_any_of(" "));
                para = tokens_sp[0];
                
                boost::erase_all(para, "\'");
                boost::erase_all(para, ",");
                //std::cout << para << std::endl;
            }

            std::vector<string> end_tokens;
            boost::split(end_tokens, strLine, boost::is_any_of(" "));
            if ( std::find(end_tokens.begin(), end_tokens.end(), "$END") != end_tokens.end() )
            {
                in_basis = false;
                in_fid = true;
                //std::cout << "In FID" << std::endl;
            }
        }
        else if ( in_fid )
        {
            if ( strLine == "$NMUSED" )
            {
                in_fid = false;
                in_nmused = true;
                //std::cout << "In NMUSED" << std::endl;

                size_t N = std_spec.size();
                cvm::cvector cvm_spec(N);
                for ( size_t n = 0; n < N; n++ )
                    cvm_spec(n+1) = std_spec[n];
                    
                cvm::cvector cvm_fid(N);
                ifft(cvm_spec, cvm_fid);
                
                // resize cvm_fid to match data
                cvm_fid.resize(fidMatch.GetNumberOfPoints());
                
	            fid.AppendFromVector(cvm_fid);
                fid.SetSamplingFrequency(fs);
                fid.SetTransmitterFrequency(ft);
                fid.SetEchoTime(te);
                std::vector<CFID> fid_vec; 
                fid_vec.push_back(fid);
                CSignal signal;
                signal.SetFids(fid_vec);
                // append signal to basis
                m_signals.push_back(signal);
            }
            else
            {
                std::vector<string> tokens;
                boost::split( tokens, strLine, boost::is_any_of(" "), boost::token_compress_on );

                for ( std::vector<string>::iterator it_token = tokens.begin(); it_token != tokens.end(); it_token++ )
                {
	                std::istringstream ins;
			        double real = 0;
                    ins.str(*it_token);
			        ins >> real;
                    //std::cout << *it_token << std::endl;
                    //std::cout << real << std::endl;
                    
                    ins.clear();
                    it_token++;
			        double imag = 0;
                    ins.str(*it_token);
			        ins >> imag;
                    //std::cout << *it_token << std::endl;
                    //std::cout << imag << std::endl;

                    std_spec.push_back(tcomplex(real, imag));
                }
            }
        }
        else
        {
            if ( strLine == "$SEQPAR" )
            {
                in_seqpar = true;
                //std::cout << "In SEQPAR" << std::endl;
            }
            else if ( strLine == "$BASIS1" )
            {
                in_basis1 = true;
                //std::cout << "In BASIS1" << std::endl;
            }
            else if ( strLine == "$NMUSED" )
            {
                in_nmused = true;
                //std::cout << "In NMUSED" << std::endl;
            }
            else if ( strLine == "$BASIS" )
            {
                in_basis = true;
                //std::cout << "In BASIS" << std::endl;
            }
        }
    }
    // don't forget the last FID
    size_t N = std_spec.size();
    cvm::cvector cvm_spec(N);
    for ( size_t n = 0; n < N; n++ )
        cvm_spec(n+1) = std_spec[n];

    cvm::cvector cvm_fid(N);
    ifft(cvm_spec, cvm_fid);

    // resize cvm_fid to match data
    cvm_fid.resize(fidMatch.GetNumberOfPoints());

    fid.AppendFromVector(cvm_fid);
    fid.SetSamplingFrequency(fs);
    fid.SetTransmitterFrequency(ft);
    fid.SetEchoTime(te);
    std::vector<CFID> fid_vec; 
    fid_vec.push_back(fid);
    CSignal signal;
    signal.SetFids(fid_vec);
    // append signal to basis
    m_signals.push_back(signal);

    // set scaling factor TODO
    // check ref values? TODO

    // append lipid and MM signals
	std::vector< basis_vector_e > metabolites;
    metabolites.push_back( BV_LIP09 );
    metabolites.push_back( BV_LIP13A );
    metabolites.push_back( BV_LIP13B );
    metabolites.push_back( BV_LIP20 );
    metabolites.push_back( BV_MM09 );
    metabolites.push_back( BV_MM12 );
    metabolites.push_back( BV_MM14 );
    metabolites.push_back( BV_MM17 );
    metabolites.push_back( BV_MM20 );
    metabolites.push_back( BV_CRCH2 );
    
    // simulate each metabolite
	for( size_t i = 0; i < metabolites.size(); ++i )
	{
		// get the description/name
		std::string fname;
	  	getMetaboliteDescription(metabolites[i], fname);

		std::vector< std::vector<double> > metab;
		// get the hardcoded parameters
		getMetaboliteMatrix(metabolites[i], metab, fidMatch.GetTransmitterFrequency());

		// store the filename 
		m_vecSignalFiles.push_back(fs::path(fname));

		// is it broad?
		m_broad_sig.push_back(has_broad_signal_name(fname));

        CSignal signal;

        // simulate
		signal_simulate_full(metab, signal, fidMatch, options, fname);
        m_signals.push_back(signal);
    }

	makeBasisMatrix();
	makeGroupMatrix();

    return true;
}
