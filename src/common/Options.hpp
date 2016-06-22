#ifndef __COPTIONS__
#define __COPTIONS__

#include "common.hpp"
#include "CConstraint.hpp"
#include "CFID.hpp"
#include "CBasis.hpp"
#include "fidio/GEOptions.hpp"

namespace tarquin 
{
	//! Supported FID formats.
	enum fid_format_e 
	{ 
		DANGER = 0, 
		SIEMENS, 
		PHILIPS, 
		GE, 
		RDA, 
		LCM, 
		SHF, 
		VARIAN, 
		BRUKER, 
		DCM, 
		PHILIPS_DCM, 
        JMRUI_TXT,
		NOTSET
	};

    //! Chemical shift referencing mode
	enum chem_shift_ref_e
	{ 
		PROTON_NAA_CR_CHO_LIP = 0,
		PROTON_NAA_CR_CHO,
		PROTON_NAA_CHO,
		PROTON_CR_CHO,
		PROTON_NAA,
		PROTON_CR,
		PROTON_CHO,
		PROTON_H2O,
		PROTON_LIP,
		PHOSPH_PCR,
		PHOSPH_PCR_GAMMAATP
	};

    //! Supported basis sets.
	enum basis_set_e 
	{ 
		PROTON_BRAIN = 0,           // Std brain metabolites
		PROTON_BRAIN_GLTH,           // Std brain metabolites Glth
		PROTON_BRAIN_GLY_GLTH,           // Std brain metabolites Gly + Glth
		PROTON_BRAIN_GLY_CIT_GLTH,       // Std brain metabolites + Gly + Cit + Glth
		PROTON_BRAIN_FULL,       // Std brain metabolites + Gly + Glth + PEth
		PROTON_BRAIN_LE,       // Std brain metabolites + Gly + Cit
		PROTON_BRAIN_NO_PCR,           // Std brain metabolites w/o PCr
		PROTON_MEGAPRESS_GABA,
		PROTON_BRAINO,
		PROTON_BRAIN_MMEXP,          // Std brain metabolites 
		PROTON_BRAIN_METAB_ONLY,   // Std brain metabolites Glth, no lip mm's
		PHOSPH_BRAIN_DECOUP
	};

    //! Supported pulse sequences
	enum pul_seq_e 
	{ 
		PRESS = 0,
		STEAM,
		LASER,
		SEMI_LASER,
		PULSE_ACQUIRE,
		CPMG,
		SPIN_ECHO,
		MEGA_PRESS,
        SHAPED_PRESS, 
        PROFILE_PRESS
	};
    
    //! Dynamic averaging scheme
	enum dyn_av_e
	{ 
		DEFAULT = 0,
		NONE,
		ALL,
		SUBTRACT,
		ODD,
		EVEN
	};
    
    /*!
     * Class that holds the options used to run the TARQUIN algorithm, such as optimisation
     * constraints, etc. Other information e.g. samping frequency, is specific to each FID, 
     * so that is held in the FID. Thus you need both a "Options" instance and a "CFID" instance 
     * to run the algorithm.
     *
     * Once initialised, e.g. by the command line parser, you need to initialise it with
     * a particular FID and basis as well, by calling the "Initialise" function. This is necessary
     * because some constraints (e.g. phi1) are functions of the sampling frequency.
     *
     */
    class Options {

	// this function initialise this class from the command line
	friend bool ParseCommandLine(int argc, char* argv[], Options& options, CFID& fid);

	public:

	Options() 
	{
	    m_nStart = 0;
	    m_nEnd = 0;
	    m_max_iters = 75;
	    
        // 2pi incase phi0 is at pi radians 
        m_phi0_lower = -2*M_PI;
	    m_phi0_upper = +2*M_PI;
	    m_phi0_typ = 0;
        
	    m_phi1_lower = -20;
	    m_phi1_upper = +20;
	    m_phi1_typ = 0;

        m_alpha_metab_lower = 0.0;
	    m_alpha_metab_upper = 10.0;
	    m_alpha_metab_typ = 2.0;

        m_alpha_broad_lower = 0.0;
	    m_alpha_broad_upper = 50.0;
	    m_alpha_broad_typ = 2.0;

        m_max_beta = 5000;
        m_beta_scale = 1;

	    m_bUsePrecompiled = false;
	    m_format = NOTSET;
	    m_water_window = 45;
	    m_df_water_window = -std::numeric_limits<treal>::infinity();
	    m_conv_window_width = 0;
	    m_bAutoPhase = true;
	    m_bAutoRef = true;

        // nb m_ref value should be linked to the default value 
        // in ClearVectors function in CFID.hpp
	    m_ref = 4.65;

	    m_ref_spec = false;
	    m_bShowPreprocessed = false;
	    m_bWaterEddy = false;
	    m_basis_comp = true;
        m_bPrintParas = false;
        m_bNofit = false;
        m_bReadOnly = false;
        m_bReadWriteOnly = false;

        m_ppm_start = 0.2; // ppm left
        m_ppm_end = 4.0; // ppm right
        m_pause = true;
        m_lb = 0.0;
        m_lb_ref = 0.0;
        m_init_beta = 0.0;
        //m_lambda = 0.05; // old value
        //m_lambda = 10;
        m_lambda = 0.2;

        m_w_att = 0.7;
        m_w_conc = 35880;
        m_au_norm = true;

        m_max_dref = 0.1;
        m_title = "";
        
        // TODO serialise
	    m_svs_only = false;
	    //m_nStartTime = 0.008;
	    m_nStartTime = 0.01;
	    m_ref_signals = PROTON_NAA_CR_CHO_LIP;
	    m_dyn_ref_signals = PROTON_NAA_CR_CHO_LIP;
        m_int_basis_set = PROTON_BRAIN_GLTH;

	    m_ref_freq = std::numeric_limits<treal>::infinity();

        m_dyn_av = DEFAULT;
        m_dyn_av_w = DEFAULT;
        m_dyn_freq_corr = false;
        m_pdfc = false;
        m_pul_seq = PRESS;
        m_fit_rows = -1;
        m_fit_cols = -1;
        m_fit_slices = -1;
        m_max_metab_shift = 0.03;
        m_max_broad_shift = 0.1;
        m_lipid_filt_freq = 1.65;
        m_nt_l2_dp = 1e-6;
        m_nt_init_mu = 1e-3;
	    m_bCombinePreproc = false;
	    m_bLipidFilter = false;
	    m_bFullEcho = false;
	    m_bSwapRowCol = false;
        m_zfill_kspace = 1;
        m_filter_kspace = false;
        m_zero_fill = 2;
        m_baseline = 25;
	    m_bPdfExt = false;
	    m_bPdfStack = false;
	    m_bAdaptSp = false;
	    m_bAdaptEp = false;
        m_press_TE1 = 0.0126;
        m_steam_TM = 0.012;
        m_cpmg_N = 10;
        m_gnuplot_cex = 2;
        m_threads = 0;
        m_pre_hsvd = false;
        m_old_phase = false;
        m_min_range = 5;
        m_lcm_basis_pts = 4096;
        m_ff = true;
        m_pre_ws_shift = true;
        m_keep_pre_ws_shift = false;

        m_prepend_pts = 0;
        m_truncate_pts = 0;
        m_replace_fp = false;

        m_ext_csv_fit = true;
        m_hsvd_comps = 50;
        m_max_hsvd_pts = 1024;
        m_crlb_td = true;
        m_crlb_optim = false;
        m_nnls = true;
        m_soft_cons = true;
	    m_bPreFitPhase = true;
	    m_bPreFitShift = false;
	    m_bPreFitBl = false;
	    m_bAppendLCMBasis = false;
	    m_bAppendNegRefBasis = false;
        m_lineshape_corr = false;
        m_invert_even_pairs = false;
	}

	//Options(){}

	~Options()
	{

	}

	//! Call when the arguments have been initialised, i.e. loaded.
	void Initialise(const CBasis& basis, const CFID& fid, CBoswell& log)
	{
        if ( m_old_phase )
        {
		    // divide by 2pi because Greg's definition of phi1 is multiplied by 2pi
            m_phi1_lower = (-M_PI/4) / ( fid.GetSamplingFrequency() / 2 ) / ( 2 * M_PI);
            m_phi1_upper = (+M_PI/4) / ( fid.GetSamplingFrequency() / 2 ) / ( 2 * M_PI);
        }

        // convert from deg/ppm to phi1
        m_phi1_lower = m_phi1_lower / ( 180/M_PI * (fid.GetTransmitterFrequency() / 1.0e6) * 2.0 * M_PI );
        m_phi1_upper = m_phi1_upper / ( 180/M_PI * (fid.GetTransmitterFrequency() / 1.0e6) * 2.0 * M_PI );
		m_phi1_typ = 0;

		// resize constraints collection to be commensurate with group matrix
		m_constraints.resize( basis.GetGroupMatrix().nsize() );

		// set some sensible default options
		SetMetabAlphaLimits(m_alpha_metab_lower, m_alpha_metab_upper, basis);
		SetBroadAlphaLimits(m_alpha_broad_lower, m_alpha_broad_upper, basis);
		SetBetaLimits(0, m_max_beta);
		log.LogMessage(LOG_INFO, "Setting metab shift limit to %f ppm", m_max_metab_shift);
		log.LogMessage(LOG_INFO, "Setting broad shift limit to %f ppm", m_max_broad_shift);
		SetMetabShiftLimits(-m_max_metab_shift, m_max_metab_shift, fid, basis);
		SetBroadShiftLimits(-m_max_broad_shift, m_max_broad_shift, fid, basis);
	}

	//! Set the starting point of the range over which the analysis is run.
	void SetRangeStart(integer nStart)
	{
	    m_nStart = nStart;
	}

	double GetRefFreq()
	{
	    return m_ref_freq;
	}

	//! Set the starting time of the range over which the analysis is run.
	void SetRangeStartTime(double nStartTime)
	{
	    m_nStartTime = nStartTime;
	}

	//! Set the ending point of the range over which the analysis is run.
	void SetRangeEnd(integer nEnd)
	{
	    m_nEnd = nEnd;
	}

	//! Set the constraints on the adjustment of phi0.
	void SetZeroPhaseLimit(treal lower, treal upper)
	{
	    m_phi0_lower = lower;
	    m_phi0_upper = upper;
	}

	void SetZeroPhaseTyp(treal typ)
	{
	    m_phi0_typ = typ;
	}

	//! Set the constraints on the adjustment of phi1.
	void SetFirstPhaseLimit(treal lower, treal upper)
	{
	    m_phi1_lower = lower;
	    m_phi1_upper = upper;
	}

	void SetFirstPhaseTyp(treal typ)
	{
	    m_phi1_typ = typ;
	}

	//! Set the alpha limits for metabolite groups
	void SetMetabAlphaLimits(treal lower, treal upper, const CBasis& basis)
	{
        std::vector<bool> broad_sig = basis.GetBroadSig();
        int n = 0;
	    for( std::vector<CConstraint>::iterator it = m_constraints.begin(); 
		    it != m_constraints.end(); it++ ) {

        if ( ! broad_sig[basis.GetBasisFromGroup(n)-1] )
        {
            it->m_minAlpha = lower;
            it->m_maxAlpha = upper;
            //it->m_typAlpha = (upper-lower)/2.0;
            // shouldn't be zero otherise won't work very well
            it->m_typAlpha = m_alpha_metab_typ;
        }

        n++;
	    }
	}
    
    //! Set the alpha limits for lipid and mm groups 
	void SetBroadAlphaLimits(treal lower, treal upper, const CBasis& basis)
	{
        std::vector<bool> broad_sig = basis.GetBroadSig();
        int n = 0;
	    for( std::vector<CConstraint>::iterator it = m_constraints.begin(); 
		    it != m_constraints.end(); it++ ) {

        if ( broad_sig[basis.GetBasisFromGroup(n)-1] )
        {
            it->m_minAlpha = lower;
            it->m_maxAlpha = upper;
            //it->m_typAlpha = (upper-lower)/2.0;
            // shouldn't be zero otherise won't work very well
            it->m_typAlpha = m_alpha_metab_typ;
        }

        n++;
	    }
	}


	void SetBetaLimits(treal lower, treal upper)
	{
	    for( std::vector<CConstraint>::iterator it = m_constraints.begin(); 
		    it != m_constraints.end(); it++ ) {

		it->m_minBeta = lower;
		it->m_maxBeta = upper;
		//it->m_typBeta = 0.0;
		//it->m_typBeta = 5.0;
	    }
	}

    void SetMaxBeta(treal max_beta)
	{
	    m_max_beta = max_beta;
	}
    
    treal GetMaxBeta() const
	{
	    return m_max_beta;
	}

    void SetBetaScale(treal beta_scale)
	{
	    m_beta_scale = beta_scale;
	}
    
    treal GetBetaScale() const
	{
	    return m_beta_scale;
	}

    void SetInitBeta(treal init_beta)
	{
	    m_init_beta = init_beta;
	}

    void AppendInitBetaUsed(treal init_beta)
	{
	    m_init_beta_used.push_back(init_beta);
	}

    void SetLambda(treal lambda)
	{
	    m_lambda = lambda;
	}
    
    void SetInitMu(treal init_mu)
	{
	    m_nt_init_mu = init_mu;
	}
    
    void SetMaxMetabShift(treal val)
	{
	    m_max_metab_shift = val;
	}

    void SetMaxBroadShift(treal val)
	{
	    m_max_broad_shift = val;
	}
    
    treal GetMaxMetabShift() const
	{
	    return m_max_metab_shift;
	}

    treal GetMaxBroadShift() const
	{
	    return m_max_broad_shift;
	}
    
    treal GetHSVDComps() const
	{
	    return m_hsvd_comps;
	}

    treal GetMaxHSVDPts() const
	{
	    return m_max_hsvd_pts;
	}

    bool GetCRLB_TD() const
	{
	    return m_crlb_td;
	}

    bool GetOptimCRLBs() const
	{
	    return m_crlb_optim;
	}

    bool GetNNLS() const
	{
	    return m_nnls;
	}

    bool GetSoftCons() const
	{
	    return m_soft_cons;
	}

    bool GetLineshapeCorr() const
	{
	    return m_lineshape_corr;
	}

    bool GetInvertEvenPairs() const
	{
	    return m_invert_even_pairs;
	}

	void SetMetabShiftLimits(treal lower, treal upper, const CFID& fid, const CBasis& basis)
	{
        
        std::vector<bool> broad_sig = basis.GetBroadSig();
        /*
        for ( size_t i = 0; i < broad_sig.size(); i++ )
        {
            std::cout << basis.GetSignalName(i) << std::endl;
            if ( broad_sig[i] )
                std::cout << "True" << std::endl;
            else
                std::cout << "False" << std::endl;
        }*/

        int n = 0;
        for( std::vector<CConstraint>::iterator it = m_constraints.begin(); 
                it != m_constraints.end(); it++ ) 
        {
            //std::cout << basis.GetBasisFromGroup(n)-1 << std::endl;
            //std::cout << basis.GetSignalName(basis.GetBasisFromGroup(n)-1) << std::endl;
            if ( ! broad_sig[basis.GetBasisFromGroup(n)-1] )
            {
                //std::cout << broad_sig[basis.GetBasisFromGroup(n)-1] << std::endl;
                it->m_minShiftHz = fid.ConvertPPMtoHz(lower);
                it->m_maxShiftHz = fid.ConvertPPMtoHz(upper);
                it->m_typShiftHz = 0.0;
            }
            n++;
        }
	}

    void SetBroadShiftLimits(treal lower, treal upper, const CFID&  fid, const CBasis& basis)
    {
        std::vector<bool> broad_sig = basis.GetBroadSig();
        int n = 0;
        for( std::vector<CConstraint>::iterator it = m_constraints.begin(); 
                it != m_constraints.end(); it++ ) 
        {
            
            if ( broad_sig[basis.GetBasisFromGroup(n)-1] )
            {
                it->m_minShiftHz = fid.ConvertPPMtoHz(lower);
                it->m_maxShiftHz = fid.ConvertPPMtoHz(upper);
                it->m_typShiftHz = 0.0;
            }
            n++;
        }
    }

	std::vector<CConstraint>& GetConstraints()
    {
        return m_constraints;
    }

	treal GetLowerLimitAlpha(std::size_t nGroup) const
	{
	    assert( nGroup < m_constraints.size() );
	    return m_constraints[nGroup].m_minAlpha;
	}

	treal GetUpperLimitAlpha(std::size_t nGroup) const
	{
	    assert( nGroup < m_constraints.size() );
	    return m_constraints[nGroup].m_maxAlpha;
	}

	treal GetTypicalAlpha(std::size_t nGroup) const
	{
	    assert( nGroup < m_constraints.size() );
	    return m_constraints[nGroup].m_typAlpha;
	}

	treal GetLowerLimitBeta(std::size_t nGroup) const
	{
	    assert( nGroup < m_constraints.size() );
	    return m_constraints[nGroup].m_minBeta;
	}

	treal GetUpperLimitBeta(std::size_t nGroup) const
	{
	    assert( nGroup < m_constraints.size() );
	    return m_constraints[nGroup].m_maxBeta;
	}

	treal GetTypicalBeta(std::size_t nGroup) const
	{
	    assert( nGroup < m_constraints.size() );
	    return m_constraints[nGroup].m_typBeta;
	}

	treal GetLowerLimitShift(std::size_t nGroup) const
	{
	    assert( nGroup < m_constraints.size() );
	    return m_constraints[nGroup].m_minShiftHz;
	}

	treal GetUpperLimitShift(std::size_t nGroup) const
	{
	    assert( nGroup < m_constraints.size() );
	    return m_constraints[nGroup].m_maxShiftHz;
	}

	treal GetTypicalShift(std::size_t nGroup) const
	{
	    assert( nGroup < m_constraints.size() );
	    return m_constraints[nGroup].m_typShiftHz;
	}

	treal GetLowerLimitPhi0() const
	{
	    return m_phi0_lower;
	}

	treal GetUpperLimitPhi0() const
	{
	    return m_phi0_upper;
	}

	treal GetTypicalPhi0() const
	{
	    return m_phi0_typ;
	}

	treal GetLowerLimitPhi1() const
	{
	    return m_phi1_lower;
	}

	treal GetUpperLimitPhi1() const
	{
	    return m_phi1_upper;
	}

	treal GetTypicalPhi1() const
	{
	    return m_phi1_typ;
	}

	treal GetPPMstart() const
	{
	    return m_ppm_start;
	}

	void SetPPMstart(treal ppm)
	{
	    m_ppm_start = ppm;
	}

	treal GetPPMend() const
	{
	    return m_ppm_end;
	}

	void SetPPMend(treal ppm)
	{
	    m_ppm_end = ppm;
	}

	treal Getlb() const
	{
	    return m_lb;
	}

	void Setlb(treal lb)
	{
	    m_lb = lb;
	}

	treal GetlbRef() const
	{
	    return m_lb_ref;
	}

	void SetlbRef(treal lb_ref)
	{
	    m_lb_ref = lb_ref;
	}
    
    size_t GetZF() const
	{
	    return m_zero_fill;
	}

    void SetZF(size_t zf)
	{
	    m_zero_fill = zf;
	}
    
    size_t GetBL() const
	{
	    return m_baseline;
	}

    void SetBL(size_t baseline)
	{
	    m_baseline = baseline;
	}

    treal GetInitBeta() const
	{
	    return m_init_beta;
	}
    
    treal GetInitBetaUsed(size_t ind) const
	{
	    return m_init_beta_used[ind];
	}
    
    treal GetLambda() const
	{
	    return m_lambda;
	}

	integer GetRangeStart() const
	{
	    return m_nStart;
	}

	double GetRangeStartTime() const
	{
	    return m_nStartTime;
	}

	integer GetRangeEnd() const
	{
	    return m_nEnd;
	}
	
	integer GetMaxIters() const
	{
	    return m_max_iters;
	}

    void SetMaxIters(integer max_iters)
	{
	    m_max_iters = max_iters;
	}

    bool GetFilterKspace() const
	{
	    return m_filter_kspace;
	}

    void SetFilterKspace(bool filter)
	{
	    m_filter_kspace = filter;
	}

    size_t GetZfillKspace() const
    {
        return m_zfill_kspace;
    }

    void SetZfillKspace(size_t zfill_kspace)
    {
        m_zfill_kspace = zfill_kspace;
    }

	std::string GetFilename() const
	{
	    return m_strFile;
	}

	void SetFilename(std::string strFile) 
	{
	    m_strFile = strFile;
	}

	std::string GetDynShiftFilename() const
	{
	    return m_dyn_shift_file;
	}

	void SetDynShiftFilename(std::string strFile) 
	{
	    m_dyn_shift_file = strFile;
	}

	std::string GetFilenameWater() const
	{
	    return m_strFileWater;
	}

	void SetFilenameWater(std::string strFileWater) 
	{
	    m_strFileWater = strFileWater;
	}

	std::string GetFilenameImag() const
	{
	    return m_strFileOutImag;
	}
    
    void SetFilenameImag(std::string strFileImag) 
	{
	    m_strFileOutImag = strFileImag;
	}
    
    std::string GetFilenamePdf() const
	{
	    return m_strFileOutPdfImag;
	}
    
    void SetFilenamePdf(std::string strFileImag) 
	{
	    m_strFileOutPdfImag = strFileImag;
	}
    
	std::string GetFilenameTxt() const
	{
	    return m_strFileOutTxt;
	}
    
    void SetFilenameTxt(std::string strFileTxt) 
	{
	    m_strFileOutTxt = strFileTxt;
	}

	std::string GetFilenameCSV() const
	{
	    return m_strFileOutCSV;
	}

	std::string GetFilenameCSVGeom() const
	{
	    return m_strFileOutCSVGeom;
	}

    std::string GetFilenameCSVSpectraAligned() const
    {
	    return m_strFileOutCSVSpectraAligned;
    }

    std::string GetFilenameCSVSpectraAlignedMag() const
    {
	    return m_strFileOutCSVSpectraAlignedMag;
    }

	std::string GetFilenameCSVFit() const
	{
	    return m_strFileOutCSVFit;
	}

	std::string GetFilenameCSVFitMag() const
	{
	    return m_strFileOutCSVFitMag;
	}
    
    void SetFilenameCSV(std::string strFileCSV) 
	{
	    m_strFileOutCSV = strFileCSV;
	}
    
    void SetFilenameCSVFit(std::string strFileCSVFit) 
	{
	    m_strFileOutCSVFit = strFileCSVFit;
	}

    void SetFilenameCSVFitAbs(std::string strFileCSVFitMag) 
	{
	    m_strFileOutCSVFitMag = strFileCSVFitMag;
	}

	bool GetUsePrecompiled() const
	{
	    return m_bUsePrecompiled;
	}

	void SetUsePrecompiled(bool bState)
	{
	    m_bUsePrecompiled = bState;
	}

	std::string GetBasisPath() const
	{
	    return m_strBasisPath;
	}

	basis_set_e GetIntBasisSet() const
	{
	    return m_int_basis_set;
	}

	void SetIntBasisSet(basis_set_e basis)
	{
	    m_int_basis_set = basis;
	}

	pul_seq_e GetPulSeq() const
	{
	    return m_pul_seq;
	}

	void SetPulSeq(pul_seq_e pul_seq)
	{
	    m_pul_seq = pul_seq;
	}

	dyn_av_e GetDynAv() const
	{
	    return m_dyn_av;
	}

	void SetDynAv(dyn_av_e dyn_av)
	{
	    m_dyn_av = dyn_av;
	}

	dyn_av_e GetDynAvW() const
	{
	    return m_dyn_av_w;
	}

	void SetDynAvW(dyn_av_e dyn_av)
	{
	    m_dyn_av_w = dyn_av;
	}

	bool GetDynFreqCorr() const
	{
	    return m_dyn_freq_corr;
	}

	bool GetPDFC() const
	{
	    return m_pdfc;
	}

	void SetDynFreqCorr(bool dyn_freq_corr)
	{
	    m_dyn_freq_corr = dyn_freq_corr;
	}

	bool GetSVSOnly() const
    {
	    return m_svs_only;
    }

	void SetBasisPath(std::string strPath)
	{
	    m_strBasisPath = strPath;
	}

	std::string GetOutputXMLPath() const
	{
	    return m_strOutputXMLPath;
	}

	void SetOutputXMLPath(std::string strPath )
	{
	    m_strOutputXMLPath = strPath;
	}

	fid_format_e GetFormat() const
	{
	    return m_format;
	}

	void SetFormat(fid_format_e format)
	{
	    m_format = format;
	}

	chem_shift_ref_e GetRefSignals() const
	{
	    return m_ref_signals;
	}

	void SetRefSignals(chem_shift_ref_e signals)
	{
	    m_ref_signals = signals;
	}

	chem_shift_ref_e GetDynRefSignals() const
	{
	    return m_dyn_ref_signals;
	}

	void SetDynRefSignals(chem_shift_ref_e signals)
	{
	    m_dyn_ref_signals = signals;
	}

	treal GetWaterWindow() const
	{
	    return m_water_window;
	}

	treal GetDFWaterWindow() const
	{
	    return m_df_water_window;
	}

	treal GetLipFilterFreq() const
    {
        return m_lipid_filt_freq;
    }

	void SetWaterWindow(treal width)
	{
	    m_water_window = width;
	}

	void SetDFWaterWindow(treal width)
	{
	    m_df_water_window = width;
	}

	treal GetRef() const
	{
	    return m_ref;
	}

	void SetRef(treal ref)
	{
	    m_ref = ref;
	}

	treal GetMaxDRef() const
	{
	    return m_max_dref;
	}
    
    void SetMaxDRef(treal max_d)
	{
	    m_max_dref = max_d;
	}
    
    treal GetNUMERICAL_TOL_INIT_MU() const
    {
        return m_nt_init_mu;
    }

    treal GetNUMERICAL_TOL_L2_DP() const
    {
        return m_nt_l2_dp;
    }

    treal GetPRESS_TE1() const
    {
        return m_press_TE1;
    }

    treal GetGnuplotCex() const
    {
        return m_gnuplot_cex;
    }

    treal GetThreads() const
    {
        return m_threads;
    }

    treal GetMinRange() const
    {
        return m_min_range;
    }

    size_t GetLCMBasisPts() const
    {
        return m_lcm_basis_pts;
    }

    void SetPRESS_TE1(treal te1)
    {
        m_press_TE1 = te1;
    }

    treal GetSTEAM_TM() const
    {
        return m_steam_TM;
    }

    void SetSTEAM_TM(treal tm)
    {
        m_steam_TM = tm;
    }

    treal GetCPMG_N() const
    {
        return m_cpmg_N;
    }

    void SetCPMG_N(int N)
    {
        m_cpmg_N = N;
    }

	std::string GetTitle() const
	{
	    return m_title;
	}

	std::string GetRefFile() const
	{
	    return m_strRefFile;
	}

	std::string GetAvListFile() const
	{
	    return m_strAvListFile;
	}

	treal GetWConc() const
	{
	    return m_w_conc;
	}

	treal GetWAtt() const
	{
	    return m_w_att;
	}

	treal GetAuNorm() const
	{
	    return m_au_norm;
	}

	void SetWAtt(treal w_att)
	{
	    m_w_att = w_att;
	}

	void SetWConc(treal w_conc)
	{
	    m_w_conc = w_conc;
	}

	bool GetRefSpec() const
	{
	    return m_ref_spec;
	}

	void SetRefSpec(bool refspec)
	{
	    m_ref_spec = refspec;
	}

	bool GetAutoPhase() const
	{
	    return m_bAutoPhase;
	}

	bool GetCombinePreproc() const
	{
	    return m_bCombinePreproc;
	}

	bool GetPreFitPhase() const
	{
	    return m_bPreFitPhase;
	}

	bool GetPreFitShift() const
	{
	    return m_bPreFitShift;
	}

	bool GetPreFitBl() const
	{
	    return m_bPreFitBl;
	}

	bool GetAppendLCMBasis() const
	{
	    return m_bAppendLCMBasis;
	}

	bool GetAppendNegRefBasis() const
	{
	    return m_bAppendNegRefBasis;
	}

	bool GetLipidFilter() const
	{
	    return m_bLipidFilter;
	}

    void SetLipidFilter(bool val)
	{
	    m_bLipidFilter = val;
	}

	bool GetFullEcho() const
	{
	    return m_bFullEcho;
	}

    void SetFullEcho(bool val)
	{
	    m_bFullEcho = val;
	}

	bool GetPdfExt() const
	{
	    return m_bPdfExt;
	}

    bool GetPdfStack() const
	{
	    return m_bPdfStack;
	}

	bool GetFastFit() const
	{
	    return m_ff;
	}
        
	bool GetExtCSVFit() const
    {
        return m_ext_csv_fit;
    }

	bool GetPreWsShift() const
	{
        return m_pre_ws_shift;
	}

	bool GetKeepPreWsShift() const
	{
        return m_keep_pre_ws_shift;
	}

	void SetKeepPreWsShift(bool val)
	{
        m_keep_pre_ws_shift = val;
	}

    void SetPdfExt(bool val)
	{
	    m_bPdfExt = val;
	}

    void SetPdfStack(bool val)
	{
	    m_bPdfStack = val;
	}

	bool GetSwapRowCol() const
	{
	    return m_bSwapRowCol;
	}

    void SetSwapRowCol(bool val)
	{
	    m_bSwapRowCol = val;
	}

	bool GetPause() const
	{
	    return m_pause;
	}

	bool GetAdaptSp() const
	{
	    return m_bAdaptSp;
	}

	bool GetPreHSVD() const
	{
	    return m_pre_hsvd;
	}

	void SetAutoPhase(bool val)
	{
		m_bAutoPhase = val;
	}

	void SetCombinePreproc(bool val)
	{
	    m_bCombinePreproc = val;
	}

	bool GetAutoReference() const
	{
	    return m_bAutoRef;
	}

	void SetAutoReference(bool val)
	{
		m_bAutoRef = val;
	}

	bool GetWaterEddy() const
	{
	    return m_bWaterEddy;
	}

	void SetWaterEddy(bool val)
	{
		m_bWaterEddy = val;
	}

	bool GetBasisComp() const
	{
	    return m_basis_comp;
	}

	bool GetPrintParas() const
	{
	    return m_bPrintParas;
	}
    
    bool GetNoFit() const
	{
	    return m_bNofit;
	}

    bool GetReadOnly() const
	{
	    return m_bReadOnly;
	}

    bool GetReadWriteOnly() const
	{
	    return m_bReadWriteOnly;
	}

	bool GetShowPreprocessed() const
	{
	    return m_bShowPreprocessed;
	}

	void SetShowPreprocessed(bool showpre)
	{
	    m_bShowPreprocessed = showpre;
	}

	integer GetConvWindowWidth() const
	{
	    return m_conv_window_width;
	}
    
    int GetFitRows() const
	{
	    return m_fit_rows;
	}

    int GetFitCols() const
	{
	    return m_fit_cols;
	}

    int GetFitSlices() const
	{
	    return m_fit_slices;
	}

    void SetFitRows(int fit_rows)
	{
	    m_fit_rows = fit_rows;
	}

    void SetFitCols(int fit_cols)
	{
	    m_fit_cols = fit_cols;
	}

    void SetFitSlices(int fit_slices)
	{
	    m_fit_slices = fit_slices;
	}

	void SetConvWindowWidth(integer nPts)
	{
	    m_conv_window_width = nPts;
	}

	std::string GetViewFile() const 
	{
		return m_strViewFile;
	}

	std::string GetPlotSigs() const 
	{
		return m_strPlotSigs;
	}

	std::string GetBasisSaveFile() const
	{
	    return m_strBasisSaveFile;
	}

	std::string GetBasisSaveFileCSV() const
	{
	    return m_strBasisSaveFileCSV;
	}

	std::string GetBasisSaveFileLCM() const
	{
	    return m_strBasisSaveFileLCM;
	}

	void SetBasisSaveFile(std::string strFile)
	{
	    m_strBasisSaveFile = strFile;
	}
    
    void SetBasisSaveFileLCM(std::string strFile)
	{
	    m_strBasisSaveFileLCM = strFile;
	}

	std::string GetOutRawFile() const
	{
	    return m_strOutRaw;
	}
	
    std::string GetOutRawFileV3() const
	{
	    return m_strOutRawV3;
	}

	std::string GetOutRawFileW() const
	{
	    return m_strOutRawW;
	}

	std::string GetOutRawFileWV3() const
	{
	    return m_strOutRawWV3;
	}
    
    std::string GetOutRawFileLcm() const
	{
	    return m_strOutRawLcm;
	}

	std::string GetOutRawFileLcmW() const
	{
	    return m_strOutRawLcmW;
	}

	std::string GetOutPreFile() const
	{
	    return m_strOutPre;
	}
	
	std::string GetOutPostFile() const
	{
	    return m_strOutPost;
	}

	GEOptions GetGEOptions() const
	{
	    return m_geOptions;
	}

	GEOptions& GetGEOptions()
	{
	    return m_geOptions;
	}
	
	void SetFitList(const coord_vec& fit_list)
	{
	    m_fit_list = fit_list;
	}

	const coord_vec& GetFitList() const
	{
	    return m_fit_list;
	}

	coord_vec& GetFitList()
	{
	    return m_fit_list;
	}

    bool GetReplaceFp()
	{
	    return m_replace_fp;
	}
    
    const int GetPrependPts() const
	{
	    return m_prepend_pts;
	}
    
    const int GetTruncatePts() const
	{
	    return m_truncate_pts;
	}

	//private:

	//! The file format of the incoming FID.
	fid_format_e m_format;

	//! A integer representation of the incoming FID for serialision purposes.
	int m_format_key; 
	
    //! Chemical shift referencing mode
    chem_shift_ref_e m_ref_signals;
    
    //! Chemical shift referencing mode for dynamic freq correction
    chem_shift_ref_e m_dyn_ref_signals;
    
    //! Internal basis set
    basis_set_e m_int_basis_set;
    
    //! dynamic averaging scheme (WS)
    dyn_av_e m_dyn_av;

    //! dynamic averaging scheme (W)
    dyn_av_e m_dyn_av_w;

    //! pre-averaging frequency correction
    bool m_dyn_freq_corr;
    
    //! pair-wise dynamic frequency correction
    bool m_pdfc;
    
    //! Only fit SVS data?
	bool m_svs_only;

    //! Pulse sequence
    pul_seq_e m_pul_seq;

	//! Name of the FID file we will be loading.
	std::string m_strFile;

	//! Name of the water reference FID file.
	std::string m_strFileWater;

	//! Name of the basis directory
	std::string m_strBasisPath;

	//! Name of the xml file to save
	std::string m_strOutputXMLPath;

	//! True if we are using precompiled basis XML files.
	bool m_bUsePrecompiled;

	//! Sample at which analysis begins.
	integer m_nStart;

	//! Time at which analysis begins.
	double m_nStartTime;
    
    //! Freq of PPM ref signal
    double m_ref_freq;

	//! Sample at which analysis ends.
	integer m_nEnd;

	//! The maximum number of iterations to perform.
	integer m_max_iters;
    
    //! The frequency of HSVD lipid filter in PPM
    treal m_lipid_filt_freq;
    
    //! factor to zero-fill kspace
    size_t m_zfill_kspace;

    //! filter kspace?
    bool m_filter_kspace;

	//! Width of time-domain convolution window.
	integer m_conv_window_width;

	//! Width of water removal window.
	treal m_water_window;

	//! Width of water removal window.
	treal m_df_water_window;

	//! Lower limit on value of phi0.
	treal m_phi0_lower;

	//! Upper limit on value of phi0.
	treal m_phi0_upper;

	//! Starting value of phi0.
	treal m_phi0_typ;

	//! Lower limit on value of phi1.
	treal m_phi1_lower;

	//! Upper limit on value of phi1.
	treal m_phi1_upper;

	//! Starting value of phi1.
	treal m_phi1_typ;

	//! Lower limit on value of metab alpha.
	treal m_alpha_metab_lower;

	//! Upper limit on value of metab alpha.
	treal m_alpha_metab_upper;

	//! Starting value of metab alpha.
	treal m_alpha_metab_typ;

	//! Lower limit on value of broad alpha.
	treal m_alpha_broad_lower;

	//! Upper limit on value of broad alpha.
	treal m_alpha_broad_upper;

	//! Starting value of broad alpha.
	treal m_alpha_broad_typ;
    
    //! Max value for beta
    treal m_max_beta;
    
    //! Optimiser beta scale factor
    treal m_beta_scale;

    //! Starting PPM value for plotting results.
	treal m_ppm_start;
    
    //! End PPM value for plotting results.
	treal m_ppm_end;
    
    //! Gaussian line broadening.
	treal m_lb;

    //! Reference signal line broadening.
	treal m_lb_ref;
	
    //! Zero filling factor
    size_t m_zero_fill;

    //! Baseline points to use
    size_t m_baseline;
    
    //! Initial beta value
	treal m_init_beta;

    //! Initial beta value actually used
	std::vector<treal> m_init_beta_used;
    
    //! Soft constraint value
	treal m_lambda;

    //! water attenuation (0.7 default)
    treal m_w_att;

    //! water attenuation (Concentration of NMR visable water)
    treal m_w_conc;

    //! should the SVS values be scaled when no water file is availalble?
    bool m_au_norm;

    //! maximum allowed shift from ref by auto_ref (ppm)
    treal m_max_dref;
    
    //! maximum metabolite shift allowed (ppm)
    treal m_max_metab_shift;

    //! maximum broad shift allowed (ppm)
    treal m_max_broad_shift;

    //! optim damping (smaller the value, the closest to steepest descent)
    treal m_nt_init_mu;

    //! optim stopping threshold for l2 norm of e
    treal m_nt_l2_dp;
    
    //! PRESS TE1 parameter
    treal m_press_TE1;
    
    //! gnuplot font expansion
    treal m_gnuplot_cex;

    //! how many threads to use?
    int m_threads;
    
    //! STEAM TM parameter
    treal m_steam_TM;
    
    //! Minimum nStart allowed when auto-detecting
    treal m_min_range;

    //! N for lcm basis output
    size_t m_lcm_basis_pts;

    //! Number of CPMG pulses
    int m_cpmg_N;
    
    //! title of pdf output
	std::string m_title;

	//! The name of the dangerplot file for the raw FID.
	std::string m_strRefFile;

	//! Should the FID be automatically phased.	
	bool m_bAutoPhase;

	//! Add up CSI spectra before auto-phasing?
    bool m_bCombinePreproc;
    
    //! Do a pre-fit fit?
    bool m_bPreFitPhase;
    bool m_bPreFitShift;
    bool m_bPreFitBl;

	bool m_bAppendLCMBasis;
	bool m_bAppendNegRefBasis;

	//! Filter out lipids using HSVD?
    bool m_bLipidFilter;

	//! Full echo data?
    bool m_bFullEcho;

	//! Extended pdf output
    bool m_bPdfExt;

    //SJW (In options.hpp): Added m_bPdfStack for stacked/waterfall plot of the individual metabolites
    bool m_bPdfStack;
    
    //! Fast fit mode?
    bool m_ff;

    //! Shift water to middle of spec before water removal?
    bool m_pre_ws_shift;
    
    //! Keep the above shift
    bool m_keep_pre_ws_shift;

	//! Swap CSI rows and columns?
    bool m_bSwapRowCol;

	//! Should we use automatic referencing.
	bool m_bAutoRef;

	//! Should the water file be used for eddy current correction.
	bool m_bWaterEddy;

	//! Should the view_fit plot be paused?
	bool m_pause;

	//! Should the basis groups be compressed (ie LCModel mode)
	bool m_basis_comp;
    
	//! Print FID paras to std out
    bool m_bPrintParas;
    
	//! Don't do the fitting part
    bool m_bNofit;
    
    //! Read data in only
    bool m_bReadOnly;
    
    bool m_bReadWriteOnly;
    
    //! Adaptive start point
	bool m_bAdaptSp;
	
    //! Pre-HSVD
    bool m_pre_hsvd;
    
    //! Revert to old limits on Phi1 para
    bool m_old_phase;
	
    //! Adaptive end point
    bool m_bAdaptEp;

	//! The constraints.
	std::vector<CConstraint> m_constraints;

	//! True if the user wishes to see the preprocessed FID before continuing.
	bool m_bShowPreprocessed;

	//! Estimated frequency (in PPM) of the middle of the fftshifted spectrum (0Hz).
	treal m_ref;

	//! Stores whether the reference has been specified from the command line.
	bool m_ref_spec;

	//! The file we are going to load to show.
	std::string m_strViewFile;

	//! The string containing additional metabolites to be plotted
	std::string m_strPlotSigs;

	//! The name of the basis file we are going to generate.
	std::string m_strBasisSaveFile;
	
	//! The name of the CSV basis file we are going to generate.
    std::string m_strBasisSaveFileCSV;
    
	//! The name of the LCM basis file we are going to generate.
    std::string m_strBasisSaveFileLCM;

	//! The name of the dangerplot file for the raw FID.
	std::string m_strOutRaw;
	std::string m_strOutRawV3;
	
	//! The name of the dangerplot file for the raw water ref FID.
    std::string m_strOutRawW;
    std::string m_strOutRawWV3;

	//! The name of the LCM file for the raw FID.
	std::string m_strOutRawLcm;
	
	//! The name of the LCM file for the raw water ref FID.
    std::string m_strOutRawLcmW;

	//! The name of the dangerplot file for the preprocess FID.
	std::string m_strOutPre;

	//! The name of the dangerplot file for the postprocess FID.
	std::string m_strOutPost;

	//! The name of the output image file
	std::string m_strFileOutImag;

	//! The name of the pdf output image file
	std::string m_strFileOutPdfImag;

	//! The name of the text output file
	std::string m_strFileOutTxt;

    //! Filename to output csv list of shifts
	std::string m_dyn_shift_file;

	//! The name of the csv output file
	std::string m_strFileOutCSV;
	
	//! The name of the csv output file for geometery 
    std::string m_strFileOutCSVGeom;
	
	//! The name of the csv fit output file
    std::string m_strFileOutCSVFit;
    
	//! The name of the csv fit output file mag
    std::string m_strFileOutCSVFitMag;
    
	//! The name of the aligned spectra output file
    std::string m_strFileOutCSVSpectraAligned;

	//! The name of the aligned spectra output file
    std::string m_strFileOutCSVSpectraAlignedMag;
    
    //! List of voxels to be fitted
    std::string m_strAvListFile;

	//! Options specific to GE files.
	GEOptions m_geOptions;

	//! List of voxels to be fit
	coord_vec m_fit_list;

    //! Number of rows to be fit
    int m_fit_rows;

    //! Number of columns to be fit
    int m_fit_cols;

    //! Number of slices to be fit
    int m_fit_slices;

    int m_prepend_pts;
    bool m_replace_fp;
    bool m_ext_csv_fit;
    int m_truncate_pts;
    int m_hsvd_comps;
    int m_max_hsvd_pts;
    bool m_crlb_td;
    bool m_crlb_optim;
    bool m_nnls;
    bool m_soft_cons;
    bool m_lineshape_corr;
    bool m_invert_even_pairs;

    };

}

#endif
