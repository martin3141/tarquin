#include "proto_buf.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include "tarquin_pb.hpp"
#include "CFID.hpp"

using namespace tarquin;
using namespace std;

bool save_results(string file_name, Workspace results)
{
    sln::results sln_results;
    
    // save fidraw
    CFID fidraw = results.GetFIDRaw(); // TODO make const?
    sln::fid* sln_fidraw = sln_results.mutable_fidraw();
    serialise_fid(sln_fidraw, fidraw);

    // save fidproc
    CFID fidproc = results.GetFIDProc(); // TODO make const?
    sln::fid* sln_fidproc = sln_results.mutable_fidproc();
    serialise_fid(sln_fidproc, fidproc);

    // save fidwater
    CFID fidwater = results.GetFIDWater(); // TODO make const?
    sln::fid* sln_fidwater = sln_results.mutable_fidwater();
    serialise_fid(sln_fidwater, fidwater);

    // save basis
    CBasis basis = results.GetBasis(); // TODO make const?
    sln::basis* sln_basis = sln_results.mutable_basis();
    serialise_basis(sln_basis, basis);

	// yhat
    const cvec_stdvec& yhat = results.GetSignalEstimate();
    for ( size_t m = 0; m < yhat.size(); m++ )
    {
		sln::comp_vec* vec = sln_results.add_yhat();
        for ( int n = 0; n < yhat[m].size(); n++ )
		{
            sln::comp_point* element = vec->add_comp_point();
			element->set_real(yhat[m](n+1).real());
			element->set_imag(yhat[m](n+1).imag());
		}
    }
   
    // options 
    Options options = results.GetOptions();
    sln::options* sln_options = sln_results.mutable_options();
    serialise_options(sln_options, options);

    // ahat
	const rvec_stdvec& ahat = results.GetAmplitudes();
    for ( size_t m = 0; m < ahat.size(); m++ )
	{
		sln::vec* vec = sln_results.add_amplitudes();
        for ( int n = 0; n < ahat[m].size(); n++ )
			vec->add_point(ahat[m](n+1));
	}

    // matgroups
    const cmat_stdvec& matgroups_vec = results.GetGroupMatrix();
	for ( size_t l = 0; l < matgroups_vec.size(); l++ )
	{
		sln::comp_mat* mat = sln_results.add_matgroups();
		for ( int m = 0; m < matgroups_vec[l].msize(); m++ )
		{
			sln::comp_vec* vec = mat->add_comp_vec();
			for ( int n = 0; n < matgroups_vec[l].nsize(); n++ )
			{
				sln::comp_point* element = vec->add_comp_point();
				element->set_real(matgroups_vec[l](m+1, n+1).real());
				element->set_imag(matgroups_vec[l](m+1, n+1).imag());
			}
		}
	}

    // matbasis
    const cmat_stdvec& matbasis_vec = results.GetBasisMatrix();
	for ( size_t l = 0; l < matbasis_vec.size(); l++ )
	{
		sln::comp_mat* mat = sln_results.add_matbasis();
		for ( int m = 0; m < matbasis_vec[l].msize(); m++ )
		{
			sln::comp_vec* vec = mat->add_comp_vec();
			for ( int n = 0; n < matbasis_vec[l].nsize(); n++ )
			{
				sln::comp_point* element = vec->add_comp_point();
				element->set_real(matbasis_vec[l](m+1, n+1).real());
				element->set_imag(matbasis_vec[l](m+1, n+1).imag());
			}
		}
	}

    // lmopts
    const std::vector<double>& lmopts = results.GetLMopts();
    for ( size_t n = 0; n < lmopts.size(); n++ )
        sln_results.add_lm_opts(lmopts[n]);
    
    // lminfo
    const std::vector<std::vector<double> >& lminfo_vec = results.GetLMinfo();
    for ( size_t m = 0; m < lminfo_vec.size(); m++ )
    {
		sln::vec* vec = sln_results.add_lm_info();
        for ( size_t n = 0; n < lminfo_vec[m].size(); n++ )
            vec->add_point(lminfo_vec[m][n]);
    }

    // norm_val
    const std::vector<double>& norm_vals = results.GetNormalisationValue();
    for ( size_t n = 0; n < norm_vals.size(); n++ )
        sln_results.add_norm_val(norm_vals[n]);

    // crlbs
    const rvec_stdvec& crlbs = results.GetCRLBs();
    for ( size_t m = 0; m < crlbs.size(); m++ )
	{
		sln::vec* vec = sln_results.add_crlbs();
        for ( int n = 0; n < crlbs[m].size(); n++ )
			vec->add_point(crlbs[m](n+1));
	}

    // fit quality 
    const std::vector<double>& Q = results.GetQ();
    for ( size_t n = 0; n < Q.size(); n++ )
        sln_results.add_q(Q[n]);
    
    // amp norm	
	const rvec_stdvec& amp_norm = results.GetAmplitudesNormalised();
    for ( size_t m = 0; m < amp_norm.size(); m++ )
	{
		sln::vec* vec = sln_results.add_amp_norm();
        for ( int n = 0; n < amp_norm[m].size(); n++ )
			vec->add_point(amp_norm[m](n+1));
	}

	// crlbs_norm
	const rvec_stdvec& crlbs_norm = results.GetCRLBsNormalised();
    for ( size_t m = 0; m < crlbs_norm.size(); m++ )
	{
		sln::vec* vec = sln_results.add_crlbs_norm();
        for ( int n = 0; n < crlbs_norm[m].size(); n++ )
			vec->add_point(crlbs_norm[m](n+1));
	}

    fstream output(file_name.c_str(), ios::out | ios::trunc | ios::binary);
    if (!sln_results.SerializeToOstream(&output)) {
        output.close();
        cerr << "Failed to write file." << endl;
        return false;
    }

    output.close();

    return true;
}

bool load_results(string file_name, Workspace& results)
{
    sln::results sln_results;
    fstream input(file_name.c_str(), ios::in | ios::binary);
    if (!sln_results.ParseFromIstream(&input)) {
        cerr << "Failed to read file." << endl;
        return false;
    }

    CFID& fidraw = results.GetFIDRaw();
    copy_fid(sln_results.fidraw(), fidraw);

    CFID& fidproc = results.GetFIDProc();
    copy_fid(sln_results.fidproc(), fidproc);

    CFID& fidwater = results.GetFIDWater();
    copy_fid(sln_results.fidwater(), fidwater);

    // load options
    Options& options = results.GetOptions();
    copy_options(sln_results.options(), options);

    // load basis
    CBasis& basis = results.GetBasis();
    copy_basis(sln_results.basis(), basis);

	// yhat
	cvec_stdvec yhat_vec;
    for ( int m = 0; m < sln_results.yhat_size(); m++ )
	{
		cvm::cvector cvm_fid(sln_results.yhat(m).comp_point_size());
		for ( int n = 0; n < sln_results.yhat(m).comp_point_size(); n++ )
		{
			cvm_fid(n+1) = tcomplex(sln_results.yhat(m).comp_point(n).real(), 
									sln_results.yhat(m).comp_point(n).imag());
		}
		yhat_vec.push_back(cvm_fid);
	}
    results.SetSignalEstimate(yhat_vec);

	// signal amplitudes
	rvec_stdvec cvm_amp_vec;
    for ( int m = 0; m < sln_results.amplitudes_size(); m++ )
    {
		cvm::rvector cvm_amp(sln_results.amplitudes(m).point_size());
		for ( int n = 0; n < sln_results.amplitudes(m).point_size(); n++ )
		{
			cvm_amp(n+1) = sln_results.amplitudes(m).point(n);
		}
		cvm_amp_vec.push_back(cvm_amp);
    }
    results.SetAmplitudes(cvm_amp_vec);

    // set matgroups
	cmat_stdvec cvm_matgroups_vec;
    for ( int l = 0; l < sln_results.matgroups_size(); l++ )
    {
		cvm::cmatrix cvm_matgroups(sln_results.matgroups(l).comp_vec_size(), sln_results.matgroups(l).comp_vec(0).comp_point_size());
        for ( int m = 0; m < sln_results.matgroups(l).comp_vec_size(); m++ )
        {
            for ( int n = 0; n < sln_results.matgroups(l).comp_vec(0).comp_point_size(); n++ )
            {
                cvm_matgroups(m+1, n+1) = tcomplex(sln_results.matgroups(l).comp_vec(m).comp_point(n).real(), sln_results.matgroups(l).comp_vec(m).comp_point(n).imag());
            }
        }
		cvm_matgroups_vec.push_back(cvm_matgroups);
    }
    results.SetGroupMatrix(cvm_matgroups_vec);

    // set matbasis
	cmat_stdvec cvm_matbasis_vec;
	for ( int l = 0; l < sln_results.matbasis_size(); l++ )
	{
		cvm::cmatrix cvm_matbasis(sln_results.matbasis(l).comp_vec_size(), sln_results.matbasis(l).comp_vec(0).comp_point_size());
		for ( int m = 0; m < sln_results.matbasis(l).comp_vec_size(); m++ )
		{
			for ( int n = 0; n < sln_results.matbasis(l).comp_vec(0).comp_point_size(); n++ )
			{
				cvm_matbasis(m+1, n+1) = tcomplex(sln_results.matbasis(l).comp_vec(m).comp_point(n).real(), sln_results.matbasis(l).comp_vec(m).comp_point(n).imag());
			}
		}
		cvm_matbasis_vec.push_back(cvm_matbasis);
	}
    results.SetBasisMatrix(cvm_matbasis_vec);

    std::vector<double> lmopts;
    for ( int n = 0; n < sln_results.lm_opts_size(); n++ )
        lmopts.push_back(sln_results.lm_opts(n));
    results.SetLMopts(lmopts);

    std::vector<std::vector<double> > lminfo_vec;
    for ( int m = 0; m < sln_results.lm_info_size(); m++ )
    {
        std::vector<double> lminfo;
        for ( int n = 0; n < sln_results.lm_info(m).point_size(); n++ )
            lminfo.push_back(sln_results.lm_info(m).point(n));

        lminfo_vec.push_back(lminfo);
    }
    results.SetLMinfo(lminfo_vec);

    // norm_val
    std::vector<double> norm_vals;
    for ( int n = 0; n < sln_results.norm_val_size(); n++ )
        norm_vals.push_back(sln_results.norm_val(n));

    results.SetNormalisationValue(norm_vals);

    // crlbs
	rvec_stdvec crlbs_vec;
    for ( int m = 0; m < sln_results.crlbs_size(); m++ )
    {
		cvm::rvector cvm_crlbs(sln_results.crlbs(m).point_size());
		for ( int n = 0; n < sln_results.crlbs(m).point_size(); n++ )
			cvm_crlbs(n+1) = sln_results.crlbs(m).point(n);

		crlbs_vec.push_back(cvm_crlbs);
    }
    results.SetCRLBs(crlbs_vec);

	// fit quality 
    std::vector<double> Q;
    for ( int n = 0; n < sln_results.q_size(); n++ )
        Q.push_back(sln_results.q(n));

    results.SetQ(Q);

	// amp_norm
	rvec_stdvec amp_norm_vec;
    for ( int m = 0; m < sln_results.amp_norm_size(); m++ )
    {
		cvm::rvector cvm_amp_norm(sln_results.amp_norm(m).point_size());
		for ( int n = 0; n < sln_results.amp_norm(m).point_size(); n++ )
			cvm_amp_norm(n+1) = sln_results.amp_norm(m).point(n);

		amp_norm_vec.push_back(cvm_amp_norm);
    }
    results.SetAmplitudesNormalised(amp_norm_vec);

	// crlbs_norm
	rvec_stdvec crlb_norm_vec;
    for ( int m = 0; m < sln_results.crlbs_norm_size(); m++ )
    {
		cvm::rvector cvm_crlb_norm(sln_results.crlbs_norm(m).point_size());
		for ( int n = 0; n < sln_results.crlbs_norm(m).point_size(); n++ )
		{
			cvm_crlb_norm(n+1) = sln_results.crlbs_norm(m).point(n);
		}
		crlb_norm_vec.push_back(cvm_crlb_norm);
    }
    results.SetCRLBsNormalised(crlb_norm_vec);

    return true;
}

bool save_basis(string file_name, CBasis basis)
{
    // create the basis object
    sln::basis sln_basis;

    // generate basis stl structs (only used for signal names)
    basis.GenerateSTLstruc_pb();
    
    sln_basis.set_strbasispath(basis.GetBasisPath());
    
    std::vector<std::string> sig_names = basis.GetSignalNamesSTL();
    for ( size_t n = 0; n < sig_names.size(); n++ )
    {
        sln_basis.add_vecsignalfiles(sig_names[n]);
    }

    std::vector<CSignal> signals = basis.GetSignals();
	typedef std::vector<CSignal>::iterator basis_iterator;
    // for each signal (basis vector) 
    //sln::fid::fidelement* element = sln_fid.add_stlfid();
	for( basis_iterator itB = signals.begin(); itB != signals.end(); itB++ ) 
    {
        // points to a vector of fids
        //sln::signal* signal = sln_basis.add_signals();
	    // for each fid in the current basis vector
        // add signal to the basis
        //sln::fid::fidelement* element = sln_fid->add_stlfid();
        sln::signal* sln_signal = sln_basis.add_signals();
	    for( CSignal::fid_iterator itF = itB->begin(); itF != itB->end(); itF++ ) 
        {
            serialise_fid(sln_signal->add_fid_vec(), *itF);

            //std::cout << itF->GetPPMRef() << std::endl;
	    }
	}

    std::vector<bool> broad_sig = basis.GetBroadSig();
    for ( size_t n = 0; n < broad_sig.size(); n++ )
    {
        sln_basis.add_broad_sig(broad_sig[n]);
    }


    fstream output(file_name.c_str(), ios::out | ios::trunc | ios::binary);
    if (!sln_basis.SerializeToOstream(&output)) {
        output.close();
        cerr << "Failed to write file." << endl;
        return false;
    }

    output.close();

    return true;
}


bool load_basis(string file_name, CBasis& basis)
{
    sln::basis sln_basis;
    fstream input(file_name.c_str(), ios::in | ios::binary);
    if (!sln_basis.ParseFromIstream(&input)) 
        return false;

    copy_basis(sln_basis, basis); 
    return true;
}


void copy_basis(sln::basis sln_basis, CBasis& basis)
{
    basis.SetBasisPath(sln_basis.strbasispath());

    std::vector<std::string> sig_files;
    for ( int n = 0; n < sln_basis.vecsignalfiles_size(); n++ )
    {
        sig_files.push_back(sln_basis.vecsignalfiles(n));
    }
    basis.SetSignalNamesSTL(sig_files);


    std::vector<CSignal> signals;
    for ( int n = 0; n < sln_basis.signals_size(); n++ )
    {
        std::vector<CFID> fids;
        for ( int p = 0; p < sln_basis.signals(n).fid_vec_size(); p++ )
        {
            sln::fid sln_fid;
            sln_fid = sln_basis.signals(n).fid_vec(p);
            CFID fid;
            // populate CFID from sln_fid 
            copy_fid(sln_fid, fid); 
            fids.push_back(fid);
            
            //cvm::cvector test = fid.GetVectorFID();

            //std::cout << fid.GetPPMRef() << std::endl;

        }
        CSignal signal;
        signal.SetFids(fids);
        signals.push_back(signal);
    }
    basis.SetSignals(signals);

    std::vector<bool> broad_vec;
    for ( int n = 0; n < sln_basis.broad_sig_size(); n++ )
    {
        broad_vec.push_back(sln_basis.broad_sig(n));
    }
    basis.SetBroadSig(broad_vec);

    basis.GenerateNonSTLstruc_pb();

}

//bool serialise_constraints(sln::constraints* sln_constraints, CConstrint& constraints)
//{
//    return true;
//}

//void copy_constraints(sln::constraints sln_constraints, CConstraint& constraints)

bool serialise_options(sln::options* sln_options, Options& options)
{
    switch ( options.GetFormat() )
    {
        case tarquin::DANGER:
            sln_options->set_format(sln::options::DANGER);
        case tarquin::SIEMENS:
            sln_options->set_format(sln::options::SIEMENS);
        case tarquin::DCM:
            sln_options->set_format(sln::options::DCM);
        case tarquin::PHILIPS:
            sln_options->set_format(sln::options::PHILIPS);
        case tarquin::PHILIPS_DCM:
            sln_options->set_format(sln::options::PHILIPS_DCM);
        case tarquin::GE:
            sln_options->set_format(sln::options::GE);
        case tarquin::RDA:
            sln_options->set_format(sln::options::RDA);
        case tarquin::LCM:
            sln_options->set_format(sln::options::LCM);
        case tarquin::SHF:
            sln_options->set_format(sln::options::SHF);
        case tarquin::VARIAN:
            sln_options->set_format(sln::options::VARIAN);
        case tarquin::BRUKER:
            sln_options->set_format(sln::options::BRUKER);
        case tarquin::JMRUI_TXT:
            sln_options->set_format(sln::options::JMRUI_TXT);
        case tarquin::NOTSET:
            sln_options->set_format(sln::options::NOTSET);
    }
    sln_options->set_strbasispath(options.GetBasisPath());
    sln_options->set_buseprecompiled(options.GetUsePrecompiled());
    sln_options->set_nstart(options.GetRangeStart());
    sln_options->set_nend(options.GetRangeEnd());
    sln_options->set_phi0_lower(options.GetLowerLimitPhi0());
    sln_options->set_phi0_upper(options.GetUpperLimitPhi0());
    sln_options->set_phi0_typ(options.GetTypicalPhi0());
    sln_options->set_phi1_lower(options.GetLowerLimitPhi1());
    sln_options->set_phi1_upper(options.GetUpperLimitPhi1());
    sln_options->set_phi1_typ(options.GetTypicalPhi1());
    
	std::vector<CConstraint> constraints = options.GetConstraints();
    for ( size_t n = 0; n < constraints.size(); n++ )
    {
        sln::constraints* sln_constraint = sln_options->add_cons();
        sln_constraint->set_minalpha(constraints[n].GetMinAlpha());
        sln_constraint->set_maxalpha(constraints[n].GetMaxAlpha());
        sln_constraint->set_typalpha(constraints[n].GetTypAlpha());
        
        sln_constraint->set_minbeta(constraints[n].GetMinBeta());
        sln_constraint->set_maxbeta(constraints[n].GetMaxBeta());
        sln_constraint->set_typbeta(constraints[n].GetTypBeta());

        sln_constraint->set_minshifthz(constraints[n].GetMinShiftHz());
        sln_constraint->set_maxshifthz(constraints[n].GetMaxShiftHz());
        sln_constraint->set_typshifthz(constraints[n].GetTypShiftHz());
    }

    sln_options->set_stroutputxmlpath(options.GetOutputXMLPath());
    sln_options->set_strfile(options.GetFilename());
    sln_options->set_strfilewater(options.GetFilenameWater());
    sln_options->set_strfileoutimag(options.GetFilenameImag());
    sln_options->set_strfileouttxt(options.GetFilenameTxt());
    sln_options->set_strfileoutcsv(options.GetFilenameCSV());
    sln_options->set_conv_window_width(options.GetConvWindowWidth());
    sln_options->set_water_window(options.GetWaterWindow());
    sln_options->set_bautophase(options.GetAutoPhase());
    sln_options->set_bautoref(options.GetAutoReference());
    sln_options->set_bshowpreprocessed(options.GetShowPreprocessed());
    sln_options->set_ref(options.GetRef());
//    sln_options->set_ref_spec(options.GetRefSpec());

	coord_vec& fit_list = options.GetFitList();
	for( coord_vec::const_iterator i = fit_list.begin(); i != fit_list.end(); ++i )	
	{
        sln::coord* sln_coord = sln_options->add_fit_list();
        sln_coord->set_row( (*i).row );
        sln_coord->set_col( (*i).col );
        sln_coord->set_slice( (*i).slice );
	}

    return true;
}

bool serialise_basis(sln::basis* sln_basis, CBasis& basis)
{
    // generate basis stl structs (only used for signal names)
    basis.GenerateSTLstruc_pb();
    
    sln_basis->set_strbasispath(basis.GetBasisPath());
    
    std::vector<std::string> sig_names = basis.GetSignalNamesSTL();
    for ( size_t n = 0; n < sig_names.size(); n++ )
    {
        sln_basis->add_vecsignalfiles(sig_names[n]);
    }

    std::vector<CSignal> signals = basis.GetSignals();
	typedef std::vector<CSignal>::iterator basis_iterator;
    // for each signal (basis vector) 
    //sln::fid::fidelement* element = sln_fid.add_stlfid();
	for( basis_iterator itB = signals.begin(); itB != signals.end(); itB++ ) 
    {
        // points to a vector of fids
        //sln::signal* signal = sln_basis.add_signals();
	    // for each fid in the current basis vector
        // add signal to the basis
        //sln::fid::fidelement* element = sln_fid->add_stlfid();
        sln::signal* sln_signal = sln_basis->add_signals();
	    for( CSignal::fid_iterator itF = itB->begin(); itF != itB->end(); itF++ ) 
        {
            serialise_fid(sln_signal->add_fid_vec(), *itF);

            //std::cout << itF->GetPPMRef() << std::endl;
	    }
	}

    std::vector<bool> broad_sig = basis.GetBroadSig();
    for ( size_t n = 0; n < broad_sig.size(); n++ )
    {
        sln_basis->add_broad_sig(broad_sig[n]);
    }

    return true;
}

bool serialise_fid(sln::fid* sln_fid, CFID fid)
{
    sln_fid->set_strfilename(fid.GetFilename());

    sln_fid->set_fs(fid.GetSamplingFrequency());
    sln_fid->set_bsamplingfrequencyknown(fid.IsKnownSamplingFrequency());

    sln_fid->set_ft(fid.GetTransmitterFrequency());
    sln_fid->set_btransmitterfrequencyknown(fid.IsKnownTransmitterFrequency());

    sln_fid->set_strsequence(fid.GetPulseSequence());
    sln_fid->set_pulsesequenceknown(fid.IsKnownPulseSequence());

    sln_fid->set_naverages(fid.GetAverages());
    sln_fid->set_baveragesknown(fid.IsKnownAverages());
    
	const pair_vec& phi0 = fid.GetPhi0();
    for ( size_t n = 0; n < phi0.size(); n++ )
    {
        sln_fid->add_phi0(phi0[n].first);
        sln_fid->add_bzeroorderphaseknown(phi0[n].second);
    }
	
    const pair_vec& phi1 = fid.GetPhi1();
    for ( size_t n = 0; n < phi1.size(); n++ )
    {
        sln_fid->add_phi1(phi1[n].first);
        sln_fid->add_bfirstorderphaseknown(phi1[n].second);
    }

	const pair_vec& ref = fid.GetPPMRef();
    for ( size_t n = 0; n < ref.size(); n++ )
    {
        sln_fid->add_ref(ref[n].first);
        sln_fid->add_breferenceknown(ref[n].second);
    }
    
    sln_fid->set_tau(fid.GetEchoTime());
    sln_fid->set_bechoknown(fid.IsKnownEchoTime());
    

	const pair_vec& snr = fid.GetSNR();
    for ( size_t n = 0; n < snr.size(); n++ )
    {
        sln_fid->add_snr(snr[n].first);
        sln_fid->add_bsnrknown(snr[n].second);
    }

    sln_fid->set_npoints(fid.GetNumberOfPoints());

    sln_fid->set_norm_val(fid.GetNormValue());

    const cvec_stdvec cvm_fids = fid.GetVectorFID();

	for ( size_t n = 0; n < cvm_fids.size(); n++ )
	{
		sln::comp_vec* vec = sln_fid->add_fids();
		for ( int m = 0; m < cvm_fids[n].size(); m++ )
		{
			sln::comp_point* element = vec->add_comp_point();
			element->set_real(cvm_fids[n](m+1).real());
			element->set_imag(cvm_fids[n](m+1).imag());
		}
	}

    sln_fid->set_rows(fid.GetRows());
    sln_fid->set_cols(fid.GetCols());
    sln_fid->set_slices(fid.GetSlices());
    
    sln_fid->set_voxel_dim_known( fid.IsKnownVoxelDim() );
    const std::vector<double>& voxel_dim = fid.GetVoxelDim();
    for ( size_t n = 0; n < voxel_dim.size(); n++ )
        sln_fid->add_voxel_dim(voxel_dim[n]);

    sln_fid->set_voi_dim_known( fid.IsKnownVoiDim() );
    const std::vector<double>& voi_dim = fid.GetVoiDim();
    for ( size_t n = 0; n < voi_dim.size(); n++ )
        sln_fid->add_voi_dim(voi_dim[n]);

    sln_fid->set_pos_known( fid.IsKnownPos() );
    const std::vector<double>& pos = fid.GetPos();
    for ( size_t n = 0; n < pos.size(); n++ )
        sln_fid->add_pos(pos[n]);
    
    sln_fid->set_row_dirn_known( fid.IsKnownRowDirn() );
    const std::vector<double>& row_dirn = fid.GetRowDirn();
    for ( size_t n = 0; n < row_dirn.size(); n++ )
        sln_fid->add_row_dirn(row_dirn[n]);

    sln_fid->set_col_dirn_known( fid.IsKnownColDirn() );
    const std::vector<double>& col_dirn = fid.GetColDirn();
    for ( size_t n = 0; n < col_dirn.size(); n++ )
        sln_fid->add_col_dirn(col_dirn[n]);

    return true;
}
bool serialise_fid(sln::fid* sln_fid, CSignal::fid_iterator& fid)
{
    sln_fid->set_strfilename(fid->GetFilename());

    sln_fid->set_fs(fid->GetSamplingFrequency());
    sln_fid->set_bsamplingfrequencyknown(fid->IsKnownSamplingFrequency());

    sln_fid->set_ft(fid->GetTransmitterFrequency());
    sln_fid->set_btransmitterfrequencyknown(fid->IsKnownTransmitterFrequency());

    sln_fid->set_strsequence(fid->GetPulseSequence());
    sln_fid->set_pulsesequenceknown(fid->IsKnownPulseSequence());

    sln_fid->set_naverages(fid->GetAverages());
    sln_fid->set_baveragesknown(fid->IsKnownAverages());
    
    const pair_vec& phi0 = fid->GetPhi0();
    for ( size_t n = 0; n < phi0.size(); n++ )
    {
        sln_fid->add_phi0(phi0[n].first);
        sln_fid->add_bzeroorderphaseknown(phi0[n].second);
    }
	
    const pair_vec& phi1 = fid->GetPhi1();
    for ( size_t n = 0; n < phi1.size(); n++ )
    {
        sln_fid->add_phi1(phi1[n].first);
        sln_fid->add_bfirstorderphaseknown(phi1[n].second);
    }

	const pair_vec& ref = fid->GetPPMRef();
	for ( size_t n = 0; n < ref.size(); n++ )
    {
        sln_fid->add_ref(ref[n].first);
        sln_fid->add_breferenceknown(ref[n].second);
    }

    sln_fid->set_tau(fid->GetEchoTime());
    sln_fid->set_bechoknown(fid->IsKnownEchoTime());
    
    const pair_vec& snr = fid->GetSNR();
    for ( size_t n = 0; n < snr.size(); n++ )
    {
        sln_fid->add_snr(snr[n].first);
        sln_fid->add_bsnrknown(snr[n].second);
    }

    sln_fid->set_npoints(fid->GetNumberOfPoints());

    sln_fid->set_norm_val(fid->GetNormValue());
	
	const cvec_stdvec cvm_fids = fid->GetVectorFID();

	for ( size_t n = 0; n < cvm_fids.size(); n++ )
	{
		sln::comp_vec* vec = sln_fid->add_fids();
		for ( int m = 0; m < cvm_fids[n].size(); m++ )
		{
			sln::comp_point* element = vec->add_comp_point();
			element->set_real(cvm_fids[n](m+1).real());
			element->set_imag(cvm_fids[n](m+1).imag());
		}
	}
    
    sln_fid->set_rows(fid->GetRows());
    sln_fid->set_cols(fid->GetCols());
    sln_fid->set_slices(fid->GetSlices());
    
    sln_fid->set_voxel_dim_known( fid->IsKnownVoxelDim() );
    const std::vector<double>& voxel_dim = fid->GetVoxelDim();
    for ( size_t n = 0; n < voxel_dim.size(); n++ )
        sln_fid->add_voxel_dim(voxel_dim[n]);

    sln_fid->set_voi_dim_known( fid->IsKnownVoiDim() );
    const std::vector<double>& voi_dim = fid->GetVoiDim();
    for ( size_t n = 0; n < voi_dim.size(); n++ )
        sln_fid->add_voi_dim(voi_dim[n]);

    sln_fid->set_pos_known( fid->IsKnownPos() );
    const std::vector<double>& pos = fid->GetPos();
    for ( size_t n = 0; n < pos.size(); n++ )
        sln_fid->add_pos(pos[n]);
    
    sln_fid->set_row_dirn_known( fid->IsKnownRowDirn() );
    const std::vector<double>& row_dirn = fid->GetRowDirn();
    for ( size_t n = 0; n < row_dirn.size(); n++ )
        sln_fid->add_row_dirn(row_dirn[n]);

    sln_fid->set_col_dirn_known( fid->IsKnownColDirn() );
    const std::vector<double>& col_dirn = fid->GetColDirn();
    for ( size_t n = 0; n < col_dirn.size(); n++ )
        sln_fid->add_col_dirn(col_dirn[n]);

    return true;
}

bool save_fid(string file_name, CFID fid)
{
    sln::fid sln_fid;
    sln_fid.set_strfilename(fid.GetFilename());

    sln_fid.set_fs(fid.GetSamplingFrequency());
    sln_fid.set_bsamplingfrequencyknown(fid.IsKnownSamplingFrequency());

    sln_fid.set_ft(fid.GetTransmitterFrequency());
    sln_fid.set_btransmitterfrequencyknown(fid.IsKnownTransmitterFrequency());

    sln_fid.set_strsequence(fid.GetPulseSequence());
    sln_fid.set_pulsesequenceknown(fid.IsKnownPulseSequence());

    sln_fid.set_naverages(fid.GetAverages());
    sln_fid.set_baveragesknown(fid.IsKnownAverages());
    
    const pair_vec& phi0 = fid.GetPhi0();
    for ( size_t n = 0; n < phi0.size(); n++ )
    {
        sln_fid.add_phi0(phi0[n].first);
        sln_fid.add_bzeroorderphaseknown(phi0[n].second);
    }
	
    const pair_vec& phi1 = fid.GetPhi1();
    for ( size_t n = 0; n < phi1.size(); n++ )
    {
        sln_fid.add_phi1(phi1[n].first);
        sln_fid.add_bfirstorderphaseknown(phi1[n].second);
    }
	
    const pair_vec& ref = fid.GetPPMRef();
	for ( size_t n = 0; n < ref.size(); n++ )
    {
        sln_fid.add_ref(ref[n].first);
        sln_fid.add_breferenceknown(ref[n].second);
    }

    sln_fid.set_tau(fid.GetEchoTime());
    sln_fid.set_bechoknown(fid.IsKnownEchoTime());
    
	const pair_vec& snr = fid.GetSNR();
    for ( size_t n = 0; n < snr.size(); n++ )
    {
        sln_fid.add_snr(snr[n].first);
        sln_fid.add_bsnrknown(snr[n].second);
    }

    sln_fid.set_npoints(fid.GetNumberOfPoints());

    sln_fid.set_norm_val(fid.GetNormValue());

	const cvec_stdvec cvm_fids = fid.GetVectorFID();

	for ( size_t n = 0; n < cvm_fids.size(); n++ )
	{
		sln::comp_vec* vec = sln_fid.add_fids();
		for ( int m = 0; m < cvm_fids[n].size(); m++ )
		{
			sln::comp_point* element = vec->add_comp_point();
			element->set_real(cvm_fids[n](m+1).real());
			element->set_imag(cvm_fids[n](m+1).imag());
		}
	}

    // Write the fid to disk.
    fstream output(file_name.c_str(), ios::out | ios::trunc | ios::binary);
    if (!sln_fid.SerializeToOstream(&output)) {
        output.close();
        cerr << "Failed to write file." << endl;
        return false;
    }

    output.close();
    return true;
}

bool load_fid(string file_name, CFID& fid)
{
    sln::fid sln_fid;
    fstream input(file_name.c_str(), ios::in | ios::binary);
    if (!sln_fid.ParseFromIstream(&input)) {
        cerr << "Failed to read file." << endl;
        return false;
    }
    
    copy_fid(sln_fid, fid); 

    return true;
}

void copy_options(sln::options sln_options, Options& options)
{
   //sln_options->set_format(options.Get());
   //
   switch ( sln_options.format() )
    {
        case sln::options::DANGER:
            options.SetFormat(tarquin::DANGER);
        case sln::options::SIEMENS:
            options.SetFormat(tarquin::SIEMENS);
        case sln::options::DCM:
            options.SetFormat(tarquin::DCM);
        case sln::options::PHILIPS:
            options.SetFormat(tarquin::PHILIPS);
        case sln::options::PHILIPS_DCM:
            options.SetFormat(tarquin::PHILIPS_DCM);
        case sln::options::GE:
            options.SetFormat(tarquin::GE);
        case sln::options::RDA:
            options.SetFormat(tarquin::RDA);
        case sln::options::LCM:
            options.SetFormat(tarquin::LCM);
        case sln::options::SHF:
            options.SetFormat(tarquin::SHF);
        case sln::options::VARIAN:
            options.SetFormat(tarquin::VARIAN);
        case sln::options::BRUKER:
            options.SetFormat(tarquin::BRUKER);
        case tarquin::JMRUI_TXT:
            options.SetFormat(tarquin::JMRUI_TXT);
        case sln::options::NOTSET:
            options.SetFormat(tarquin::NOTSET);
    }

    options.SetBasisPath(sln_options.strbasispath());
    options.SetUsePrecompiled(sln_options.buseprecompiled());
    options.SetRangeStart(sln_options.nstart());
    options.SetRangeEnd(sln_options.nend());
    options.SetZeroPhaseLimit(sln_options.phi0_lower(), sln_options.phi0_upper());
    options.SetZeroPhaseTyp(sln_options.phi0_typ());
    options.SetFirstPhaseLimit(sln_options.phi1_lower(), sln_options.phi1_upper());
    options.SetFirstPhaseTyp(sln_options.phi1_typ());
    
    std::vector<CConstraint>& constraints = options.GetConstraints();

    for ( int n = 0; n < sln_options.cons_size(); n++ )
    {

        CConstraint constraint;
        constraint.SetMinAlpha(sln_options.cons(n).minalpha());
        constraint.SetMaxAlpha(sln_options.cons(n).maxalpha());
        constraint.SetTypAlpha(sln_options.cons(n).typalpha());

        constraint.SetMinBeta(sln_options.cons(n).minbeta());
        constraint.SetMaxBeta(sln_options.cons(n).maxbeta());
        constraint.SetTypBeta(sln_options.cons(n).typbeta());

        constraint.SetMinShiftHz(sln_options.cons(n).minshifthz());
        constraint.SetMaxShiftHz(sln_options.cons(n).maxshifthz());
        constraint.SetTypShiftHz(sln_options.cons(n).typshifthz());

        constraints.push_back(constraint);
    }

    options.SetOutputXMLPath(sln_options.stroutputxmlpath());  
    options.SetFilename(sln_options.strfile());  
    options.SetFilenameWater(sln_options.strfilewater());  
    options.SetFilenameImag(sln_options.strfileoutimag());  
    options.SetFilenameTxt(sln_options.strfileouttxt());  
    options.SetFilenameCSV(sln_options.strfileoutcsv());  
    options.SetConvWindowWidth(sln_options.conv_window_width());  
    options.SetWaterWindow(sln_options.water_window());  
    options.SetAutoPhase(sln_options.bautophase());  
    options.SetAutoReference(sln_options.bautoref());  
    options.SetShowPreprocessed(sln_options.bshowpreprocessed());  
    options.SetRef(sln_options.ref());  

	coord_vec& fit_list = options.GetFitList();
    for ( int n = 0; n < sln_options.fit_list_size(); n++ )
	{
        coord fit_coord(sln_options.fit_list(n).row(), 
						sln_options.fit_list(n).col(), 
						sln_options.fit_list(n).slice());

        fit_list.push_back(fit_coord);
	}

}

void copy_fid(sln::fid sln_fid, CFID& fid)
{
    fid.SetFilename(sln_fid.strfilename());

    fid.SetSamplingFrequency(sln_fid.fs(), sln_fid.bsamplingfrequencyknown());

    fid.SetTransmitterFrequency(sln_fid.ft(), sln_fid.btransmitterfrequencyknown());

    fid.SetPulseSequence(sln_fid.strsequence(), sln_fid.pulsesequenceknown());

    fid.SetAverages(sln_fid.naverages(), sln_fid.baveragesknown());

    pair_vec phi0;
	for ( int n = 0; n < sln_fid.phi0_size(); n++ )	
        phi0.push_back(std::make_pair(sln_fid.phi0(n), sln_fid.bzeroorderphaseknown(n)));
    
    fid.SetPhi0(phi0);

    pair_vec phi1;
	for ( int n = 0; n < sln_fid.phi1_size(); n++ )	
        phi1.push_back(std::make_pair(sln_fid.phi1(n), sln_fid.bfirstorderphaseknown(n)));
    
    fid.SetPhi1(phi1);

	pair_vec ref;
	for ( int n = 0; n < sln_fid.ref_size(); n++ )	
        ref.push_back(std::make_pair(sln_fid.ref(n), sln_fid.breferenceknown(n)));

    fid.SetPPMRef(ref);

    fid.SetEchoTime(sln_fid.tau(), sln_fid.bechoknown());
    
    pair_vec snr;
	for ( int n = 0; n < sln_fid.snr_size(); n++ )	
        snr.push_back(std::make_pair(sln_fid.snr(n), sln_fid.bsnrknown(n)));

    fid.SetSNR(snr);

    fid.SetNumberOfPoints(sln_fid.npoints());

    fid.SetNormValue(sln_fid.norm_val());

	for ( int m = 0; m < sln_fid.fids_size(); m++ )	
	{
		int len = sln_fid.fids(m).comp_point_size();
		cvm::cvector cvm_fid(len);
		for ( int n = 0; n < len; n++)
		{
			cvm_fid(n+1) = tcomplex(sln_fid.fids(m).comp_point(n).real(), sln_fid.fids(m).comp_point(n).imag());
		}
		fid.AppendFromVector(cvm_fid);
	}

    fid.SetRows(sln_fid.rows());
    fid.SetCols(sln_fid.cols());
    fid.SetSlices(sln_fid.slices());

    std::vector<double> vox_dim;
	for ( int n = 0; n < sln_fid.voxel_dim_size(); n++ )	
        vox_dim.push_back(sln_fid.voxel_dim(n));
    fid.SetVoxelDim(vox_dim, sln_fid.voxel_dim_known());
    
    std::vector<double> voi_dim;
	for ( int n = 0; n < sln_fid.voi_dim_size(); n++ )	
        voi_dim.push_back(sln_fid.voi_dim(n));
    fid.SetVoiDim(voi_dim, sln_fid.voi_dim_known());

    std::vector<double> pos;
	for ( int n = 0; n < sln_fid.pos_size(); n++ )	
        pos.push_back(sln_fid.pos(n));
    fid.SetPos(pos, sln_fid.pos_known());

    std::vector<double> row_dirn;
	for ( int n = 0; n < sln_fid.row_dirn_size(); n++ )	
        row_dirn.push_back(sln_fid.row_dirn(n));
    fid.SetRowDirn(row_dirn, sln_fid.row_dirn_known());

    std::vector<double> col_dirn;
	for ( int n = 0; n < sln_fid.col_dirn_size(); n++ )	
        col_dirn.push_back(sln_fid.col_dirn(n));
    fid.SetColDirn(col_dirn, sln_fid.col_dirn_known());
}
