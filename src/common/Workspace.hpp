#ifndef TARQUIN_RESULTS_INCLUDED 
#define TARQUIN_RESULTS_INCLUDED 

#include "CFID.hpp"
#include "CBasis.hpp"
#include "Options.hpp"
#include <complex>

namespace tarquin 
{
	/*!
	 * Represet the options, final and intermediate results of the TARQUIN algorithm.
	 */
	class Workspace
	{
		public:

			//! Output to a stream for a useful summary display.
			friend std::ostream& operator << (std::ostream& os, const Workspace& rhs);

		public:

			//! See if all the necessary parameters have been specified.
			bool CheckParameters();
            
            //! Constructor
            Workspace()
            {
                 m_fidWater.SetWRef(true);
            }

			//! Return a reference to the FID we are fitting.
			CFID& GetFIDRaw()
			{
				return m_fidRaw;
			}

			/*! Return a reference to the raw FID or the processed FID, defaulting to the
			 * processed FID if both are avaialble. */
			CFID& GetFID()
			{
				if( m_fidProc.GetNumberOfPoints() )
					return m_fidProc;
				else
					return m_fidRaw;
			}

			const CFID& GetFID() const
			{
				if( m_fidProc.GetNumberOfPoints() )
					return m_fidProc;
				else
					return m_fidRaw;
			}

			const CFID& GetFIDRaw() const 
			{
				return m_fidRaw;
			}

			CFID& GetFIDProc()
			{
				return m_fidProc;
			}

			const CFID& GetFIDWater() const
			{
				return m_fidWater;
			}

			CFID& GetFIDWater()
			{
				return m_fidWater;
			}

			const CFID& GetFIDProc() const
			{
				return m_fidProc;
			}

			Options& GetOptions()
			{
				return m_options;
			}

			const Options& GetOptions() const
			{
				return m_options;
			}

			//! Return a reference to the modified group matrix object.
			cmat_stdvec& GetGroupMatrix()
			{
				return m_matGroups;
			}

			//! Return a reference to the modified basis matrix object.
			cmat_stdvec& GetBasisMatrix()
			{
				return m_matBasis;
			}

			const cmat_stdvec& GetBasisMatrix() const
			{
				return m_matBasis;
			}

			void SetBasisMatrix(cmat_stdvec basis_mat)
			{

				//m_matBasis.resize(basis_mat.msize(),basis_mat.nsize());
				m_matBasis = basis_mat;
			}

			void SetGroupMatrix(cmat_stdvec group_mat)
			{

				//m_matGroups.resize(group_mat.msize(),group_mat.nsize());
				m_matGroups = group_mat;
			}

			//! Get the basis object.
			CBasis& GetBasis()
			{
				return m_basis;
			}

			//! Get the basis object const version.
			const CBasis& GetBasis() const
			{
				return m_basis;
			}

			//! Get reference to amplitude vector so we can always store it here
			rvec_stdvec& GetAmplitudes()
			{
				return m_amplitudes;
			}

			void SetAmplitudes(const rvec_stdvec amps)
			{
				m_amplitudes = amps;
			}

            rvec_stdvec& GetAmplitudesComb()
			{
				return m_amplitudes_comb;
			}

			void SetAmplitudesComb(const rvec_stdvec amps)
			{
				m_amplitudes_comb = amps;
			}

            void AppendParas(const cvm::rvector paras)
            {
                m_paras.push_back(paras);
            }

            const cvm::rvector& GetParas(size_t n) const
            {
                return m_paras[n];
            }

            void AppendWaterWidth(const treal width)
            {
                m_water_width.push_back(width);
            }
            
            treal GetWaterWidth(size_t n) const
            {
                assert( n < m_water_width.size() );
                return m_water_width[n];
            }

            void AppendDynShift(const treal shift)
            {
                m_dyn_shift.push_back(shift);
            }

            void ClearDynShift()
            {
                m_dyn_shift.clear();
            }
            
            treal GetDynShift(size_t n) const
            {
                assert( n < m_dyn_shift.size() );
                return m_dyn_shift[n];
            }
            
            size_t GetDynShiftSize() const
            {
                return m_dyn_shift.size();
            }


            void AppendWaterFreq(const treal freq)
            {
                m_water_freq.push_back(freq);
            }
            
            treal GetWaterFreq(size_t n) const
            {
                assert( n < m_water_freq.size() );
                return m_water_freq[n];
            }

			const rvec_stdvec& GetAmplitudes() const
			{
				return m_amplitudes;
			}

			rvec_stdvec& GetCRLBs()
			{
				return m_crlbs;
			}

			const rvec_stdvec& GetCRLBs() const
			{
				return m_crlbs;
			}
			
            rvec_stdvec& GetCRLBsComb()
			{
				return m_crlbs_comb;
			}

			const rvec_stdvec& GetCRLBsComb() const
			{
				return m_crlbs_comb;
			}

			std::vector<std::string>& GetMetabNamesComb()
			{
				return m_metab_names_comb;
			}
			
            const std::vector<std::string>& GetMetabNamesComb() const
			{
				return m_metab_names_comb;
			}

			void SetCRLBs(rvec_stdvec crlbs)
			{
				m_crlbs = crlbs;
			}
			
            void SetCRLBsComb(rvec_stdvec crlbs)
			{
				m_crlbs_comb = crlbs;
			}

			const std::vector<double>& GetQ() const
			{
				return m_Q;
			}
            
            std::vector<double>& GetQ()
			{
				return m_Q;
			}

			const std::vector<double>& GetQ_rel() const
			{
				return m_Q_rel;
			}
            
            std::vector<double>& GetQ_rel()
			{
				return m_Q_rel;
			}

			const std::vector<double>& GetMetabRat() const
			{
				return m_metab_rat;
			}
            
            std::vector<double>& GetMetabRat()
			{
				return m_metab_rat;
			}

			const std::vector<double>& GetPeakMetabRat() const
			{
				return m_peak_metab_rat;
			}
            
            std::vector<double>& GetPeakMetabRat()
			{
				return m_peak_metab_rat;
			}

            std::vector<double>& GetMetabFWHM()
			{
				return m_metab_fwhm;
			}

            const std::vector<double>& GetMetabFWHM() const
			{
				return m_metab_fwhm;
			}

			void SetQ(std::vector<double> Q)
			{
				m_Q = Q;
			}

            const std::vector<double>& GetBLV() const
			{
				return m_blv;
			}
            
            std::vector<double>& GetBLV()
			{
				return m_blv;
			}

			void SetBLV(std::vector<double> blv)
			{
				m_blv = blv;
			}
            
            const std::vector<double>& GetBLS() const
			{
				return m_bls;
			}
            
            std::vector<double>& GetBLS()
			{
				return m_bls;
			}

			void SetBLS(std::vector<double> bls)
			{
				m_bls = bls;
			}
            
            const std::vector<double>& GetSpecNoise() const
			{
				return m_spec_noise;
			}
            
            std::vector<double>& GetSpecNoise()
			{
				return m_spec_noise;
			}

			void SetSpecNoise(std::vector<double> spec_noise)
			{
				m_spec_noise = spec_noise;
			}

            const std::vector<double>& GetMetabSNR() const
			{
				return m_metab_snr;
			}
            
            std::vector<double>& GetMetabSNR()
			{
				return m_metab_snr;
			}

			void SetMetabSNR(std::vector<double> metab_snr)
			{
				m_metab_snr = metab_snr;
			}


			/*!
			 * \brief Compute the normalised amplitude of the argument, i.e. scaled
			 * by water signal, any multiplication we did for stability, etc.
			 */

            /*
            Notes on water scaling:
            For 1H a signal simulated with one proton has an amplitude of 1/2
            */
			treal NormaliseValue(int voxel, double val) const
			{
				assert( 0 != m_norm_val[voxel]);

				// the unnormalised water amplitude = 
				// (the amplitude of the water signal * the value 
                // we normalised by for stability)
                assert(m_fidWater.GetNormValue() == 1.0);
				treal aw = m_norm_val[voxel] * m_fidWater.GetNormValue();

				// undo the normalisation we did for stability
                assert(m_fidProc.GetNormValue() == 1.0);
				val *= m_fidProc.GetNormValue();

                // divide the amplitude of the value by 2
                // to undo the simulator 1H scaling factor of 1/2
                // ie simulator outputs signal with amp=0.5
                // for molecule containing one proton
                val /= 2.0;
                
                // find out if a water file has been specified
                if ( m_options.GetFilenameWater().size() > 0 )
                {
                    // factor of 2.0 is for 2 protons in water
                    val *= m_options.GetWConc() * m_options.GetWAtt() * 2.0 / aw;
                }
                
				return val;
			}

			//! Get the normalised amplitudes.
			const rvec_stdvec& GetAmplitudesNormalised() const
			{
				return m_amp_norm;
			}

			//! Get the normalised amplitudes.
			rvec_stdvec& GetAmplitudesNormalised()
			{
				return m_amp_norm;
			}
			
			//! Set the normalised amplitudes.
			void SetAmplitudesNormalised(const rvec_stdvec& amp_norm)
			{
				//m_amp_norm.resize(amp_norm.size());
				m_amp_norm = amp_norm;
			}

			//! Get the normalised CRLBs.
			const rvec_stdvec& GetCRLBsNormalised() const
			{
				return m_crlbs_norm;
			}
	
			//! Get the normalised CRLBs.
			rvec_stdvec& GetCRLBsNormalised()
			{
				return m_crlbs_norm;
			}

			//! Set the normalised CRLBs.
			void SetCRLBsNormalised(const rvec_stdvec& crlbs_norm)
			{
				//m_crlbs_norm[0].resize(crlbs_norm.size());
				m_crlbs_norm = crlbs_norm;
			}

			//! Get the normalised amplitudes.
			const rvec_stdvec& GetAmplitudesNormalisedComb() const
			{
				return m_amp_norm_comb;
			}

			//! Get the normalised amplitudes.
			rvec_stdvec& GetAmplitudesNormalisedComb()
			{
				return m_amp_norm_comb;
			}
			
			//! Set the normalised amplitudes.
			void SetAmplitudesNormalisedComb(const rvec_stdvec& amp_norm)
			{
				//m_amp_norm.resize(amp_norm.size());
				m_amp_norm_comb = amp_norm;
			}

			//! Get the normalised CRLBs.
			const rvec_stdvec& GetCRLBsNormalisedComb() const
			{
				return m_crlbs_norm_comb;
			}
	
			//! Get the normalised CRLBs.
			rvec_stdvec& GetCRLBsNormalisedComb()
			{
				return m_crlbs_norm_comb;
			}

			//! Set the normalised CRLBs.
			void SetCRLBsNormalisedComb(const rvec_stdvec& crlbs_norm)
			{
				//m_crlbs_norm[0].resize(crlbs_norm.size());
				m_crlbs_norm_comb = crlbs_norm;
			}

			//! Get reference to signal estimate.
			cvec_stdvec& GetSignalEstimate()
			{
				return m_yhat;
			}

			void SetSignalEstimate(const cvec_stdvec& yhat)
			{
				//m_yhat.resize(yhat.size());
				m_yhat = yhat;
			}

			const cvec_stdvec& GetSignalEstimate() const
			{
				return m_yhat;
			}

			//! Get the options passed to the optimiser.
			std::vector<double>& GetLMopts()
			{
				return m_LM_opts;
			}

			//! Set the options passed to the optimiser.
			void SetLMopts( std::vector<double> LM_opts )
			{
				m_LM_opts = LM_opts;
			}

			//! Get the convergence results from the optimiser.
			std::vector<std::vector<double> >& GetLMinfo()
			{
				return m_LM_info;
			}

			const std::vector<std::vector<double> >& GetLMinfo() const
			{
				return m_LM_info;
			}

			//! Set the convergence results from the optimiser.	
			void SetLMinfo( std::vector<std::vector<double> > LM_info )
			{
				m_LM_info = LM_info;
			}

			//! Set the value we using when giving water-scaled outputs.
			void SetNormalisationValue(std::vector<double> norm_val)
			{
				m_norm_val = norm_val;
			}

			//! Get the value used for the water scaled outputs.
			const std::vector<double>& GetNormalisationValue() const
			{
				return m_norm_val;
			}

			std::vector<double>& GetNormalisationValue()
			{
				return m_norm_val;
			}

		protected:

			//! FID we are analysing (straight from the file).
			CFID m_fidRaw;

			//! Preprocessed FID (water removed, phased, referenced).
			CFID m_fidProc;

			//! Water reference FID (is not necessarily used).
			CFID m_fidWater;

			//! Basis set used for the analysis.
			CBasis m_basis;

			//! Options for the algorithm.
			Options m_options;

			//! Output of fit - modified group matrix.
			cmat_stdvec m_matGroups;

			//! Output of fit - modified metabolite matrix.
			cmat_stdvec m_matBasis;

			//! Output of fit - amplitudes of each metabolite, the raw numbers
			rvec_stdvec m_amplitudes;

			//! Output of fit - amplitudes of each metabolite, the raw numbers for combined signals (TNAA=NAA+NAAG)
			rvec_stdvec m_amplitudes_comb;
		    
            //! Paramters
            rvec_stdvec m_paras;

			//! Normalised amplitudes, by water or available, or by total amplitude if not.
			rvec_stdvec m_amp_norm;
			
            //! Normalised amplitudes, for combined signals
			rvec_stdvec m_amp_norm_comb;

			//! Estimate of CRLBS (one per metabolite).
			rvec_stdvec m_crlbs;
			
			//! Estimate of CRLBS (one per metabolite) for combined signals.
            rvec_stdvec m_crlbs_comb;

			//! Normalised CLRBS, by water or available, or by total amplitude if not.
			rvec_stdvec m_crlbs_norm;

			//! Normalised CLRBS, for combined signals
			rvec_stdvec m_crlbs_norm_comb;

			//! Estimate of signal.
			cvec_stdvec m_yhat;

			//! Normalisation value (divider of amplitudes).
			std::vector<double> m_norm_val;

			//! Width in Hz of water ref peak
			std::vector<double> m_water_width;
		    
            //! Shift in Hz of dynamic shift correction
            std::vector<double> m_dyn_shift;
            
            //! Frequency in Hz of water ref peak
			std::vector<double> m_water_freq;

			//! options passed to LM
			std::vector<double> m_LM_opts;

			//! convergence info from LM optimisation
			std::vector<std::vector<double> > m_LM_info;

			//! Fit quality
			std::vector<double> m_Q;

			std::vector<double> m_Q_rel;
			
            std::vector<double> m_metab_rat;
            
            std::vector<double> m_peak_metab_rat;
			
			//! Metabolite FWHM
            std::vector<double> m_metab_fwhm;
			
            //! Baseline variability
            std::vector<double> m_blv;
            
            //! Baseline shape
            std::vector<double> m_bls;
            
            //! Spectral noise
            std::vector<double> m_spec_noise;
            
            //! Metabolite SNR
            std::vector<double> m_metab_snr;
            
            //! Combined metabolite names found (TNAA, TCho...)
            std::vector<std::string> m_metab_names_comb;
	};
}

#endif
