#ifndef __CFID__
#define __CFID__

#include "common.hpp"
#include "cvm_util.hpp"
#include <stdint.h>
#include <fstream>

namespace tarquin 
{
/*
 * This wraps the FID, associated variables and flags to indicate whether variables have been
 * initialised to sensible values. The FID stored does not have the member phasing parameters applied to it.
 *
 * - m_cvmFID is the time domain FID signal
 *
 *   TODO: Consider moving the FID loading functions into a separate object.
 */

// 1-based structure for specifying coords
struct coord 
{
	int32_t row;
	int32_t col;
	int32_t slice;
	
	coord() : row(-1), col(-1), slice(-1) { }

	coord(int32_t row_, int32_t col_, int32_t slice_) :
		row(row_), col(col_), slice(slice_) { }
};

inline bool operator==(const tarquin::coord& lhs, const tarquin::coord& rhs)
{
    return (lhs.row == rhs.row) && (lhs.col == rhs.col) && (lhs.slice == rhs.slice);
}

inline bool operator<(const tarquin::coord& lhs, const tarquin::coord& rhs)
{
    if ( lhs.slice < rhs.slice )
        return true;
    else if (lhs.slice == rhs.slice )
    {
        if ( lhs.col < rhs.col )
            return true;
        else if (lhs.col == rhs.col )
        {
            if ( lhs.row < rhs.row )
                return true;
            else
                return false;
        }
        else
            return false;
    }
    else
        return false;
}


// type to respresent a list of fids to be analysed
typedef std::vector<coord> coord_vec;

// Reperesents a series of complex vectors eg a multi-vovel acquisition
typedef std::vector<cvm::cvector> cvec_stdvec;

// Reperesents a series of real vectors
typedef std::vector<cvm::rvector> rvec_stdvec;

// Reperesents a series of complex matrices
typedef std::vector<cvm::cmatrix> cmat_stdvec;

// Reperesents a series of real matrices
typedef std::vector<cvm::rmatrix> rmat_stdvec;

typedef std::vector<std::pair<treal, bool> > pair_vec;

class Options;
class CBoswell;
class Workspace;

class CFID 
{
	private:

		//! Name of file from which this was loaded.
		std::string m_strFilename;

		//! std vector of cvm vectors for the FID samples (unphased).
		cvec_stdvec m_cvmFID;

		//! Sampling frequency.
		std::pair<treal, bool> m_sampling_frequency;

		//! Transmitter frequency.
		std::pair<treal, bool> m_transmitter_frequency;

		//! Abitrary string specifying the sequence.
		std::pair<std::string, bool> m_pulse_sequence;

		//! Number of averages used to make the signal.
		std::pair<integer, bool> m_num_averages;

		//! Zero-order phase correction.
		pair_vec m_phi0;

		//! First-order phase correction.
		pair_vec m_phi1;

		//! Reference variable
		pair_vec m_ref;

		//! Echo time if this is a spin echo (time in seconds).
		std::pair<treal, bool> m_echo_time;	

		//! SNR estimate.
		pair_vec m_snr;
		
		//! Set number of points (e.g. for a dummy FID).
		integer m_nPoints;

		//! The value used to normalise this FID (i.e. give its vector infinity norm of 1 in FD)
		treal m_norm_val;
	    
        //! Does this data file contain water reference data also? TODO put into options?
        bool m_cwf;
        
        //! Am I water reference data? TODO add me to serialised data
        bool m_wref;

        //! Am I dynamic data? TODO add me to serialised data
        bool m_dyn;

		//! Rows of data
		integer m_rows;

		//! Columns of data
		integer m_cols;

		//! Slices of data
		integer m_slices;
        
        //! Geom info (rows, col, slice)
		std::pair<std::vector<double>, bool> m_voxel_dim;
        std::pair<std::vector<double>, bool> m_voi_dim;
        std::pair<std::vector<double>, bool> m_pos;
        std::pair<std::vector<double>, bool> m_row_dirn;
        std::pair<std::vector<double>, bool> m_col_dirn;

	public:

		CFID();

		~CFID();

		// TODO: Eliminate the need for this.
		friend bool ParseCommandLine(int argc, char* argv[], Options& options, CFID& fid);


		//! Set the sampling frequency parameters, etc. from the FID that is the argument.
		inline void SetParametersFromFID(const CFID& rhs)
		{
			SetSamplingFrequency(rhs.GetSamplingFrequency());
			SetTransmitterFrequency(rhs.GetTransmitterFrequency());
			SetEchoTime(rhs.GetEchoTime());
			SetPPMRef(rhs.GetPPMRef());
			SetPulseSequence(rhs.GetPulseSequence());
		}

		//! Load a FID, depending on the options in options.
		void Load(std::string strFilename, Options& options, Workspace& workspace, CBoswell& log);
		
        //! Load the GE water ref FID, depending on the options in options.
		void LoadW(std::string strFilename, Options& options, CBoswell& log);

		void AppendFromVector(const cvm::cvector& y)
		{
			m_cvmFID.push_back(y);
			m_nPoints = y.size(); //TODO
		}

		void ClearVector()
		{
			m_cvmFID.clear();
			m_nPoints = 0; //TODO
		}

		void ClearVectors()
		{
			m_cvmFID.clear();
			m_nPoints = 0; //TODO
            m_phi0.clear();
            m_phi1.clear();
            m_ref.clear();
            // this should be linked to the default options value somehow
			m_ref.push_back(std::make_pair(4.65, false)); 
			m_phi0.push_back(std::make_pair(0, false));
			m_phi1.push_back(std::make_pair(0, false));
		}

		//! Save as dangerplot format.
		bool SaveToFile(std::string strFilename, size_t num = 0);

		//bool SaveToBinFile(std::string filename);

		//! Save as LCModel format.
		bool SaveToFileLCM(std::string strFilename);
        
        //! shift the CSI grid in units of voxels
        //! ie x = 0.5 is a shift in the column direction of +1/2 a voxel
        void ShiftGrid(double row_shift, double col_shift, double slice_shift, CBoswell& log);
        
        void zfill_kspace(size_t factor, Options& options, CBoswell& log);

        void AverageData(Options& options);
        
        void FreqCorrData(const Options& options);

		//! Get the frequency scale in Hz, for plotting this FID's DFT [-fs/2,...,+fs/2).
		cvm::rvector GetFreqScale() const;
		cvm::rvector GetFreqScale(int zf) const;

        //! get the spectrum
        cvm::rvector GetSpec(int zf);

		//! Get the time scale in seconds
		cvm::rvector GetTimeScale() const;

		//! Get the frequency scale in ppm, for plotting this FID's DFT.
		cvm::rvector GetPPMScale(const coord& voxel, int zf=1) const;

        //! Get some nice values to help pick the fit region in the GUI
        cvm::rvector GetRawMap();

		//! Convert a value from PPM to Hz specific to this FID's transmitter frequency.
		inline treal ConvertPPMtoHz(treal ppm) const
		{
			assert( true == m_transmitter_frequency.second );
			return ((ppm * m_transmitter_frequency.first) / 1e6);
		}

		//! Convert a value from Hz to PPM, specific to this FID's transmitter frequency.
		inline treal ConvertHzToPPM(treal hz) const
		{
			assert( true == m_transmitter_frequency.second );
			return ((hz * 1e6) / m_transmitter_frequency.first);
		}

	private:

		//! Some scaling to make all FIDs, regardless of original scale, in the same scale.
		void Normalise();

	protected:

	public:
		

		inline void ShiftRef(treal new_ref, size_t voxel = 0)
		{
			assert( m_ref[voxel].second == true );
			assert( m_sampling_frequency.second == true );
			assert( m_transmitter_frequency.second == true );

			treal shift_ref = m_ref[voxel].first - new_ref;
			integer N = m_nPoints;

			treal dt = 1.0 / m_sampling_frequency.first;
			tcomplex j(0, 1);
			tcomplex omega;
			omega = 2.0 * M_PI * -shift_ref * m_transmitter_frequency.first / 1e6;

			// shift the FID by shift_ref 
			for(integer n = 0; n < N; n++)
			{
				treal t = n*dt;
				m_cvmFID[voxel](n+1) = m_cvmFID[voxel](n+1)*exp(j*omega*t);
			}

			m_ref[voxel] = std::make_pair(new_ref, true);
			//std::cout << "shifting basis by:" << omega << std::endl;
			//std::cout << "dt:" << dt << std::endl;
			//std::cout << "N:" << N << std::endl;
		}
        
        inline void ShiftRef(treal new_ref, const coord& voxel)
		{

	        int v = vox2ind(voxel);
			assert( m_ref[v].second == true );
			assert( m_sampling_frequency.second == true );
			assert( m_transmitter_frequency.second == true );

			treal shift_ref = m_ref[v].first - new_ref;
			integer N = m_nPoints;

			treal dt = 1.0 / m_sampling_frequency.first;
			tcomplex j(0, 1);
			tcomplex omega;
			omega = 2.0 * M_PI * -shift_ref * m_transmitter_frequency.first / 1e6;

			// shift the FID by shift_ref 
			for(integer n = 0; n < N; n++)
			{
				treal t = n*dt;
				m_cvmFID[v](n+1) = m_cvmFID[v](n+1)*exp(j*omega*t);
			}

			m_ref[v] = std::make_pair(new_ref, true);
			//std::cout << "shifting basis by:" << omega << std::endl;
			//std::cout << "dt:" << dt << std::endl;
			//std::cout << "N:" << N << std::endl;
		}


		inline void SetNumberOfPoints(integer n)
		{
			m_nPoints = n;
		}

		inline integer GetNumberOfPoints() const
		{
			// this is when the FID is intialised from a file
			if( 0 == m_nPoints && m_cvmFID.size() > 0 )
				return m_cvmFID[0].size(); //TODO
			else
				return m_nPoints;
		}

		inline treal GetEchoTime() const
		{
			return m_echo_time.first;
		}

		const inline pair_vec& GetSNR() const
		{
			return m_snr;
		}

        inline pair_vec& GetSNR()
		{
			return m_snr;
		}

		inline void SetSNR(pair_vec& snr)
		{
			m_snr = snr;
		}
        
		inline std::string GetPulseSequence() const
		{
			return m_pulse_sequence.first;
		}

		inline void SetPulseSequence(std::string str, bool known=true)
		{
			m_pulse_sequence = std::make_pair(str, known);
		}

		inline void SetSamplingFrequency(treal fs, bool known=true)
		{
			m_sampling_frequency = std::make_pair(fs, known);
		}

		inline treal GetSamplingFrequency() const
		{
			//assert( m_sampling_frequency.second );
			return m_sampling_frequency.first;
		}

		inline void SetTransmitterFrequency(treal ft, bool known=true)
		{
			m_transmitter_frequency = std::make_pair(ft, known);
		}

		inline treal GetTransmitterFrequency() const
		{
			return m_transmitter_frequency.first;
		}
    
        inline bool GetCWF() const
		{
			return m_cwf;
		}

        inline bool GetWRef() const
		{
			return m_wref;
		}

		inline void SetAverages(int nAverages, bool known=true)
		{
			m_num_averages = std::make_pair(nAverages, known);
		}

		inline int GetAverages() const
		{
			//assert( true == m_bAveragesKnown );
			return m_num_averages.first;
		}

		inline void SetFilename(std::string strFilename)
		{
			m_strFilename = strFilename;
		}

		inline std::string GetFilename() const
		{
			return m_strFilename;
		}

		//! return a reference to the nth cvm fid vector.
		//! zero based
		inline cvm::cvector& GetVectorFID(int32_t n)
		{
			assert( static_cast<size_t>(n) < m_cvmFID.size() );
			assert( n > -1 );
			return m_cvmFID[n];
		}

		inline cvm::cvector& GetVectorFID(const coord& voxel)
		{
			int n = vox2ind(voxel);
        	assert( static_cast<size_t>(n) < m_cvmFID.size() );
			assert( n > -1 );
			return m_cvmFID[n];
		}


		inline const cvm::cvector& GetVectorFID(const coord& voxel) const
		{
			int n = vox2ind(voxel);
        	assert( static_cast<size_t>(n) < m_cvmFID.size() );
			assert( n > -1 );
			return m_cvmFID[n];
		}

		//! return a reference to the nth cvm fid vector (const version).
		inline const cvm::cvector& GetVectorFID(int32_t n) const
		{
            assert( static_cast<size_t>(n) < m_cvmFID.size() );
			assert( n > -1 );
			return m_cvmFID[n];
		}
		
		//! Return a reference to the cvm FID vector.
		inline cvec_stdvec& GetVectorFID()
		{
			return m_cvmFID;
		}

		//! Return a reference to the cvm FID vector (const version).
		inline const cvec_stdvec& GetVectorFID() const
		{
			return m_cvmFID;
		}
        
		inline void SetPPMRef(const coord& voxel, treal ppm, bool known=true)
		{
			int n = vox2ind(voxel);
            assert( static_cast<size_t>(n) < m_ref.size() );
			m_ref[n] = std::make_pair(ppm,known);
		}
        
		inline void SetPPMRef(pair_vec ref)
		{
			m_ref = ref;
		}

		inline void SetPPMRef(double ref)
        {
            // iterate through the elements
            for (pair_vec::iterator it = m_ref.begin(); it != m_ref.end(); ++it)
                *it = std::make_pair(ref, true);
        }

		inline treal GetPPMRef(const coord& voxel) const
		{
			int n = vox2ind(voxel);
			return m_ref[n].first;
		}
		
		inline pair_vec GetPPMRef() const
		{
			return m_ref;
		}
        
        inline treal GetPPMRef(size_t n) const
		{
            assert ( n < m_ref.size() ); // TODO find out why this doesn't work!
			return m_ref[n].first;
		}

		pair_vec& GetPhi0()
		{
			return m_phi0;
		}

		pair_vec& GetPhi1()
		{
			return m_phi1;
		}
        
        inline treal GetPhi0(const coord& voxel) const
		{
			int n = vox2ind(voxel);
			return m_phi0[n].first;
		}

        inline treal GetPhi1(const coord& voxel) const
		{
			int n = vox2ind(voxel);
			return m_phi1[n].first;
		}
        
        inline double LimitPhi0(double phi0)
        {
            double div = phi0 / M_PI;
            if ( div > 1 )
            {
                phi0 = phi0 - 2*M_PI;
                LimitPhi0(phi0);
            }
            else if ( div <= -1 )
            {
                phi0 = phi0 + 2*M_PI;
                LimitPhi0(phi0);
            }
            // result looks good
            return phi0;
        }
        
        inline void SetPhi0(const coord& voxel, double phi0)
		{
            //std::cout << "Phi0 size :" << m_phi0.size() << std::endl << std::flush;
            //std::cout << "ph0 in  : " << 180*phi0/M_PI << std::endl;
            // make sure that pi => phi0 > -pi
            phi0 = LimitPhi0(phi0); 
            //std::cout << "ph0 out : " << 180*phi0/M_PI << std::endl;

			int n = vox2ind(voxel);
            assert( static_cast<size_t>(n) < m_phi0.size() );
            m_phi0[n].first = phi0;
            m_phi0[n].second = true;
		}

        inline void SetPhi1(const coord& voxel, double phi1)
		{
            //std::cout << "Phi1 size :" << m_phi1.size() << std::endl << std::flush;
			int n = vox2ind(voxel);
            assert( static_cast<size_t>(n) < m_phi1.size() );
            m_phi1[n].first = phi1;
            m_phi1[n].second = true;
		}


        const pair_vec& GetPhi0() const
		{
			return m_phi0;
		}

		const pair_vec& GetPhi1() const
		{
			return m_phi1;
		}

		inline void SetPhi0(pair_vec phi0)
		{
			m_phi0 = phi0;
		}

		inline void SetPhi1(pair_vec phi1)
		{
			m_phi1 = phi1;
		}

		inline void SetEchoTime(treal tau, bool known=true)
		{
			m_echo_time = std::make_pair(tau, known);
		}

		inline bool IsKnownSamplingFrequency() const
		{
			return m_sampling_frequency.second;
		}

		inline bool IsKnownTransmitterFrequency() const
		{
			return m_transmitter_frequency.second;
		}

		inline bool IsKnownPulseSequence() const
		{
			return m_pulse_sequence.second;
		}

		inline bool IsKnownAverages() const
		{
			return m_num_averages.second;
		}

		inline bool IsKnownPhi0() const
		{
			return m_phi0[0].second; // TODO Martin
		}

		inline bool IsKnownPhi1() const
		{
			return m_phi1[0].second; // TODO Martin
		}

		inline bool IsKnownReference() const
		{
			return m_ref[0].second; // TODO Martin
		}

		inline bool IsKnownEchoTime() const
		{
			return m_echo_time.second;
		}

		inline bool IsKnownSNR()
		{
			return m_snr[0].second; // TODO Martin
		}
		
		inline treal GetNormValue() const
		{
			return m_norm_val;
		}

		inline void SetNormValue(treal normv)
		{
			m_norm_val = normv;
		}
        
        inline void SetCWF(bool cwf)
		{
			m_cwf = cwf;
		}

        inline void SetWRef(bool wref)
		{
			m_wref = wref;
		}

        inline void SetDyn(bool dyn)
		{
			m_dyn = dyn;
		}

        inline void SetRows(int rows)
		{
			m_rows = rows;
		}

        inline void SetCols(int cols)
		{
			m_cols = cols;
		}

        inline void SetSlices(int slices)
		{
			m_slices = slices;
		}

		int GetRows() const
		{
			return m_rows;
		}

		int GetCols() const
		{
			return m_cols;
		}

		int GetSlices() const
		{
			return m_slices;
		}
		
		int GetVoxelCount() const
		{
			return m_rows * m_cols * m_slices;
		}

		// function for obtaing a vector index from a voxel location
		inline int vox2ind(const coord& voxel)
		{
			int n = voxel.row + m_rows * (voxel.col - 1) + m_cols * m_rows * (voxel.slice - 1) - 1;
			return n;
		}

		inline int vox2ind(const coord& voxel) const
		{
			int n = voxel.row + m_rows * (voxel.col - 1) + m_cols * m_rows * (voxel.slice - 1) - 1;
			return n;
		}

        void trim_echo()
        {
            std::cout << "trimming echo" << std::endl;
            // do some preprocessing so this data can be analysed as normal
            // assumes the max echo point is in the middle of the vector
            for ( size_t n = 0; n < m_cvmFID.size(); n++ )
            {
                cvm::cvector &fid = m_cvmFID[n];
                size_t zero_pt = static_cast<size_t>(std::ceil(fid.size()/2.0f));
                size_t cut_pts = 5;
                //std::cout << zero_pt << std::endl;
                // create a new fid
                cvm::cvector fid_temp(zero_pt - cut_pts);
                for ( size_t n = 1; n < (zero_pt + 1 - cut_pts); n++ )
                {
                    //fid_temp(n) = fid(zero_pt + n - 1);
                    if ( n == 1 )
                        fid_temp(n) = fid(zero_pt + n - 1) + fid(zero_pt + n - 1);
                    else if ( n == zero_pt )
                        fid_temp(n) = 0;
                    else
                        fid_temp(n) = fid(zero_pt + n - 1) + conj(fid(zero_pt - n ));
                }
                fid.resize(zero_pt - cut_pts);
                fid = fid_temp;
                SetNumberOfPoints(zero_pt - cut_pts);
            }
        }
        
        void swap_row_col()
        {
            //std::cout << std::endl << m_cols << std::endl;
            //std::cout << m_rows << std::endl;
            // create a new fid list
            cvec_stdvec new_list;
            for ( int col = 1; col < m_rows + 1; col++ )
            {
                for ( int row = 1; row < m_cols + 1; row++ )
                {
	                coord spec_num(col, row, 1); 
                    //std::cout << vox2ind(spec_num) << std::endl;
                    new_list.push_back(m_cvmFID[vox2ind(spec_num)]);
                }
            }
            m_cvmFID = new_list;
            int temp_rows = m_rows;
            m_rows = m_cols;
            m_cols = temp_rows;
        }

		inline bool IsKnownVoxelDim() const
		{
			return m_voxel_dim.second;
		}

		inline void SetVoxelDim(const std::vector<double>& voxel_dim, bool known=true)
		{
            assert (voxel_dim.size() == 3);
			m_voxel_dim = std::make_pair(voxel_dim, known);
		}

		inline const std::vector<double>& GetVoxelDim() const
        {
            return m_voxel_dim.first;
        }
        
        inline bool IsKnownVoiDim() const
		{
			return m_voi_dim.second;
		}

		inline void SetVoiDim(const std::vector<double>& voi_dim, bool known=true)
		{
            assert (voi_dim.size() == 3);
			m_voi_dim = std::make_pair(voi_dim, known);
		}

		inline const std::vector<double>& GetVoiDim() const
        {
            return m_voi_dim.first;
        }
        
        inline bool IsKnownPos() const
		{
			return m_pos.second;
		}

		inline void SetPos(const std::vector<double>& pos, bool known=true)
		{
            assert (pos.size() == 3);
			m_pos = std::make_pair(pos, known);
		}

		inline const std::vector<double>& GetPos() const
        {
            return m_pos.first;
        }
        
        inline bool IsKnownRowDirn() const
		{
			return m_row_dirn.second;
		}

		inline void SetRowDirn(const std::vector<double>& row_dirn, bool known=true)
		{
            assert (row_dirn.size() == 3);
			m_row_dirn = std::make_pair(row_dirn, known);
		}

		inline const std::vector<double>& GetRowDirn() const
        {
            return m_row_dirn.first;
        }

        inline bool IsKnownColDirn() const
		{
			return m_col_dirn.second;
		}

		inline void SetColDirn(const std::vector<double>& col_dirn, bool known=true)
		{
            assert (col_dirn.size() == 3);
			m_col_dirn = std::make_pair(col_dirn, known);
		}

		inline const std::vector<double>& GetColDirn() const
        {
            return m_col_dirn.first;
        }

        inline void cross(const cvm::rvector &a, const cvm::rvector &b, cvm::rvector &vec_out)
        {
            vec_out.resize(3);
            vec_out(1) = a(2)*b(3)-a(3)*b(2);
            vec_out(2) = a(3)*b(1)-a(1)*b(3);
            vec_out(3) = a(1)*b(2)-a(2)*b(1);

            // normalise at the end
            vec_out = vec_out/vec_out.norm2();
        }

};


}


#endif
