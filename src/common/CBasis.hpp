#ifndef __CBASIS__
#define __CBASIS__

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <map>
#include <vector>
#include "common.hpp"
#include "CFID.hpp"
#include "CSignal.hpp"
#include "CBoswell.hpp"

namespace tarquin 
{
	/*!
	 * A class that defines the basis used to analyse a signal. 
	 *
	 * This contains a collection of CSignal objects, where each one corresponds to an analysis vector, e.g.
	 * a metabolite.
	 */

	class CBasis 
	{
		public:

			typedef std::vector<CSignal>::iterator basis_iterator;
			typedef std::vector<CSignal>::const_iterator const_basis_iterator;

			void Initialise(float ppm_ref);

			/*
			 * ! Generate a basis, using the files specified, to match the specified FID.
			 * \param strBasisPath is the directory containing the CSV files.
			 * \param fidMatch is the FID we are generating a basis for.
			 */ 

			void GenerateSTLstruc_pb();

			void GenerateNonSTLstruc_pb();

			bool Simulate(const CFID& fidMatch, const Options& opts, CBoswell& log);
			
			bool Simulate(std::string strBasisPath, const CFID& fidMatch, const Options& opts, CBoswell& log);

			bool check(const CFID& fid, CBoswell& log);

			//bool Simulate(std::string strBasisPath, const CFID& fidMatch);

			bool SimulateCSV(std::string strBasisPath, const Options& options, const CFID& fidMatch);

            bool ReadLCMBasis(std::string strBasisPath, const CFID& fidMatch, const Options& options, CBoswell& log);

			//! Return the signal structure 
			const std::vector<CSignal>& GetSignals() const
			{
				return m_signals;
			}

			void SetSignals(const std::vector<CSignal>& signals)
			{
				m_signals = signals;
			}

			//! Return the (summed) matrix of basis vectors.
			const cvm::cmatrix& GetBasisMatrix() const
			{
				assert( m_matBasis.msize() > 0 );
				return m_matBasis;
			}

			treal GetSamplingFrequency()
			{
				CFID& fid = *(m_signals[0].begin()); 
				return fid.GetSamplingFrequency();
			}

			treal GetTransmitterFrequency()
			{
				CFID& fid = *(m_signals[0].begin()); 
				return fid.GetTransmitterFrequency();
			}

			treal GetEchoTime()
			{
				CFID& fid = *(m_signals[0].begin()); 
				return fid.GetEchoTime();
			}

			integer GetNumberOfPoints()
			{
				CFID& fid = *(m_signals[0].begin()); 
				return fid.GetNumberOfPoints();
			}

			std::string GetPulseSequence()
			{
				CFID& fid = *(m_signals[0].begin()); 
				return fid.GetPulseSequence();
			}

			//! Return the (non-summed) matrix of group basis vectors.
			const cvm::cmatrix& GetGroupMatrix() const
			{
				assert( m_matGroups.msize() > 0 );
				return m_matGroups;
			}

			//! Return the matrix that maps groups to metabolites.
			const cvm::cmatrix& GetSummationMatrix() const
			{
				assert( m_matSum.msize() > 0 );
				return m_matSum;
			}

			const std::vector<bool>& GetBroadSig() const
			{
				return m_broad_sig;
			}

			void SetBroadSig(std::vector<bool> broad_sig)
			{
				m_broad_sig = broad_sig;
			}

			//! Map from column of group matrix to column of basis matrix.
			integer GetBasisFromGroup(integer nGroupCol) const 
			{
				assert( nGroupCol < (integer)m_indexTracker.size() );
				return m_indexTracker[nGroupCol];
			}

			//! Return the name, worked out from the filename, of the metabolite.
			std::string GetSignalName(std::size_t n) const
			{
				assert( n < m_vecSignalFiles.size() );

				// name with extension
				std::string strFile = m_vecSignalFiles[n].filename().string();

				// without extension 
				return GetFilenameBase(strFile);
			}

			//! Return a vector of strings which correspond to signal names
			std::vector<std::string> GetSignalNamesSTL() const
			{
				return m_vecSignalFilesSTL;
			}

			void SetSignalNamesSTL(std::vector<std::string> sig_files)
			{
				m_vecSignalFilesSTL = sig_files;
			}

			//! Return a vector of strings which correspond to signal names
			std::vector<std::string> GetSignalNames() const
			{
				std::vector<std::string> strFiles;

				//std::cout << m_vecSignalFiles.size() << std::endl;

				// get vector of names 
				for( size_t n = 0 ; n < m_vecSignalFiles.size() ; n++ )
				{	
					std::string signal_name = GetFilenameBase(m_vecSignalFiles[n].filename().string());
					strFiles.push_back(signal_name);
				}

				// without extension 
				return strFiles;
			}

			
			//! Adjust m_matGroups, m_matSum and m_indexTracker such that each metabolite has only one group (the summation of any subgoups defined in the spreadsheet) 
			void CompressGroupMatrix();
			
            void SaveTxt(std::string file_name, CFID fid);
            
            void SaveLCM(std::string file_name, CFID fid);

			const std::string& GetBasisPath() const
			{
				return m_strBasisPath;
			}

			void SetBasisPath(std::string basis_path)
			{
				m_strBasisPath = basis_path;
			}

private:

			void applyRef(treal new_ref);

			//! Populates the basis matrix with the data from the signals produced at simulation.
			bool makeBasisMatrix();

			//! Populates the group matrix with the data from the signals produced at simulation.
			void makeGroupMatrix();



		private:

			//! Location of directory from which this basis is generated.
			std::string m_strBasisPath;

			//! Collection of boost filesystem path objects, one for each file in the basis.
			std::vector<boost::filesystem::path> m_vecSignalFiles; 

			// For serialisation 
			std::vector<std::string> m_vecSignalFilesSTL;

			//! Collection of signals (the results of simulation).
			std::vector<CSignal> m_signals;

			//! A vector describing which signals are broad
			std::vector<bool> m_broad_sig;

			//! Basis vector matrix.
			cvm::cmatrix m_matBasis;

			//! Group matrix.
			cvm::cmatrix m_matGroups;

			//! Summation matrix: basis vector matrix = summation matrix * group matrix.
			cvm::cmatrix m_matSum;

			//! Map from group col to basis col.
			std::vector<integer> m_indexTracker;

	};

}

#endif
