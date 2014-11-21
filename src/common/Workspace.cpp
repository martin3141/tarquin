#include "Workspace.hpp"
#include "levmar.h" 


bool tarquin::Workspace::CheckParameters()
{
	bool bCondition = 
		m_fidRaw.IsKnownEchoTime()             && 
		m_fidRaw.IsKnownSamplingFrequency()    &&
		m_fidRaw.IsKnownTransmitterFrequency() && 
		m_fidRaw.GetNumberOfPoints() > 0;
		//m_options.GetBasisPath() != "";

	if( bCondition )
		return true;
	else
		return false;
}

namespace tarquin
{

std::ostream& operator<< (std::ostream& os, const tarquin::Workspace& rhs)
{
	const Options& options = rhs.GetOptions();
	const CFID& fid         = rhs.GetFID();
	const CFID& fidwater    = rhs.GetFIDWater();

	os << "\nInput file:                 " << options.GetFilename();

	if( options.GetFilenameWater() != "" )
		os << "\nWater reference file:       " << options.GetFilenameWater();
	else
		os << "\nWater reference file:       " << "none";


	/*os << "\nInput file format:          ";
	if( tarquin::SIEMENS == options.GetFormat() )
		os << "siemens";
	else if( tarquin::PHILIPS == options.GetFormat() )
		os << "philips";
	else if( tarquin::GE == options.GetFormat() )
		os << "ge";
	else if( tarquin::DANGER == options.GetFormat() )
		os << "dangerplot";
	else if( tarquin::VARIAN == options.GetFormat() )
		os << "varian";*/

	//os << "\nInfinity Norm of WS vector: " << fid.GetNormValue();
	//if( options.GetFilenameWater() != "" )
	//	os << "\nInfinity Norm of W vector:  " << fidwater.GetNormValue();

	os <<     "\nData points:                " << fid.GetNumberOfPoints();
	os <<     "\nRows:                       " << fid.GetRows();
	os <<     "\nColumns:                    " << fid.GetCols();
	os <<     "\nSlices:                     " << fid.GetSlices();
	if( fid.IsKnownEchoTime() )
		os << "\nEcho time:                  " << fid.GetEchoTime() << " s";
	else 
		os << "\nEcho time:                  " << "unknown (needed!)"; 

	if( fid.IsKnownSamplingFrequency() )
		os << "\nSampling frequency:         " << fid.GetSamplingFrequency() << " Hz";
	else
		os << "\nSampling frequency:         " << "unknown (needed!)";

	if( fid.IsKnownTransmitterFrequency() )
		os << "\nTransmitter frequency:      " << fid.GetTransmitterFrequency() << " Hz";
	else
		os << "\nTransmitter frequency:      " << "unknown (needed!)";

	if( 0 != options.GetRangeStart() )
		os << "\nStarting sample:            " << options.GetRangeStart();
	else 
		os << "\nStarting sample:            " << "unknown (implies bad FID data!)";

	if( 0 != options.GetRangeEnd() )
		os << "\nEnding sample:              " << options.GetRangeEnd();
	else
		os << "\nEnding sample:              " << "unknown (implies bad FID data!)";

	if( 0 != options.GetWaterWindow() ) {
		os << "\nRemoving signals in:        " << "[-" << options.GetWaterWindow()
			<< ",+" << options.GetWaterWindow() << "] Hz";
	}

	if( 0 != options.GetConvWindowWidth() ) {
		os << "\nConvolution window width:   " << options.GetConvWindowWidth()
			<< " points";
	}

	if( true == options.GetAutoPhase() ) 
		os << "\nAutophasing yes/no?:        " << "yes";
	else
		os << "\nAutophasing yes/no?:        " << "no";

	if( true == options.GetWaterEddy() ) 
		os << "\nWater eddy correction?:     " << "yes";
	else
		os << "\nWater eddy correction?:     " << "no";

    /*
	if( true == fid.IsKnownPhi0() )
		os << "\nPhi 0:                      " << fid.GetPhi0();
	else
		os << "\nPhi 0:                      " << "unknown (assuming zero)";

	if( true == fid.IsKnownPhi1() )
		os << "\nPhi 1:                      " << fid.GetPhi1();
	else
		os << "\nPhi 1:                      " << "unknown (assuming zero)";
        */ // TODO

	if( true == options.GetAutoReference() )
		os << "\nAutoreferencing yes/no?:    " << "yes";
	else
		os << "\nAutoreferencing yes/no?:    " << "no";

	pair_vec ref = fid.GetPPMRef();
	if( true == ref[0].second )
		os << "\nReference:                  " << ref[0].first << " ppm";
	else
		os << "\nReference:                  " << "unknown (needed!)";

	/*if( options.GetUsePrecompiled() )
		os << "\nPrecompiled basis file:     " << options.GetBasisPath();
	else
		os << "\nBasis directory:            " << options.GetBasisPath();
        */
	//os << "\nOutput file:                " << options.GetOutputXMLPath();

	os << "\n" << std::endl;

	return os;
}

}

