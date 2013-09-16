#include "CFIDReaderPhilipsDCM.hpp"
#include "CFID.hpp"
#include "CBoswell.hpp"
#include "CDICOMFile.hpp"

#include <iostream>

tarquin::CFIDReaderPhilipsDCM::CFIDReaderPhilipsDCM(CFID& fid, CBoswell& log) : 
		CFIDReader(fid, log)
{
}

void tarquin::CFIDReaderPhilipsDCM::Load(std::string strFilename, const Options& opts, CBoswell& log)
{
	size_t file_size = GetFileSize(strFilename);

	if( !file_size )
		throw Exception("file size is zero");

	// attempt to open the file in DICOM format
	CDICOMFile file;

	file.Open(strFilename);
    
    // field strength tag
    long ft_size = file.MoveToTag("2001", "1083"); 
	std::string ft(ft_size/8, '\0'); 
	file.GetFileStream().read(&ft[0], ft_size);
    std::string ft_str = ft;
    ft_str = ft_str.substr(0,16);
    float ft_float;
    std::istringstream iss(ft_str, std::istringstream::in);
    iss >> ft_float;
    m_fid.SetTransmitterFrequency(ft_float*1e6);

    // find sampling frequency tag
    file.MoveToTag("2005", "1030"); 
	float fs1 = 0;
	float fs2 = 0;
	file.GetFileStream().read((char*)&fs1, 4);
	file.GetFileStream().read((char*)&fs2, 4);
    m_fid.SetSamplingFrequency(fs1);

	// move to point where FID starts
	long nBytesInFID = file.MoveToTag("2005", "1270"); 

	if( -1 == nBytesInFID ) 
	{
		throw Exception("file did not contain the magic tag (2005, 1270)");
	}

	ReadFIDData(file.GetFileStream(), nBytesInFID);
    
    // find echo time tag
    file.MoveToTag("2005", "1310"); 
	float echo;
	file.GetFileStream().read((char*)&echo, 4);
    m_fid.SetEchoTime(echo/1000);
}

void tarquin::CFIDReaderPhilipsDCM::ReadFIDData(std::ifstream& file, std::size_t nLength)
{
	// floats are in IEEE 32 bit, little endian
	// (Siemens scanners use glorified PCs)
	float real_part_ieee;
	float imag_part_ieee;

	long nSamples = nLength / 8;
	int n = 0;

	cvm::cvector FID(nSamples);

	for( std::size_t nBytesRead = 0; nBytesRead < nLength; nBytesRead+= 8 ) 
	{
		file.read((char*)&real_part_ieee, 4);
		file.read((char*)&imag_part_ieee, 4);

		if( true == file.eof() ) 
			throw Exception("ran out file before reading in target number of bytes");

		// only for big endian boxes (Sun?)
#ifdef __BIG_ENDIAN__
		BYTESWAP(real_part_ieee);
		BYTESWAP(imag_part_ieee);
#endif

		// store this sample
		FID[n+1] =  tcomplex(real_part_ieee, -imag_part_ieee);
		n++;
	}

	m_fid.AppendFromVector(FID);
}

template <class T>
bool tarquin::CFIDReaderPhilipsDCM::from_string(T& t, const std::string& s, std::ios_base& (*f)(std::ios_base&))
{
	std::istringstream iss(s);
	return !(iss >> f >> t).fail();
}

