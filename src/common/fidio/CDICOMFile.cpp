#include "CDICOMFile.hpp"
#include "common.hpp"
#include "CFIDReader.hpp"

#include <iostream>
#include <sstream>
//#include "pstdint.h"
#include <stdint.h>

tarquin::CDICOMFile::CDICOMFile()
{
	// short but padded
	m_padded_vrs.push_back("OB");
	m_padded_vrs.push_back("OW");
	m_padded_vrs.push_back("OF");
	m_padded_vrs.push_back("SQ");
	m_padded_vrs.push_back("UN");
	m_padded_vrs.push_back("UT");

	// add all possible value representation types
	m_all_vrs.push_back("AE");
	m_all_vrs.push_back("AS");
	m_all_vrs.push_back("AT");
	m_all_vrs.push_back("CS");
	m_all_vrs.push_back("DA");
	m_all_vrs.push_back("DS");
	m_all_vrs.push_back("DT");
	m_all_vrs.push_back("FL");
	m_all_vrs.push_back("FD");
	m_all_vrs.push_back("IS");
	m_all_vrs.push_back("LO");
	m_all_vrs.push_back("LT");
	m_all_vrs.push_back("OB");
	m_all_vrs.push_back("OF");
	m_all_vrs.push_back("OW");
	m_all_vrs.push_back("PN");
	m_all_vrs.push_back("SH");
	m_all_vrs.push_back("SL");
	m_all_vrs.push_back("SQ");
	m_all_vrs.push_back("ST");
	m_all_vrs.push_back("SS");
	m_all_vrs.push_back("TM");
	m_all_vrs.push_back("UI");
	m_all_vrs.push_back("UL");
	m_all_vrs.push_back("UN");
	m_all_vrs.push_back("US");
	m_all_vrs.push_back("UT");
}

void tarquin::CDICOMFile::Close()
{
    m_file.close();
}

// open the file and skip to where the proper DICOM data begins
// this should work for both part 10 and none part 10 files
void tarquin::CDICOMFile::Open(std::string strFilename)
{
	m_strFilename = strFilename;

	// attempt to open the file
	m_file.open(strFilename.c_str(), std::ios::binary);

	if( m_file.fail() )
		throw CFIDReader::Exception("failed to open file: %s", strFilename.c_str());

	// see if this is a part-10 style file first - skip the first 128 bytes
	m_file.seekg(128, std::ios_base::beg);

	// can we find 'DICM' ?
	char szPart10Tag[5] = {0};
	m_file.read(szPart10Tag, 4);
	szPart10Tag[4] = '\0';

    //std::cout << std::endl;
    //std::cout << szPart10Tag << std::endl;

	// just in case we got passed a really small file
	if( m_file.eof() )
		throw CFIDReader::Exception("reached end of file without finding magic bytes, probably not a DICOM file");

	// it was a a part 10 file, so we return now
	if( std::string(szPart10Tag) == "DICM" )
		return;

	// we couldn't find the part 10 marker, so we should skip back to the begininng of the file
	m_file.seekg(0, std::ios_base::beg);
}

// advance m_file to the (group, element) specified
// returns the length of the field, if the field is not found then the length is -1
long tarquin::CDICOMFile::MoveToTag(std::string strGroup, std::string strElement, bool from_start)
{
	// is the user being a nitwit?
	if( !m_file.is_open() ) 
		throw CFIDReader::Exception("file not open before moving to tag (%s, %s)", strGroup.c_str(), strElement.c_str());

    // go back to the start of the file	
    if ( from_start )
    {
        m_file.seekg(128, std::ios_base::beg);
        char szPart10Tag[5] = {0};
        m_file.read(szPart10Tag, 4);
        szPart10Tag[4] = '\0';
        if( std::string(szPart10Tag) != "DICM" )
            m_file.seekg(0, std::ios_base::beg);
    }

	while( !m_file.eof() ) 
	{
		uint16_t wGroup   = 0;
		uint16_t wElement = 0;
		char szCurrentVR[3] = {0};

		// get the tag and VR at this point in the file
		m_file.read((char*)&wGroup,   sizeof(uint16_t));
		m_file.read((char*)&wElement, sizeof(uint16_t));
		m_file.read(szCurrentVR, 2);


		// get current group and element as strings
		std::string strCurrentGroup   = convertToStringHex(wGroup);
		std::string strCurrentElement = convertToStringHex(wElement);

        //std::cout << "(" << strCurrentGroup << "," << strCurrentElement << ")" << std::endl;

		// this will be set to the length of the field
		uint32_t value_length = 0;

		//
		// an explicit VR
		//
		if( is_explicit_vr(szCurrentVR) ) 
		{
			// 32-bit length with two bytes of padding
			if( is_vr_with_padding(szCurrentVR) )
			{
				// skip the two bytes
				m_file.seekg(2, std::ios_base::cur);

				// read 32-bits length
				m_file.read((char*)&value_length, sizeof(uint32_t));
			}
			// 16-bit length with no padding
			else
			{
				uint16_t length = 0;
				m_file.read((char*)&length, sizeof(uint16_t));
				value_length = length;
			}
		}

		//
		// an implicit VR, the value length is 32 bytes
		// 
		else
		{
            //std::cout << std::endl << "Implicit VR";
			// adjust back to bytes the thing we just read in tentatively as being the VR
			m_file.seekg(-2, std::ios_base::cur);

			// read in length of this field
			m_file.read((char*)&value_length, sizeof(uint32_t));

            // changed MW Apr 2011, changed back again on 20th Feb 2012
            // to fix problem with data from PukkaJ
            //value_length = 0;
		}

		//fprintf(stdout, "\n(G,E) = (%s,%s), VR = %s, length = 0x%x", strCurrentGroup.c_str(), strCurrentElement.c_str(), szCurrentVR, value_length);
		//fprintf(stdout, "\n(G,E) = (%s,%s), VR = %s, length = 0x%x", strCurrentGroup.c_str(), strCurrentElement.c_str(), szCurrentVR, value_length);

		// If this is a sequence element, the 'length' may be undefined, i.e. have this value,

        // If we're on a sequence then decend
		if( std::string(szCurrentVR) == "SQ" ) 
			continue;

		if( value_length == uint32_t(0xFFFFFFFF) ) 
			continue;

		if ( m_file.tellg() == static_cast<std::ifstream::pos_type>(-1) )
		    return -1;
        
        if( strCurrentGroup == "FFFE" && strCurrentElement == "E000" ) 
			continue;
        
		// is this the field we are looking for?
		if( strGroup == strCurrentGroup && strElement == strCurrentElement ) 
			return value_length;
		else  // no, so skip this
			m_file.seekg(value_length, std::ios_base::cur);
	}

	// if we got here then the tag was not found
	return -1;
}

bool tarquin::CDICOMFile::is_explicit_vr(std::string strVR)
{
	return std::find(m_all_vrs.begin(), m_all_vrs.end(), strVR) != m_all_vrs.end();
}

bool tarquin::CDICOMFile::is_vr_with_padding(std::string strVR)
{
	return std::find(m_padded_vrs.begin(), m_padded_vrs.end(), strVR) != m_padded_vrs.end();
}

template<typename in_t>
std::string tarquin::CDICOMFile::convertToStringHex(in_t number)
{
	std::ostringstream strm;
	strm.width(4);
	strm.fill('0');
	strm.flags( std::ios_base::hex | std::ios_base::uppercase );

	strm << number;
	return strm.str();
}

