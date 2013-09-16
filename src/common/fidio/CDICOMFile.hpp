#ifndef __CDICOMFILE__
#define __CDICOMFILE__

#include <string>
#include <fstream>
#include <vector>

namespace tarquin 
{

/*!
 * Class for working with files in DICOM format.
 */
class CDICOMFile 
{

	public:

		//! Constructor.
		CDICOMFile();

		//! Destructor.
		~CDICOMFile()
		{

		}

		//! Open a file.
		void Open(std::string strFilename);
		
		//! Close a file.
        void Close();

		//! Move the file to a particular (group, element) and return the length of the field, -1 on fail. 
		long MoveToTag(std::string strGroup, std::string strElement, bool from_start=true );
		
		//! Get the ifstream inside this file, the one that has been opened and moved.
		inline std::ifstream& GetFileStream()
		{
			return m_file;
		}

		//! True if the VR passed in is explicit.
		bool is_explicit_vr(std::string strVR);

		//! True if the length field for the VR passed in is 2 bytes (rather than 4).
		bool is_short_vr(std::string strVR);

		//! True if the VR passed is short but has two bytes of padding.
		bool is_vr_with_padding(std::string strVR);

		//! Convert a number type to a hex string of format %04X
		template<typename in_t>
			std::string convertToStringHex(in_t number);

	protected:

		//! Name of the file we have opened.
		std::string m_strFilename;

		//! File we have opened.
		std::ifstream m_file;

		//! VRs of length 4 with 2 bytes padding.
		std::vector<std::string> m_padded_vrs;

		//! Collection of all possibile VRs.
		std::vector<std::string> m_all_vrs;

};

}

#endif
