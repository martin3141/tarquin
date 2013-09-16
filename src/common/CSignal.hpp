#ifndef __CSIGNAL__
#define __CSIGNAL__

#include <vector>
#include "CFID.hpp"

//#include <boost/archive/text_oarchive.hpp>
//#include <boost/archive/text_iarchive.hpp>
//#include <boost/archive/xml_iarchive.hpp>
//#include <boost/archive/xml_oarchive.hpp>
//#include <boost/archive/binary_oarchive.hpp>
//#include <boost/archive/binary_iarchive.hpp>
//#include <boost/serialization/nvp.hpp>
//#include <boost/serialization/vector.hpp>

namespace tarquin {

	/*! 
	 * The class that represents a collection of FIDs. For example, this would represent a metabolite
	 * used in the basis set; a metabolite that consists of a collection of groups. Each group would be a CFID.
	 */

	class CSignal {

		public:

			typedef std::vector<CFID>::iterator fid_iterator;
			typedef std::vector<CFID>::const_iterator const_fid_iterator;

			CSignal() {};

			/*friend class boost::serialization::access;
			template<class Archive>
			void serialize(Archive & ar, const unsigned int)
			{
				ar & BOOST_SERIALIZATION_NVP(m_fids);
			}*/
            
            inline void SetFids(std::vector<CFID> fids)
            {
                m_fids = fids;
            }

			//! The collection of FIDs, that when summed, form this signal.
			std::vector<CFID> m_fids;

			inline fid_iterator begin()
			{
				return m_fids.begin();
			}

			inline fid_iterator end()
			{
				return m_fids.end();
			}

			inline std::size_t size()
			{
				return m_fids.size();
			}
	};


}

#endif
