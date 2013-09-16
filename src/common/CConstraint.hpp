#ifndef __CConstraint__
#define __CConstraint__

#include "common.hpp"

//#include <boost/archive/xml_iarchive.hpp>
//#include <boost/archive/xml_oarchive.hpp>
//#include <boost/archive/text_oarchive.hpp>
//#include <boost/archive/text_iarchive.hpp>
//#include <boost/serialization/nvp.hpp>
//#include <boost/serialization/vector.hpp>

namespace tarquin {

    /*!
     * Class describing the maximum and minimum parameter adjustments of a single basis vector.
     */
    class CConstraint {

	friend class Options;

	public:

	    CConstraint()
	    {
		m_minAlpha = 0;
		m_maxAlpha = 6;
		m_typAlpha = 0;

		m_minBeta = 0;
		m_maxBeta = 100;
		m_typBeta = 0;

		m_minShiftHz = 0;
		m_maxShiftHz = 0;
		m_typShiftHz = 0;
	    }
        
        inline void SetMinAlpha(treal min_alpha)
        {
            m_minAlpha = min_alpha;
        }
        inline void SetMaxAlpha(treal max_alpha)
        {
            m_maxAlpha = max_alpha;
        }
        inline void SetTypAlpha(treal typ_alpha)
        {
            m_typAlpha = typ_alpha;
        }
        inline void SetMinBeta(treal min_beta)
        {
            m_minBeta = min_beta;
        }
        inline void SetMaxBeta(treal max_beta)
        {
            m_maxBeta = max_beta;
        }
        inline void SetTypBeta(treal typ_beta)
        {
            m_typBeta = typ_beta;
        }
        inline void SetMinShiftHz(treal min_shift)
        {
            m_minShiftHz = min_shift;
        }
        inline void SetMaxShiftHz(treal max_shift)
        {
            m_maxShiftHz = max_shift;
        }
        inline void SetTypShiftHz(treal typ_shift)
        {
            m_typShiftHz = typ_shift;
        }

        
        inline treal GetMinAlpha()
        {
            return m_minAlpha;
        }
        inline treal GetMaxAlpha()
        {
            return m_maxAlpha;
        }
        inline treal GetTypAlpha()
        {
            return m_typAlpha;
        }
        inline treal GetMinBeta()
        {
            return m_minBeta;
        }
        inline treal GetMaxBeta()
        {
            return m_maxBeta;
        }
        inline treal GetTypBeta()
        {
            return m_typBeta;
        }
        inline treal GetMinShiftHz()
        {
            return m_minShiftHz;
        }
        inline treal GetMaxShiftHz()
        {
            return m_maxShiftHz;
        }
        inline treal GetTypShiftHz()
        {
            return m_typShiftHz;
        }

	private:

		/*friend class boost::serialization::access;
		template<class Archive>
			void serialize(Archive & ar, const unsigned int)
			{
				ar & BOOST_SERIALIZATION_NVP(m_minAlpha);
				ar & BOOST_SERIALIZATION_NVP(m_maxAlpha);
				ar & BOOST_SERIALIZATION_NVP(m_typAlpha);
				ar & BOOST_SERIALIZATION_NVP(m_minBeta);
				ar & BOOST_SERIALIZATION_NVP(m_maxBeta);
				ar & BOOST_SERIALIZATION_NVP(m_typBeta);
	    		ar & BOOST_SERIALIZATION_NVP(m_minShiftHz);
	    		ar & BOOST_SERIALIZATION_NVP(m_maxShiftHz);
	    		ar & BOOST_SERIALIZATION_NVP(m_typShiftHz);
			}
            */

		//! Minimum Alphaian damping adjustment.
		treal m_minAlpha;

		//! Maximum Alphaian damping adjustment.
	    treal m_maxAlpha;

	    //! Starting value for Alphaian adjustment.
	    treal m_typAlpha;

	    //! Minimum Betaian damping adjustment.
	    treal m_minBeta;

	    //! Maximum Betaian damping adjustment.
	    treal m_maxBeta;

	    //! Starting value for Betaian adjustment.
	    treal m_typBeta;

	    //! Minimum frequency shift in Hz.
	    treal m_minShiftHz;

	    //! Maximum frequency shift in Hz.
	    treal m_maxShiftHz;

	    //! Starting value of frequency shift.
	    treal m_typShiftHz;
    };
}

#endif
