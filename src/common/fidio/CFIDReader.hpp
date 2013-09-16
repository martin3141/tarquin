#ifndef __CFIDREADER__
#define __CFIDREADER__

#include <boost/tokenizer.hpp>
#include "common.hpp"
#include "exception.hpp"
#include "Options.hpp"

namespace tarquin 
{

class CFID;
class CBoswell;

/*!
 * Abstract base class of all FID parsers. Parsers specific to each parser derive
 * from this object.
 */
class CFIDReader 
{
	public:

		/*!
		 * \brief Exceptions thrown during reading are of this type.
		 */
		class Exception : public tarquin::Exception
		{
			public:

				Exception(const char* szfmt, ...) 
				{
					va_list lst;
					va_start(lst, szfmt);
					Message(szfmt, lst);
					va_end(lst);
				}
		};

protected:

    //! The pair of token and line number. 
    typedef std::pair<std::string, std::size_t> Token;

    //! A collection of tokens.
    typedef std::vector<Token> TokenList;

protected:

    //! The FID we are going to read into.
    CFID& m_fid;

    //! Where messages go.
    CBoswell& m_log;

    //! The list of tokens in the file.
    TokenList m_tokens;

protected:

    CFIDReader(CFID& fid, CBoswell& log);

    //! Populates the token list in the reader object, given the separator list.
    virtual void Lex(std::string strFilename, const boost::escaped_list_separator<char>& sep);

    virtual ~CFIDReader() { };

public:

    //virtual void Load(std::string strFilename, CBoswell& log) = 0;
    virtual void Load(std::string strFilename, const Options& opts, CBoswell& log) = 0;
};

}
#endif
