#include "CFIDReader.hpp"
#include "CBoswell.hpp"
#include "CFID.hpp"
#include <boost/tokenizer.hpp>
#include <fstream>
#include <string>

tarquin::CFIDReader::CFIDReader(CFID& fid, CBoswell& log) : 
	m_fid(fid), 
	m_log(log)
{
}

void tarquin::CFIDReader::Lex(std::string strFilename, const boost::escaped_list_separator<char>& sep)
{
    try
    {
        // attempt to open file
        std::ifstream fin(strFilename.c_str());

        if( false == fin.is_open() )
        {
            Exception e("failed to open file: %s", strFilename.c_str()); 
            throw e;
        }

        size_t nLine = 0;

        // read line by line 
        for( std::string strLine; getline(fin, strLine); ) 
        {
            nLine++;

            // skip over empty lines
            if( 0 == strLine.length() )
                continue;

            // skip over comments
            if( '#' == strLine.at(0) )
                continue;

            // tokenize this line
            boost::tokenizer< boost::escaped_list_separator<char> > tokens_line(strLine, sep);

            // variables used when loading 
            boost::tokenizer< boost::escaped_list_separator<char> >::iterator it_token = tokens_line.begin();

            // check for bad file format, i.e. failure to tokenize
            if( it_token == tokens_line.end() ) 
            {
                Exception e("failed to tokenize line %d of %s", nLine, strFilename.c_str()); 
                throw e;
            }


            // add line tokens to top level collection
            for( ; it_token != tokens_line.end(); it_token++ ) 
                if( it_token->length() > 0 ) 
                    m_tokens.push_back(make_pair(*it_token, nLine));
        }
    }
    catch( const std::exception& e )
    {
        throw Exception("Error parsing file, check file format : %s", e.what());
    }
}

