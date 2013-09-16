#include "CBoswell.hpp"
#include <time.h>
#include <stdarg.h>
#include <boost/date_time.hpp>
#include "version/version.h"

namespace bp = boost::posix_time;

tarquin::CBoswell::CBoswell(const LOG_DESTINATION destination) : 
	m_currentLevel(LOG_INFO),
	m_debug_level(DEBUG_LEVEL_OFF)
{
	if( LOG_STDOUT == destination )
		m_files.push_back(stdout);

	LogMessage(LOG_INFO, "%s %s Started", tarquin::version::display_name().c_str(), tarquin::version::version_string().c_str());
}

tarquin::CBoswell::~CBoswell()
{
	LogMessage(LOG_INFO, "%s Finished\n", tarquin::version::display_name().c_str());
	Flush();
}

/*!
 * Write a message to the application log. This is a message that will occur in the normal
 * execution, i.e. it is not a debugging message, so only write things using this function
 * that most users will care about.
 */
void tarquin::CBoswell::LogMessage(LOG_LEVEL level, const char* fmt, ...)
{
	std::string out_msg;

	switch( level )
	{
		case LOG_ERROR:
			out_msg += "ERROR  :";
			break;

		case LOG_WARNING:
			out_msg += "WARNING:";
			break;

		case LOG_INFO:
			out_msg += "INFO   :";
			break;
	}


	// print the message to the buffer
	va_list args;
	va_start(args, fmt);
	PrintMessage(out_msg, fmt, args);
	va_end(args);
}

void tarquin::CBoswell::DebugMessage(DEBUG_LEVEL level, const char* fmt, ...)
{
	// do nothing if we are not outputting at the specified debug level
	if( level < m_debug_level )
		return;
		
	
	// todo, test if we are below the debugging level here
	std::string prefix = "DEBUG  :";

	// print the message to the buffer
	va_list args;
	va_start(args, fmt);
	PrintMessage(prefix, fmt, args);
	va_end(args);
}

void tarquin::CBoswell::PrintMessage(std::string prefix, const char* fmt, va_list args)
{
	// what is the time now?
	bp::ptime time_now = bp::second_clock::local_time();

	// the string that we will be displaying
	std::string out_msg = bp::to_simple_string(time_now);

	char buff[4096];
	vsnprintf(buff, sizeof(buff)-1, fmt, args); 

	out_msg += " ";
	out_msg += prefix;
	out_msg += " ";
	out_msg += buff;

	// print to all the files
	for( size_t i = 0; i < m_files.size(); ++i )
	{
		fprintf(m_files[i], "\n%s", out_msg.c_str());
		fflush(m_files[i]);
	}
}


void tarquin::CBoswell::BeginTask(std::string strMessage)
{
	LogMessage(LOG_INFO, strMessage.c_str());
}

void tarquin::CBoswell::UpdateTask(std::string strMessage)
{
	// no output, so don't generate the progress dots
	if( !m_files.size() )
		return;
	
	fprintf(stderr, "%s", strMessage.c_str());
	fflush(stderr);
}

void tarquin::CBoswell::EndTask(std::string strMessage)
{
	LogMessage(LOG_INFO, strMessage.c_str());
}

std::string tarquin::CBoswell::GetTime()
{
	time_t raw;
	time(&raw);

	std::string str(ctime(&raw));
	return str;
}

std::ostringstream& tarquin::CBoswell::Out(LOG_LEVEL level)
{
	// grab the current level and in case of compound expressions
	m_currentLevel = level;

	// stuff that goes to the log have useful prefixes
	if( LOG_ERROR == level ) {

		m_os << "\n" << "ERROR:   ";
	}
	else if( LOG_WARNING == level ) {

		m_os << "\n" << "WARNING: ";
	}
	else { // LOG_INFO 

		// no prefixes, user text should come out natively
	}

	return m_os;
}

void tarquin::CBoswell::Flush(FILE* fp)
{
	fprintf(fp, "%s", m_os.str().c_str());
	fflush(fp);

	// the log is now empty
	m_os.str("");
}

void tarquin::CBoswell::Flush(std::string& strLog)
{
	strLog += m_os.str();
	m_os.str("");
}

void tarquin::CBoswell::Flush()
{
	// write it to all the streams
	for( size_t i = 0; i < m_files.size(); ++i )
	{
		fprintf(m_files[i], "%s", m_os.str().c_str());
	}

	m_os.str("");
}

