#ifndef __CBOSWELL__
#define __CBOSWELL__

#include <sstream>
#include <string>
#include <stdio.h>
#include <vector>


namespace tarquin 
{
	enum LOG_DESTINATION
	{
		LOG_NOWHERE,
		LOG_STDOUT
	};

	enum LOG_LEVEL 
	{ 
		LOG_ERROR, 
		LOG_WARNING, 
		LOG_INFO 
	};

	//! Debugging levels, from least to most verbose.
	enum DEBUG_LEVEL 
	{
		DEBUG_LEVEL_OFF,
		DEBUG_LEVEL_1,
		DEBUG_LEVEL_2,
		DEBUG_LEVEL_3
	};

/*! 
 * Defines the logging object that records useful messages, etc. 
 */
class CBoswell 
{
	public:

		CBoswell(const LOG_DESTINATION destination);
		virtual ~CBoswell();

		void LogMessage(LOG_LEVEL level, const char* fmt, ...);

		void DebugMessage(DEBUG_LEVEL level, const char* fmt, ...);

		//! Get the current time as a string for nice printing.
		std::string GetTime();

		/*! 
		 * Returns the string stream to which messages are written.
		 * \param level controls what kind of message this is.
		 */
		std::ostringstream& Out(LOG_LEVEL level);

		//! Writes the log to all file pointers.
		void Flush();

		//! Writes the log to the specified file pointer (e.g. stdout).
		void Flush(FILE* fp);

		//! Writes the log to the specified string (useful for GUI).
		void Flush(std::string& strLog);

		/*! 
		 * Marks the beginning of a long process that we don't know how long it will take.
		 * This will go to the log and stdout, be default.
		 */
		virtual void BeginTask(std::string strMessage);

		/*! 
		 * Marks intermediate points in a long process. This only goes to stdout.
		 */
		virtual void UpdateTask(std::string strMessage);

		/*! 
		 * Marks the end of a long process. This will go to the log
		 * and stdout.
		 */
		virtual void EndTask(std::string strMessage);

	private:

		// prevent these from being used accidentally
		CBoswell(const CBoswell& in);
		CBoswell& operator=(const CBoswell& rhs);

		void PrintMessage(std::string prefix, const char* fmt, va_list ap);

	protected:

		//! Where the messages go.
		std::ostringstream m_os;

		//! Records the current logging level (for compound expressions).
		LOG_LEVEL m_currentLevel;

		//! The list of things we should flush messages.
		std::vector<FILE*> m_files;

		DEBUG_LEVEL m_debug_level;

};

/*!
 * \brief RAII style log message.
 *
 * Create one of these at the beginning of the task, giving it the first
 * and final message. Call operator() to increment the task. When it goes
 * out of scope the final message will be written.
 */
class Task
{
	public:

		Task(CBoswell& log, std::string task_name) :
			m_log(log),
			m_task_name(task_name)
		{
			m_log.LogMessage(LOG_INFO, "Starting '%s'", task_name.c_str());
		}

		void operator() (std::string msg=".")
		{
			fprintf(stderr, "%s", msg.c_str());
		}

		~Task()
		{
			m_log.LogMessage(LOG_INFO, "Finished '%s'", m_task_name.c_str());
		}

	private:

		CBoswell& m_log;

		std::string m_task_name;

};

}

#endif
