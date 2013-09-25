#include "version.h"
#include "version_constants.h"

// 
// Reposistory revision
//
std::string tarquin::version::revision()
{
	return SVN_REVISION;
}

//
// TARQUIN
//
std::string tarquin::version::display_name()
{
	return "TARQUIN";
}

std::string tarquin::version::major_version()
{
	return TARQUIN_MAJOR_VERSION;
}

std::string tarquin::version::minor_version()
{
	return TARQUIN_MINOR_VERSION;
}

std::string tarquin::version::version_string()
{
	return major_version() + "." + minor_version(); //+ "." + revision();
}

std::string tarquin::version::copyright()
{
	return "(c) Copyright 2006-2013 Greg Reynolds and Martin Wilson.";
}
