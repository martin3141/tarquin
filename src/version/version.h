#ifndef TARQUIN_VERSION_INCLUDED
#define TARQUIN_VERSION_INCLUDED

#include <string>

namespace tarquin
{
	namespace version
	{
		std::string revision();

		std::string display_name();

		std::string copyright();

		std::string major_version();

		std::string minor_version();

		std::string version_string();

	} // namespace version

} // namespace tarquin


#endif
