project(version)

# make a file that contains the version string
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/version_constants.h.in ${CMAKE_CURRENT_SOURCE_DIR}/version_constants.h)

# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# All files in this directory to be compiled into the library
file(GLOB version_SOURCES *.cpp)

# Create an executable file called helloworld from sources:
add_library(version ${version_SOURCES})
