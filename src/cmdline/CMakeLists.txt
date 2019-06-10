project(tarquincmd)

# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# All files in this directory to be compiled into the binary
file(GLOB tarquincmd_SOURCES *.cpp)

# Needs TARQUIN headers
include_directories(${TARQUIN_INCLUDE})
include_directories(${CVM_INCLUDE})

# Target for the console version
add_executable(tarquin ${tarquincmd_SOURCES})

# This library uses cvmlib
target_link_libraries(tarquin 
	tarquincore
	${Boost_FILESYSTEM_LIBRARY} 
	${Boost_DATE_TIME_LIBRARY} 
	${Boost_THREAD_LIBRARY} 
	${Boost_CHRONO_LIBRARY} 
	${Boost_SYSTEM_LIBRARY} 
)
