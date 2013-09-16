project(tarquincmd)

# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# All files in this directory to be compiled into the binary
file(GLOB tarquincmd_SOURCES *.cpp)

# Needs TARQUIN headers
include_directories(${TARQUIN_INCLUDE})
include_directories(${CVM_INCLUDE})

# Target for the console version
add_executable(simulate ${tarquincmd_SOURCES})

# This library uses cvmlib
target_link_libraries(simulate 
	tarquincore
	debug ${PROTOBUF_LIB_DEBUG} optimized ${PROTOBUF_LIB_RELEASE}
	${Boost_FILESYSTEM_LIBRARY} 
	${Boost_DATE_TIME_LIBRARY} 
	${Boost_THREAD_LIBRARY} 
	${Boost_SYSTEM_LIBRARY} 
)
