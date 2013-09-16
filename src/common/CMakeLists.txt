project(tarquincore)

# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# All files in this directory to be compiled into the library
file(GLOB core_SOURCES *.cpp)
file(GLOB fidio_SOURCES fidio/*.cpp)

# Needs some headers
include_directories(${CVM_INCLUDE})
include_directories(${LEVMAR_INCLUDE})
include_directories(${FAST_SIM_INCLUDE})
include_directories(${EIGEN3_INCLUDE})

# This is the library
add_library(tarquincore 
	${core_SOURCES} 
	${fidio_SOURCES}
)

target_link_libraries(tarquincore 
    ${PLATFORM_THREAD_LIBRARY}
	fast_sim
	${CVM_LIB}
	levmar
	version
	${FFTW_LIB}
)

