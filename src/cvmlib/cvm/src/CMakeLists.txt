project(cvm)

# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# All files in this directory to be compiled into the library
file(GLOB cvm_SOURCES *.cpp)

# This is the library
add_library(cvm ${cvm_SOURCES})

# The CVM library needs these dependencies.
target_link_libraries(cvm 
	ftn
	${BLAS_LIB}
	${LAPACK_LIB}
	${F2C_LIB}
)

