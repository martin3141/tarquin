project(tarquin)


cmake_minimum_required(VERSION 2.8)

#
# Setup platform specific paths for Linux
#
if( ${CMAKE_SYSTEM} MATCHES "Linux" )

	add_definitions(-DPLATFORM_LINUX)

	#add_definitions(-ffloat-store)
	add_definitions(-msse2)
	add_definitions(-mfpmath=sse)
	add_definitions(-Wall)

	set( BLAS_LIB   -lblas )
	#set( BLAS_LIB   /usr/lib/libblas.a )

	set( LAPACK_LIB -llapack )
	#set( LAPACK_LIB /usr/lib/liblapack.a )

    set( PLATFORM_THREAD_LIBRARY -lpthread )
    
	set( FFTW_INCLUDE /usr/include )

	set( FFTW_LIB     -lfftw3 )
	#set( FFTW_LIB     /usr/lib/i386-linux-gnu/libfftw3.a )

	set( QWT_INCLUDE /usr/include/qwt-qt4 )
	set( QWT_LIB     -lqwt-qt4 )

	#set( QWT_INCLUDE /usr/local/qwt-5.2.3/include )
	#set( QWT_LIB     /usr/local/qwt-5.2.3/lib/libqwt.a )

	set( QWT_LIB_RELEASE ${QWT_LIB} )
	set( QWT_LIB_DEBUG   ${QWT_LIB} )

	# this find package works on Linux
	find_package(Protobuf 2.3.0 REQUIRED)
    
    # dynamic lib
	set( PROTOBUF_LIB_DEBUG   ${PROTOBUF_LIBRARY} )
	set( PROTOBUF_LIB_RELEASE ${PROTOBUF_LIBRARY} )
    
    # static lib
    #set( PROTOBUF_LIB_DEBUG   /usr/lib/libprotobuf.a )
	#set( PROTOBUF_LIB_RELEASE /usr/lib/libprotobuf.a )


endif()

#
# Setup platform specific paths for OSX
#
if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )

	add_definitions(-DPLATFORM_LINUX)
	
    add_definitions(-DTIME_UTC_CHANGE)

	#add_definitions(-isysroot /Developer/SDKs/MacOSX10.5.sdk)
	#add_definitions(-mmacosx-version-min=10.5)
	#add_definitions(-m32)
	#add_definitions(-arch i386)
	add_definitions(-arch x86_64)
	#add_definitions(-no_compact_linkedit)

	#add_definitions(-ffloat-store)
	add_definitions(-msse2)
	add_definitions(-mfpmath=sse)
	add_definitions(-Wall)

	set( BLAS_LIB   -lblas )
	set( LAPACK_LIB -llapack )

	set( FFTW_INCLUDE /usr/include )
	set( FFTW_LIB     -lfftw3 )

	set( QWT_INCLUDE /usr/local/qwt-5.2.1/include )
	set( QWT_LIB     /usr/local/qwt-5.2.1/lib/libqwt.a )
	#set( QWT_LIB     /usr/local/qwt-5.2.1/lib/libqwt.dylib )
	#set( QWT_LIB     -lqwt )
	set( QWT_LIB_RELEASE ${QWT_LIB} )
	set( QWT_LIB_DEBUG   ${QWT_LIB} )

	# this find package works on OSX
	find_package(Protobuf 2.3.0 REQUIRED)

	#set( PROTOBUF_LIB_DEBUG   ${PROTOBUF_LIBRARY} )
	#set( PROTOBUF_LIB_RELEASE ${PROTOBUF_LIBRARY} )
    
	set( PROTOBUF_LIB_DEBUG   /usr/local/lib/libprotobuf.a )
	set( PROTOBUF_LIB_RELEASE /usr/local/lib/libprotobuf.a )

endif()

#
# Setup platform specific paths for Windows
#
if( WIN32 )

	add_definitions(-DPLATFORM_WINDOWS)
    add_definitions(-DTIME_UTC_CHANGE)
	add_definitions(-DNOMINMAX)
	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
	add_definitions(-DEIGEN_SSE2_ON_MSVC_2008_OR_LATER)

	set( BOOST_ROOT                "c:/opt/boost-tarquin")

	set( TARQUIN_DEPS_DIR          "d:/libs")

	set( LAPACK_ROOT               "${TARQUIN_DEPS_DIR}/clapack-3.2.1-CMAKE/build" )
	set( BLAS_LIB                  "${LAPACK_ROOT}/BLAS/SRC/Release/blas.lib" )
	set( LAPACK_LIB                "${LAPACK_ROOT}/SRC/Release/lapack.lib" )
	set( F2C_LIB                   "${LAPACK_ROOT}/F2CLIBS/libf2c/Release/libf2c.lib" )

	set( CVM_BASE                  "${TARQUIN_DEPS_DIR}/cvmlib" )
	set( CVM_INCLUDE               "${CVM_BASE}" )
	set( CVM_LIB                   "${CVM_BASE}/cvm_ia32.lib" )

	set( FFTW_ROOT                 "${TARQUIN_DEPS_DIR}/fftw-3.2.2/" ) 
	set( FFTW_INCLUDE              "${FFTW_ROOT}/api" )
	set( FFTW_LIBRARY_DIR          "${FFTW_ROOT}/fftw-3.2.2-libs/Release/" )
	set( FFTW_LIB                  "${FFTW_LIBRARY_DIR}/libfftw-3.2.2.lib" "${FFTW_LIBRARY_DIR}/libfftwf-3.2.2.lib" )

	set( PROTOBUF_BASE_DIR         "${TARQUIN_DEPS_DIR}/protobuf-2.3.0/")
	set( PROTOBUF_INCLUDE_DIRS     "${PROTOBUF_BASE_DIR}/src/" )
	set( PROTOBUF_LIB_RELEASE  "${PROTOBUF_BASE_DIR}/vsprojects/Release/libprotobuf.lib" )
	set( PROTOBUF_LIB_DEBUG    "${PROTOBUF_BASE_DIR}/vsprojects/Debug/libprotobuf.lib" )

	set( QWT_ROOT                  "${TARQUIN_DEPS_DIR}/qwt-5.2")
	set( QWT_INCLUDE               "${QWT_ROOT}/src" )
	set( QWT_LIB_RELEASE           "${QWT_ROOT}/lib/qwt.lib" )
	set( QWT_LIB_DEBUG             "${QWT_ROOT}/lib/qwtd.lib" )

endif()

#
# We are using Boost
#
set(Boost_USE_STATIC_LIBS   ON)
set(Boost_USE_MULTITHREADED ON)

# Changed by Nitram to allow a build on Ubuntu 8.04
find_package( 
	Boost 1.46.0 REQUIRED COMPONENTS 
	unit_test_framework 
	date_time 
	filesystem 
	system 
	thread
	)

add_definitions(-DBOOST_FILESYSTEM_VERSION=3)

# we will want boost for all projects
include_directories(${Boost_INCLUDE_DIRS})
add_definitions(-DBOOST_ALL_NO_LIB)

#
# We are using Google's protobuf
#
include_directories(${PROTOBUF_INCLUDE_DIRS})

#
# We are still using FFTW 'cos Greggy's FFT still isn't fast enough
#
include_directories(${FFTW_INCLUDE})

# include paths relative to top of the tree
include_directories(${CMAKE_SOURCE_DIR})

# we are going to link CVMLIB statically
add_definitions(-DCVM_STATIC)

# Some include variables
set(TARQUIN_INCLUDE  ${CMAKE_SOURCE_DIR}/common)
set(LEVMAR_INCLUDE   ${CMAKE_SOURCE_DIR}/levmar/levmar-2.5)
set(FAST_SIM_INCLUDE ${CMAKE_SOURCE_DIR}/fast_sim/fast_sim)

#
# For civilised platforms, turn on warnings
#
if( CMAKE_COMPILER_IS_GNUCXX )
add_definitions(-Wall -Wextra)
endif()

# link location for boost
link_directories(${Boost_LIBRARY_DIRS})

###################################
#
# Setup automatic version numbering
#
#find_package(Subversion REQUIRED)

# grab the current repository version
#Subversion_WC_INFO(${CMAKE_CURRENT_SOURCE_DIR} ER)
#set(SUBVERSION_REVISION ${ER_WC_REVISION})

set(GLOBAL_REV notset)

# The one and only place version numbers are set
set(TARQUIN_MAJOR_VERSION "4")
set(TARQUIN_MINOR_VERSION "3.6")

#
# end of product version numbers
##################################

# testing
enable_testing()

# This is where stuff gets put when it is built 
set(REDIST_TREE_ROOT ${CMAKE_BINARY_DIR}/redist)
file(MAKE_DIRECTORY ${REDIST_TREE_ROOT})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${REDIST_TREE_ROOT})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${REDIST_TREE_ROOT})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${REDIST_TREE_ROOT})

#
# External dependencies that we ship in the tree
#
add_subdirectory(version)

# on Windows we need the prebuilt version of cvmlib, because we
# can't afford the Fortran compiler
if( NOT WIN32 )
	add_subdirectory(cvmlib/cvm)
	set( CVM_INCLUDE ${CMAKE_SOURCE_DIR}/cvmlib/cvm/src)
	set( CVM_LIB     cvm )
endif()

set( EIGEN3_INCLUDE ${CMAKE_SOURCE_DIR}/fast_sim/eigen3_1_1 )

add_subdirectory(levmar/levmar-2.5)

#
# Parts of TARQUIN
#
add_subdirectory(fast_sim)
add_subdirectory(common)
add_subdirectory(cmdline)
add_subdirectory(qtgui)
add_subdirectory(unit_tests)
add_subdirectory(generate_spectra)
