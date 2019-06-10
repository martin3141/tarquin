project(tarquingui)

# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Qt5 package is mandatory
find_package(Qt5 COMPONENTS Core Widgets Concurrent PrintSupport REQUIRED)

# include QT stuff
include_directories(${Qt5Widgets_INCLUDE_DIRS})
include_directories(${Qt5Concurrent_INCLUDE_DIRS})
include_directories(${Qt5PrintSupport_INCLUDE_DIRS})

# include the output directory, because we need the header files generated from the UI files
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# we are using QWT
include_directories(${QWT_INCLUDE})

include_directories(${TARQUIN_INCLUDE})

# I don't know how to make CMake propagate includes when it uses libraries
include_directories(${CVM_INCLUDE})
include_directories(${LEVMAR_INCLUDE})
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE})

# Files for the GUI proper
file(GLOB tarquingui_SRCS *.cpp)

# Files for the plotting widget
file(GLOB plot_SRCS plot/*.cpp)

# Headers that need to go through the QT MOC
set( tarquingui_MOC_HEADERS 
	session.h
	mainwindow.h 
	aboutdlg.h
	inputdlg.h
	resultsdlg.h
	tarquinplotwidget.h
	quickfitdlg.h
	formatdlg.h
	voxel_select.h
)

# Resource files
set( tarquingui_QRC resources/tarquingui.qrc )

# UI files that need turning into C++ for MOCing
file(GLOB tarquingui_UIs ui/*.ui )

# handle the UI files
qt5_wrap_ui( tarquingui_SRCS ${tarquingui_UIs} )

# Resources for this application get turned into C++ files
qt5_add_resources( tarquingui_SRCS ${tarquingui_QRC} )

# moc all the appopriate headers
qt5_wrap_cpp( tarquingui_SRCS ${tarquingui_MOC_HEADERS} )

# Create the executable file (windows RC file is for application icon)
add_executable(tarquingui WIN32 
	${tarquingui_SRCS} 
	${plot_SRCS}
	#	tarquingui.rc
	)

# we need all these libraries
# note that we need ${QT_QTMAIN_LIBRARY} for the entry point
# management on Win32
target_link_libraries(tarquingui 
	tarquincore
	version
	debug ${QWT_LIB_DEBUG} optimized ${QWT_LIB_RELEASE}
	${QT_LIBRARIES} 
	${QT_QTMAIN_LIBRARY}
	${Boost_FILESYSTEM_LIBRARY} 
	${Boost_DATE_TIME_LIBRARY} 
	${Boost_THREAD_LIBRARY} 
	${Boost_CHRONO_LIBRARY} 
	${Boost_SYSTEM_LIBRARY} 
)
