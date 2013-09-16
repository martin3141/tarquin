project(ftn)

# Why won't FORTRAN go away?
enable_language(Fortran)

# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# All files in this directory to be compiled into the library
file(GLOB ftn_SOURCES *.f)

# This is the library
add_library(ftn ${ftn_SOURCES})


