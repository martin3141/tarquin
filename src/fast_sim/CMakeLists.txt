PROJECT(fast_sim)

INCLUDE_DIRECTORIES(${CVM_INCLUDE})
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE})

ADD_SUBDIRECTORY(fast_sim)
ADD_SUBDIRECTORY(bin)

# these should be controlled by the top level file
#ADD_DEFINITIONS(-Wall -Wextra -O2 -DNDEBUG)

