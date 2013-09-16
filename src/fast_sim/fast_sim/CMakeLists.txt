project(fast_sim)

include_directories(${FAST_SIM_INCLUDE})

file(GLOB fast_sim_SOURCES *.cpp)

add_library(fast_sim ${fast_sim_SOURCES})

