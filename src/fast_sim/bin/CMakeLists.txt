project(run_sim)

include_directories(${FAST_SIM_INCLUDE})

add_executable(run_sim main.cpp)

target_link_libraries(run_sim 
	fast_sim
	${CVM_LIB}
	)
