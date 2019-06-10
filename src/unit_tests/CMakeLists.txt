project(tarquin_tests)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

# glob all the test sources
file(GLOB tarquin_tests_SOURCES *.cpp)

# glob all the test data
file(GLOB tarquin_tests_DATA data/*.dpt data/*.SDAT data/*.SPAR data/*.ref data/*.ima data/*.rda data/*.P data/*.dcm data/*.xml data/NAA.csv)

include_directories(${TARQUIN_INCLUDE})
include_directories(${CVM_INCLUDE})
include_directories(${LEVMAR_INCLUDE})
include_directories(${FAST_SIM_INCLUDE})
include_directories(${EIGEN3_INCLUDE})

add_executable(tarquin_tests ${tarquin_tests_SOURCES})

target_link_libraries(tarquin_tests 
	tarquincore
	${Boost_FILESYSTEM_LIBRARY} 
	${Boost_DATE_TIME_LIBRARY} 
	${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
	${Boost_THREAD_LIBRARY} 
	${Boost_CHRONO_LIBRARY} 
	${Boost_SYSTEM_LIBRARY} 
	)

add_test(tarquin_tests "tarquin_tests") 

message( status "the redist root is ${REDIST_TREE_ROOT}" )

# copy some files to the testing directory
FILE(COPY ${tarquin_tests_DATA} DESTINATION ${REDIST_TREE_ROOT}/data)
