cmake_minimum_required (VERSION 2.8)
project(generalized_information_filter)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -march=native")
set(GTEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/googletest CACHE STRING "gtest directory")

##################### Find/Include #####################
find_package(kindr REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(glog REQUIRED)

include_directories(include)
include_directories(${kindr_INCLUDE_DIRS})
include_directories(${EIGEN3_INCLUDE_DIR})

##################### Catkin #####################
if(DEFINED CATKIN_DEVEL_PREFIX)
	find_package(catkin REQUIRED COMPONENTS
	)
	catkin_package(
		INCLUDE_DIRS include ${catkin_INCLUDE_DIRS}
	    LIBRARIES ${PROJECT_NAME}
		CATKIN_DEPENDS
	)
endif()

##################### Create library #####################
add_library(${PROJECT_NAME}
	src/measurement.cc
	src/element-vector.cc
	src/element-vector-definition.cc
)
target_link_libraries(${PROJECT_NAME} pthread glog)

##################### Check for gtest #####################
if(EXISTS "${GTEST_DIR}" AND NOT DEFINED CATKIN_DEVEL_PREFIX)
	message(STATUS "Building GTests!")
	option(BUILD_GTEST "build gtest" ON)
	add_subdirectory(${GTEST_DIR} gtest_custom)
	enable_testing()
	include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
	
	add_executable(test_simple test/test-simple.cc)
	target_link_libraries(test_simple ${PROJECT_NAME} gtest_main gtest)
	add_test(test_simple test_simple)
endif()

##################### Cmake find generation #####################
# Generate FindGIF.cmake file
file(WRITE cmake/FindGIF.cmake
"# This file was automatically generated during the installation of the genralized information filter library
# and can be used through cmake to find the corresponding header files. A copy of this
# file was created in ${CMAKE_ROOT}/Modules.

set(GIF_INCLUDE_DIRS
${CMAKE_INSTALL_PREFIX}/include/GIF/include
)
set(GIF_FOUND TRUE)
message(\"-- GIF found (include: ${CMAKE_INSTALL_PREFIX}/include/GIF/include)\")
"
)

##################### Install stuff if required #####################
# Setting for make install
install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_INSTALL_PREFIX}/include/GIF)")
install(DIRECTORY include DESTINATION ${CMAKE_INSTALL_PREFIX}/include/GIF)
install(FILES cmake/FindGIF.cmake DESTINATION ${CMAKE_ROOT}/Modules)
