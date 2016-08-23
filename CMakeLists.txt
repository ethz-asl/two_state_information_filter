cmake_minimum_required (VERSION 2.8)
project(generalized_information_filter)
add_definitions(-std=c++11)

##################### Find, include, and compile library #####################

find_package(kindr REQUIRED)
find_package(Eigen3 REQUIRED)

include_directories(${kindr_INCLUDE_DIRS})
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(include)

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/gtest/" AND (NOT DEFINED CATKIN_DEVEL_PREFIX))
	message(STATUS "Building GTests!")
	option(BUILD_GTEST "build gtest" ON)
	add_subdirectory(gtest gtest)
	enable_testing()
	include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
	
	add_library(${PROJECT_NAME}  
	 src/measurement.cc 
	 src/state.cc
	 src/state-definition.cc 
  )
	target_link_libraries(GIF pthread)
	
	
	add_executable(test_simple test/test-simple.cc)
	target_link_libraries(test_simple ${PROJECT_NAME} gtest_main gtest pthread)
	add_test(test_simple test_simple)
endif()

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

# Setting for make install
install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_INSTALL_PREFIX}/include/GIF)")
install(DIRECTORY include DESTINATION ${CMAKE_INSTALL_PREFIX}/include/GIF)
install(FILES cmake/FindGIF.cmake DESTINATION ${CMAKE_ROOT}/Modules)
