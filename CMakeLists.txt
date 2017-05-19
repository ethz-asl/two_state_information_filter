cmake_minimum_required (VERSION 2.8)
project(two_state_implicit_filter)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -march=native")
set(TSIF_VERBOSE 0 CACHE STRING "Verbose flag")
add_definitions(-DTSIF_VERBOSE=${TSIF_VERBOSE})

include_directories(include)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

add_executable(test_state src/test_state.cpp)
add_executable(test_quaternion src/test_quaternion.cpp)
add_executable(test_unit_vector src/test_unit_vector.cpp)
add_executable(test_timeline src/test_timeline.cpp)
add_executable(test_imu_gps src/test_imu_gps.cpp src/simulator.cpp src/camera.cpp)

add_executable(mytest src/test.cpp)

if(ENABLE_OPENCV)
	set(OpenCV_DIR /opt/ros/kinetic/share/OpenCV-3.2.0-dev)
 	find_package(OpenCV REQUIRED)
	add_executable(test_vio src/test_vio.cpp src/simulator.cpp src/camera.cpp)
	target_link_libraries(test_vio ${OpenCV_LIBS})
endif()

