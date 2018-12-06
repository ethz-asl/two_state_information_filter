cmake_minimum_required (VERSION 2.8)
project(two_state_information_filter)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -march=native")
set(TSIF_VERBOSE 0 CACHE STRING "Verbose flag")
add_definitions(-DTSIF_VERBOSE=${TSIF_VERBOSE})

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
find_package(catkin REQUIRED)
find_package(Eigen REQUIRED)

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${EIGEN_INCLUDE_DIR}
)

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

catkin_package(
  INCLUDE_DIRS include ${EIGEN_INCLUDE_DIR}
  LIBRARIES ${PROJECT_NAME}
  CATKIN_DEPENDS
)

add_library(${PROJECT_NAME}
  src/camera.cpp
  src/simulator.cpp
)

target_link_libraries(${PROJECT_NAME} ${catkin_LIBRARIES} pthread)

# Install header
install(DIRECTORY include/tsif/
  DESTINATION ${CATKIN_GLOBAL_INCLUDE_DESTINATION}/tsif/
)

# Install executables and library
install(TARGETS 
          ${PROJECT_NAME}
          test_state 
          test_quaternion 
          test_unit_vector 
          test_timeline 
          test_imu_gps mytest
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)

# Install config 
install(DIRECTORY cfg
  DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}
)
