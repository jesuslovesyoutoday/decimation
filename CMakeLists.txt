cmake_minimum_required(VERSION 3.12)

project(decimation)

find_package(VTK REQUIRED)
find_package(xtensor REQUIRED)
find_package(Eigen3 REQUIRED)

if (NOT VTK_FOUND)
  message(FATAL_ERROR "ReadAllUnstructuredGridTypes: Unable to find the VTK build folder.")
endif()

add_executable(decimation main.cpp decimation.cpp)
target_link_libraries(decimation ${VTK_LIBRARIES})
target_link_libraries(decimation Eigen3::Eigen)
