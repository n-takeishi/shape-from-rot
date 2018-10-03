cmake_minimum_required(VERSION 3.0)

project(sfr)

set(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_BUILD_TYPE Release)

find_package(Boost 1.43 COMPONENTS system program_options)
find_package(OpenCV 3.0.0 REQUIRED)

include_directories(BEFORE "src/3rdparty/")

add_library(sfr STATIC src/sfr.cpp src/PropagateFactor.cpp)
target_link_libraries(sfr ${OpenCV_LIBS} gtsam tbbmalloc tbb boost_system m lapack blas)

add_executable(main src/main.cpp)
target_link_libraries(main sfr boost_program_options)
