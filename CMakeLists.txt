cmake_minimum_required(VERSION 2.8)
project(akima_spline)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
add_executable(akima_spline ./main.cpp)

install(TARGETS akima_spline
      DESTINATION ${PROJECT_SOURCE_DIR}/bin)

