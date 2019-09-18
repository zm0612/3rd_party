cmake_minimum_required(VERSION 2.8.3)
project(s2geometry)
add_definitions(--std=c++14 -Wno-unknown-pragmas -w)
find_package(catkin_simple REQUIRED)

catkin_simple(ALL_DEPS_REQUIRED)

cs_add_library(${PROJECT_NAME}_test
	src/mathutil.cc	
	src/s1angle.cc
	src/s2cell_id.cc
	src/s2coords.cc
	src/s2latlng.cc
)
cs_add_executable(s2geometry example/main.cpp)

cs_install()
cs_export()

