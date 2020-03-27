#
# Main CMakeLists
#
cmake_minimum_required(VERSION 3.1)

project(zn)

set (CMAKE_CXX_STANDARD 17)
if (WIN32)
	add_definitions(-DNOMINMAX)
endif(WIN32)

add_subdirectory(tests)

