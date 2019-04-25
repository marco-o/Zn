cmake_minimum_required(VERSION 3.1)

set (CMAKE_CXX_STANDARD 11)
if (WIN32)
	add_definitions(-DNOMINMAX)
endif(WIN32)

add_subdirectory(tests)

