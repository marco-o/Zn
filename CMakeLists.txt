cmake_minimum_required(VERSION 2.8)

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")

if (BOOST_PATH)
    include_directories(${BOOST_PATH})
    add_definitions(-DHAVE_BOOST)
else (BOOST_PATH)
    message("Define BOOST_PATH") 
endif (BOOST_PATH)

project(app_project)

add_executable(myapp main.cpp)

#install(TARGETS myapp DESTINATION bin)
