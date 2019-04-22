cmake_minimum_required(VERSION 3.1)

set (CMAKE_CXX_STANDARD 11)

if (BOOST_PATH)
    include_directories(${BOOST_PATH})
    add_definitions(-DHAVE_BOOST)
else (BOOST_PATH)
    message("Define BOOST_PATH") 
endif (BOOST_PATH)

project(zn)


add_subdirectory(tests)

#install(TARGETS myapp DESTINATION bin)
