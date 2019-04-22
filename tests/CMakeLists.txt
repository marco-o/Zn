#
# Test executables
#
set (TEST_SOURCES main.cpp)
set (TEST_HEADERS ../include/znbasic.h ../include/zngroup.h)
include_directories(../include)

find_package(Boost COMPONENTS unit_test_framework)

#find_package(Boost COMPONENTS unit_test_framework REQUIRED)

if (Boost_UNIT_TEST_FRAMEWORK_FOUND)
    include_directories(${BOOST_PATH})
    add_definitions(-DHAVE_BOOST)
    add_executable(zntest zntest.cpp)
    target_link_libraries(zntest ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
else (Boost_UNIT_TEST_FRAMEWORK_FOUND)
    message("BOOST not found")
endif (Boost_UNIT_TEST_FRAMEWORK_FOUND)

project(zn)


add_executable(main ${TEST_SOURCES} ${TEST_HEADERS})

