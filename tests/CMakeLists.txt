#
# Test executables
#
set (TEST_SOURCES main.cpp)
set (TEST_HEADERS ../include/znbasic.h ../include/zngroup.h)
include_directories(../include)

add_executable(zntest ${TEST_SOURCES} ${TEST_HEADERS})
