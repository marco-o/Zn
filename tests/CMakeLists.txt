#
# Test executables
#

set (COMMON_HEADERS ../include/znbasic.h ../include/zngroup.h ../include/znlinear_solver.h
                  ../include/zneratosthenes_sieve.h ../include/znquadratic_residue.h
				  ../include/znquadratic_sieve_base.h ../include/znquadratic_sieve.h
				  ../include/znpolynomial.h ../include/znmpqs.h ../include/znsiqs.h 
				  ../include/znqueue.h ../include/znelliptic_curve_fact.h)
set (SMALL_HEADERS ../include/basic.h ../include/zneratosthenes_sieve.h 
                   ../include/znelliptic_curve_fact.h)
include_directories(../include)
set (EXECUTABLE_OUTPUT_PATH ../bin)

if (WIN32)
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
endif (WIN32)

if (BOOST_ROOT)
    set(Boost_UNIT_TEST_FRAMEWORK_FOUND on)
	set(BOOST_PATH ${BOOST_ROOT})
    add_definitions(-DHAVE_BOOST)
    include_directories(${BOOST_PATH})
else (BOOST_ROOT)
	find_package(Boost COMPONENTS unit_test_framework)
endif (BOOST_ROOT)


add_executable(main main.cpp logger.cpp ${COMMON_HEADERS})
add_executable(numbox numbox.cpp logger.cpp ${COMMON_HEADERS})
add_executable(small small.cpp logger.cpp ${COMMON_HEADERS})


if (Boost_UNIT_TEST_FRAMEWORK_FOUND)
	if (GMP_ROOT)
		target_link_libraries(main mpir.lib)
	endif (GMP_ROOT)
    if(UNIX)
		target_link_libraries(main pthread)
		target_link_libraries(numbox pthread)
    endif(UNIX)
endif (Boost_UNIT_TEST_FRAMEWORK_FOUND)
