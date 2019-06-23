#
# Test executables
#

set (COMMON_HEADERS ../include/znbasic.h ../include/zngroup.h ../include/znlinear_solver.h
                  ../include/zneratosthenes_sieve.h ../include/znquadratic_residue.h
				  ../include/znquadratic_sieve0.h ../include/znquadratic_sieve.h
				  ../include/znquadratic_sieve_base.h  ../include/znpolynomial.h
				  ../include/znmpqs.h ../include/znsiqs.h ../include/znqueue.h)
include_directories(../include)

if (WIN32)
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
endif (WIN32)

if (BOOST_ROOT)
    set(Boost_UNIT_TEST_FRAMEWORK_FOUND on)
	set(BOOST_PATH ${BOOST_ROOT})
else (BOOST_ROOT)
	find_package(Boost COMPONENTS unit_test_framework)
endif (BOOST_ROOT)

if (Boost_UNIT_TEST_FRAMEWORK_FOUND)
    include_directories(${BOOST_PATH})
    add_definitions(-DHAVE_BOOST)
    add_executable(zntest zntest.cpp logger.cpp ${COMMON_HEADERS})
    target_link_libraries(zntest ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
	add_executable(numbox logger.cpp numbox.cpp ${COMMON_HEADERS})
	if (GMP_ROOT)
		add_definitions(-DHAVE_GMP)
		include_directories(${GMP_ROOT})
		link_directories(${GMP_ROOT})
	endif (GMP_ROOT)
else (Boost_UNIT_TEST_FRAMEWORK_FOUND)
    message("BOOST not found")
endif (Boost_UNIT_TEST_FRAMEWORK_FOUND)

add_executable(main main.cpp logger.cpp ${COMMON_HEADERS})

if (Boost_UNIT_TEST_FRAMEWORK_FOUND)
	if (GMP_ROOT)
		target_link_libraries(main mpir.lib)
	endif (GMP_ROOT)
        if(UNIX)
		target_link_libraries(main pthread)
		target_link_libraries(numbox pthread)
        endif(UNIX)
endif (Boost_UNIT_TEST_FRAMEWORK_FOUND)
