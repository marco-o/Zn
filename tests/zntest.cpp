//#define BOOST_TEST_MODULE zntest
#include <boost/test/included/unit_test.hpp>
#include "zngroup.h"
#include <boost/multiprecision/cpp_int.hpp>



using namespace boost::unit_test;
using namespace zn ;

template <> int zn::module_var_t<int, 0>::value_ = 0 ;
using namespace boost::multiprecision ;
template <> cpp_int module_var_t<cpp_int, 0>::value_ = cpp_int(0) ;

template <class T>
void test_zn(T a, T b, T m)
{
    typedef zn_t<T, module_var_t<T, 0>> zn ;
    zn za(a), zb(b), zc ;
    module_var_t<T, 0>::set(m) ;

    
    zn s = za + zb ;
    zn d = za - zb ;
    bool equal = (s + d == za + za) ; 
    BOOST_TEST(equal) ;
}

void free_test_function()
{
    test_zn<int>(10, 15, 31) ;
    test_zn<cpp_int>(10, 15, 31) ;
}

test_suite* init_unit_test_suite( int /*argc*/, char* /*argv*/[] )
{
    framework::master_test_suite().add( BOOST_TEST_CASE( &free_test_function ) );
    return 0;
}

