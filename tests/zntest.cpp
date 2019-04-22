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
int test_zn(T a, T b, T m)
{
    typedef zn_t<T, module_var_t<T, 0>> zn ;
    module_var_t<T, 0>::set(m) ;
    zn za(a), zb(b), zc ;

    
    zn s = za + zb ;
    zn d = za - zb ;
    bool equal = s + d == 2 * za ;
    BOOST_TEST(equal) ;
    equal = (s - d == 2 * zb) ;
    BOOST_TEST(equal) ;
    equal = (s + d == 2 * a) ;
    BOOST_TEST(equal) ;
    equal = ((za + 2) == (a + 2)) ;
    BOOST_TEST(equal) ;
    equal = ((2 + za) == (a + 2)) ;
    BOOST_TEST(equal) ;
    equal = ((2 - za) == (2 - a)) ;
    BOOST_TEST(equal) ;
    equal = ((za - 2) == (a - 2)) ;
    BOOST_TEST(equal) ;
    equal = (za * zb == a * b) ;
    BOOST_TEST(equal) ;
    equal = (za * b == a * b) ;
    BOOST_TEST(equal) ;
    equal = (a * zb == a * b) ;
    BOOST_TEST(equal) ;
    equal = b * (za / zb) == a  ;
    BOOST_TEST(equal) ;
    
    
    return 0 ;
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

