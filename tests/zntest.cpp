//#define BOOST_TEST_MODULE zntest
#include <boost/test/included/unit_test.hpp>
#include "zngroup.h"
#include "zneratosthenes_sieve.h"
#include <boost/multiprecision/cpp_int.hpp>
#include "znquadratic_residue.h"



using namespace boost::unit_test;
using namespace zn ;

template <> int zn::module_var_t<int, 0>::value_ = 0 ;
using namespace boost::multiprecision ;
template <> cpp_int module_var_t<cpp_int, 0>::value_ = cpp_int(0) ;

template <class T>
int test_zngroup_1(T a, T b, T m)
{
    module_var_t<T, 0>::set(m) ;
    zn_t<T, module_var_t<T, 0>> za(a), zb(b) ;

    bool equal = za + zb == a + b ;
    BOOST_TEST(equal) ;
    equal = ((za + b) == (a + b)) ;
    BOOST_TEST(equal) ;
    equal = ((b + za) == (a + b)) ;
    BOOST_TEST(equal) ;
    equal = ((b - za) == (b - a)) ;
    BOOST_TEST(equal) ;
    equal = ((za - b) == (a - b)) ;
    BOOST_TEST(equal) ;
    equal = (za * zb == a * b) ;
    BOOST_TEST(equal) ;
    equal = (za * b == a * b) ;
    BOOST_TEST(equal) ;
    equal = (a * zb == a * b) ;
    BOOST_TEST(equal) ;
    equal = b * (za / zb) == a  ;
    BOOST_TEST(equal) ;
    equal = b * (za / b) == a  ;
    BOOST_TEST(equal) ;
    equal = b * (a / zb) == a  ;
    BOOST_TEST(equal) ;
    
    
    return 0 ;
}

template<class T>
void test_zngroup()
{
    test_zngroup_1<T>(10, 15, 31) ;
}

template <class T>
void test_eratosthenes(void)
{
    auto primes = eratosthenes_sieve<T>(100) ;
    BOOST_TEST(primes.size() == 25) ;
    BOOST_TEST(*primes.rbegin() == 97) ;
}

template <class T>
void test_quadratic_residue(void)
{
	auto primes = eratosthenes_sieve<T>(40000);
	primes.erase(primes.begin(), primes.begin() + primes.size() * 2 / 3);
	for (auto p : primes)
		for (int i = 0; i < 100; i++)
		{
			T r = (rand() % (p - 1)) + 1;
			T q = quadratic_residue(r, p);
			T r1 = (q * q) % p;
			if (q != 0)
				BOOST_TEST(r1 == r);
			else
				for (i = 1; i < p; i++)
					if ((i * i) % p == r)
						BOOST_TEST(false);
		}
}

test_suite* init_unit_test_suite( int /*argc*/, char* /*argv*/[] )
{
    framework::master_test_suite().add(BOOST_TEST_CASE( &test_zngroup<int>));
    framework::master_test_suite().add(BOOST_TEST_CASE( &test_zngroup<cpp_int>));
    framework::master_test_suite().add(BOOST_TEST_CASE( &test_eratosthenes<int>));
	framework::master_test_suite().add(BOOST_TEST_CASE(&test_eratosthenes<cpp_int>));
	framework::master_test_suite().add(BOOST_TEST_CASE(&test_quadratic_residue<int>));
	framework::master_test_suite().add(BOOST_TEST_CASE(&test_quadratic_residue<cpp_int>));
	return 0;
}

