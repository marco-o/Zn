//---------------------------------------------------------------------------------
//
//  Zn 
//  Copyright Marco Oman 2019
//
// Distributed under the Boost Software License, Version 1.0. 
// (See accompanying file LICENSE_1_0.txt or copy at 
// http://www.boost.org/LICENSE_1_0.txt)
//
#include <iostream>
#include <string.h>
#ifdef HAVE_BOOST
#include <boost/multiprecision/cpp_int.hpp>
#endif
#include "zngroup.h"
#include "zneratosthenes_sieve.h"
#include "znquadratic_residue.h"
#include "znquadratic_sieve.h"

#ifdef HAVE_BOOST
using namespace boost::multiprecision;
#endif

using namespace zn ;

#ifdef HAVE_BOOST

template <class large_int, class small_int = int>
void test_quadratic_sieve(const large_int &n, small_int base_size)
{
	auto result = quadratic_sieve(n, base_size);
	std::cout << result << std::endl;
}
#endif

void test_eratosthenes(int bound)
{
    auto p = eratosthenes_sieve<int>(bound) ;
    for (auto k : p)
        std::cout << k << ", " ;
    std::cout << std::endl ;
}

void test_zn1(void)
{
    typedef module_t<int, 17> mod17_t ;
    typedef zn_t<int, module_t<int, 17> > z17_t ;
    z17_t a(15), b(7) ;
    std::cout << (a + b).value() << std::endl ;
    std::cout << (a - b).value() << std::endl ;
    std::cout << (a * b).value() << std::endl ;
    std::cout << (a / b).value() << std::endl ;
    std::cout << a.inverse().value() << ", " << b.inverse().value() << std::endl ;
    std::cout << gcd(27, 84) << std::endl ;
}

template<> int module_var_t<int, 0>::value_ = 0 ;
typedef zn_t<int, module_var_t<int, 0> > zn_var_t ;


#define BOOST_TEST(x)  x

template <class T>
int test_zn(T a, T b, T m)
{
    typedef zn_t<T, module_var_t<T, 0>> zn ;
    zn za(a), zb(b), zc ;
    module_var_t<T, 0>::set(m) ;

    
    zn s = za + zb ;
    zn d = za - zb ;
    bool equal = s + d == 2 * za ;
    BOOST_TEST(equal) ;
    equal = (s - d == 2 * zb) ;
    BOOST_TEST(equal) ;
    equal = (s + d == 2 * a) ;
    BOOST_TEST((za * zb == a * b)) ;
    BOOST_TEST((za + 2) == (a + 2)) ;
    BOOST_TEST((2 + za) == (a + 2)) ;
    BOOST_TEST((2 - za) == (2 - a)) ;
    BOOST_TEST((2 - za) == (2 - a)) ;
    
    
    return 0 ;
}

#ifdef HAVE_BOOST
template <> cpp_int module_var_t<cpp_int, 0>::value_ = cpp_int(0) ;
void test_boost(void)
{
    cpp_int a(2), b(3) ;
    zn_t<cpp_int, module_var_t<cpp_int, 0>> x ;
    module_var_t<cpp_int, 0>::set(b) ;
    std::cout << a * b << std::endl ;
}
#endif

template <class T>
void test_power1(T a, T b, T m)
{
    module_var_t<T, 0>::set(m) ;
    zn_t<T, module_var_t<T, 0> > a1(a) ;
    auto pwr = power(a1, b) ;
    std::cout << a1.value() << "^" << b << " = " << pwr.value() 
              << "(mod " << m << ")" << std::endl ;    
}

void test_power(int a, int b, int m)
{
    test_power1<int>(a, b, m) ;
#ifdef HAVE_BOOST
    test_power1<cpp_int>(a, b, m) ;
#endif
}

void test_zn_var(int m)
{
    module_var_t<int, 0>::set(m) ;
    std::cout << "Testing for m = " << m << std::endl ;
    zn_var_t a(3), b(5) ;
    std::cout << "a / b = " << (a / b).value() << "(mod " << m << ")" << std::endl ;
}

template <class Int>
void test_quadratic_residue(Int a, Int m, Int ps)
{
	Int r = quadratic_residue(a, m, ps);
	std::cout << "Input = " << a << ", " << m << "\n"
			  << " r = " << r << " (" << (r * r) % m << ")" << std::endl;

}

int main(int argc, char *argv[])
{
    int a = 2 ;
    int b = 7 ; 
    int m = 19 ;
	int ps = 1;
	int base_size = 5;
#ifdef HAVE_BOOST
	cpp_int n = 29;
#endif
    try
    {
		for (int i = 0; i < argc; i++)
			if (strncmp(argv[i], "--eratosthenes=", 15) == 0)
				test_eratosthenes(atoi(argv[i] + 15));
			else if (strcmp(argv[i], "--zn") == 0)
			{
				test_zn<int>(a, b, m);
#ifdef HAVE_BOOST
				test_zn<cpp_int>(a, b, m);
#endif
			}
#ifdef HAVE_BOOST
			else if (strncmp(argv[i], "--base-size=", 12) == 0)
				base_size = atoi(argv[i] + 12);
			else if (strncmp(argv[i], "--n=", 4) == 0)
				n = cpp_int(argv[i] + 4);
			else if (strcmp(argv[i], "--qs") == 0)
				test_quadratic_sieve<long long, int>(n.convert_to<int>(), base_size);
				//test_quadratic_sieve<cpp_int>(n, base_size);
#endif
		else if (strncmp(argv[i], "--zv=", 5) == 0)
            test_zn_var(atoi(argv[i] + 5)) ;
        else if (strcmp(argv[i], "--power") == 0)
            test_power(a, b, m) ;
		else if (strcmp(argv[i], "--qr") == 0)
		{
			test_quadratic_residue<long long>(a, m, ps);
#ifdef HAVE_BOOST
			test_quadratic_residue<cpp_int>(a, m, ps);
#endif
		}
#ifdef HAVE_BOOST
        else if (strcmp(argv[i], "--boost") == 0)
            test_zn<cpp_int>(a, b, m) ;
#endif
        else if (strncmp(argv[i], "--a=", 4) == 0)
            a = atoi(argv[i] + 4) ;
        else if (strncmp(argv[i], "--b=", 4) == 0)
            b = atoi(argv[i] + 4) ;
		else if (strncmp(argv[i], "--m=", 4) == 0)
			m = atoi(argv[i] + 4);
		else if (strncmp(argv[i], "--ps=", 5) == 0)
			ps = atoi(argv[i] + 5);
	}
    catch (std::exception &exc)
    {
        std::cout << "Exception: " << exc.what() << std::endl ;
    }
    return 0 ;
}
