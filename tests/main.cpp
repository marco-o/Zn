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
#include <fstream>
#include <string.h>
#ifdef HAVE_BOOST
#include <boost/multiprecision/cpp_int.hpp>
#ifdef HAVE_GMP
#include <boost/multiprecision/gmp.hpp>
#endif
#endif
#include "zngroup.h"
#include "znelliptic_curve_fact.h"
#include "zneratosthenes_sieve.h"
#include "znquadratic_residue.h"
#include "znquadratic_sieve0.h"
#include "znquadratic_sieve.h"
#include "znmpqs.h"
#include "znsiqs.h"

#ifdef HAVE_BOOST
using namespace boost::multiprecision;
#endif

using namespace zn ;


template <class large_int, class small_int = int>
void test_quadratic_sieve0(const large_int &n, small_int base_size)
{
	auto p1 = quadratic_sieve0(n, base_size);
	auto p2 = n / p1;
	std::cout << p1 << " * " << p2 << " = " << n << std::endl;
}

template <class large_int, class small_int = int>
void test_quadratic_sieve(const large_int &n, small_int base_size)
{
	auto p1 = quadratic_sieve(n, base_size);
	auto p2 = n / p1;
	std::cout << p1 << " * " << p2 << " = " << n << std::endl;
}

template <class large_int, class small_int = int, class real = float>
void test_multiple_polynomial_quadratic_sieve(const large_int &n, const large_int &m, small_int base_size, int k = 0)
{
	auto p1 = multiple_polynomial_quadratic_sieve<large_int, small_int, real>(n, m, base_size, k);
	auto p2 = n / p1;
	std::cout << p1 << " * " << p2 << " = " << n << std::endl;
}

template <class large_int, class small_int = int, class real = float>
void test_self_initializing_quadratic_sieve(const large_int &n, const large_int &m, small_int base_size, int order = 0)
{
	auto p1 = self_initializing_quadratic_sieve<large_int, small_int, real>(n, m, base_size, order);
	auto p2 = n / p1;
	std::cout << p1 << " * " << p2 << " = " << n << std::endl;
}

template <class large_int>
void test_pollard_rho(const large_int &n, int count)
{
	large_int p1 = pollards_rho(n, count);
	auto p2 = n / p1;
	std::cout << p1 << " * " << p2 << " = " << n << std::endl;
	exit(0);
}

template <class large_int, class small_int = long>
void test_pollard_p1(const large_int &n, int count)
{
	pollard_p1_t<large_int, small_int> pollard_p1;
	pollard_p1.init(safe_cast<small_int>(n));

	large_int p1 = pollard_p1.fact(n, count);
	auto p2 = n / p1;
	std::cout << p1 << " * " << p2 << " = " << n << std::endl;
	exit(0);
}


template <class large_int>
void test_elliptic_curve_homo(const large_int &n)
{
	typedef typename elliptic_curve_projective_t<large_int>::point_t point_t;
	point_t pt{ 1, 1 };
	elliptic_curve_projective_t<large_int> ec(n);
	ec.init(6, pt);
#ifdef _DEBUG
	ec.test(pt);
#endif
}

template <class large_int, class small_int>
void test_elliptic_curve(const large_int &n, small_int range)
{
	typedef typename elliptic_curve_t<large_int>::point_t point_t;
	point_t pt{ 1, 1 };
	elliptic_curve_t<large_int> ec(n);
	auto primes = eratosthenes_sieve<small_int>(static_cast<int>(range));
	small_int exp = 1;
	small_int limit = static_cast<small_int>(sqrt(std::numeric_limits<small_int>::max()));
	for (small_int p = 2 ; p < range ; p++)
	{
		small_int exp1 = exp * p;
		if (exp1 > limit)
		{
			if (ec.run(5, pt, exp))
				break;
			exp = p;
		}
		else
			exp = exp1;
	}
	std::cout << "Factor = " << ec.factor() << std::endl;
}

void test_eratosthenes(int bound)
{
    auto p = eratosthenes_sieve<int>(bound) ;
    for (auto k : p)
        std::cout << k << ", " ;
    std::cout << std::endl ;
}

template <class large_int, class small_int = int>
void test_polynomial_generation(const large_int n, small_int m, small_int base_size)
{
	small_int range = quadratic_sieve_base_t<large_int, small_int, short>::primes_range(base_size * 2);
	auto primes = eratosthenes_sieve<small_int>(static_cast<int>(range));
	std::vector<prime_info_t<small_int>> residual_primes;
	small_int limit = static_cast<small_int>(std::sqrt(std::numeric_limits<small_int>::max()));
	for (auto p : primes)
	{
		if (p >= limit)
			break;
		small_int p2 = p * p;
		small_int n1 = safe_cast<small_int>(n % p2);
		small_int residue = quadratic_residue<small_int>(n1, p2, p); // actually a power of prime
		if (residue > 0)
			residual_primes.push_back(prime_info_t<small_int>(p, p2, residue));
	}
	polynomial_generator_t<large_int, small_int> generator(n, m, residual_primes);
	generator.order_init(8);
	for (int i = 0; i < 100; i++)
	{
		polynomial_seed_t seed = generator();
		polynomial_siqs_t<large_int, small_int> poly(seed.index, residual_primes, n);
		std::sort(poly.b1.begin(), poly.b1.end());
		auto g1 = residual_primes[seed.index[3]];
		auto n2 = safe_cast<small_int>(n % g1.prime);
		std::cout << "n2 = " << n2 << std::endl;
		std::cout << "a = " << poly.a << "\n";
		for (auto &b1 : poly.b1)
		{
			auto b2 = safe_cast<small_int>(b1 % g1.prime);
			std::cout << "b = " << b1 << ", b^2 = " << (b2 * b2) % g1.prime << std::endl;
		}
	}
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
	const char *count = "100";
	long long a = 2 ;
	long long b = 7 ;
	long long m = 160000;
	long long ps = 1;
	long long base_size = 3900;
	int k = 0;
#ifdef HAVE_BOOST
	const char *n = "43169554144061480807721762059907068496313438381696909238551841";
	const char *a1 = "2";
	const char *m1 = "160000";
	const char *b1 = "7";
	const char *ps1 = "1";
#endif

    try
    {
		log_base_t::init(argc, argv);
		if (argc == 1)
			test_multiple_polynomial_quadratic_sieve<cpp_int, long long, short>(cpp_int(n), cpp_int(m1), base_size);
		for (int i = 0; i < argc; i++)
			if (strncmp(argv[i], "--eratosthenes=", 15) == 0)
				test_eratosthenes(atoi(argv[i] + 15));
		else if (strcmp(argv[i], "--zn") == 0)
		{
			test_zn<int>(static_cast<int>(a), static_cast<int>(b), static_cast<int>(m));
#ifdef HAVE_BOOST
			test_zn<cpp_int>(a, b, m);
#endif
		}
		else if (strncmp(argv[i], "--count=", 8) == 0)
			count = argv[i] + 8;
#ifdef HAVE_BOOST
		else if (strncmp(argv[i], "--base-size=", 12) == 0)
			base_size = atoi(argv[i] + 12);
		else if (strncmp(argv[i], "--n=", 4) == 0)
			n = argv[i] + 4;
		else if (strncmp(argv[i], "--k=", 4) == 0)
			k = atoi(argv[i] + 4);
		else if (strcmp(argv[i], "--qs0") == 0)
			test_quadratic_sieve0<long long, int>(atoll(n), static_cast<int>(base_size));
		else if (strcmp(argv[i], "--qs") == 0)
			test_quadratic_sieve<long long, long long>(atoll(n), static_cast<int>(base_size));
		else if (strcmp(argv[i], "--qsc0") == 0)
			test_quadratic_sieve0<cpp_int, long long>(cpp_int(n), base_size);
		else if (strcmp(argv[i], "--qsc") == 0)
			test_quadratic_sieve<cpp_int, long long>(cpp_int(n), base_size);
		else if (strcmp(argv[i], "--mpqs") == 0)
			test_multiple_polynomial_quadratic_sieve<cpp_int, long long, short>(cpp_int(n), cpp_int(m1), base_size, k);
		else if (strcmp(argv[i], "--mpqsl") == 0)
			test_self_initializing_quadratic_sieve<long long, long long>(atoll(n), atoll(m1), base_size);
		else if (strcmp(argv[i], "--siqs") == 0)
			test_self_initializing_quadratic_sieve<cpp_int, long long, unsigned char>(cpp_int(n), cpp_int(m1), base_size, k);
		else if (strcmp(argv[i], "--siqsl") == 0)
			test_multiple_polynomial_quadratic_sieve<long long, long long>(atoll(n), atoll(m1), base_size);
		else if (strcmp(argv[i], "--polytest") == 0)
			test_polynomial_generation<cpp_int, long long>(cpp_int(n), atoll(m1), base_size);
#ifdef HAVE_GMP
		else if (strcmp(argv[i], "--qsg") == 0)
			test_quadratic_sieve<mpz_int, long long>(mpz_int(n), base_size);
		else if (strcmp(argv[i], "--mpqsg") == 0)
			test_multiple_polynomial_quadratic_sieve<mpz_int, long long, short>(mpz_int(n), mpz_int(m1), base_size);
#endif
#endif
		else if (strcmp(argv[i], "--ec") == 0)
			test_elliptic_curve<cpp_int, long long>(cpp_int(n), base_size);
		else if (strcmp(argv[i], "--ecl") == 0)
			test_elliptic_curve<long long, int>(atoll(n), static_cast<int>(base_size));
		else if (strcmp(argv[i], "--ech") == 0)
			test_elliptic_curve_homo<long long>(atoll(n));
		else if (strcmp(argv[i], "--rhol") == 0)
			test_pollard_rho(atoll(n), atoi(count));
		else if (strcmp(argv[i], "--rho") == 0)
			test_pollard_rho(cpp_int(n), atoi(count));
		else if (strcmp(argv[i], "--p1l") == 0)
			test_pollard_p1<long long, long>(atoll(n), atoi(count));
		else if (strcmp(argv[i], "--p1") == 0)
			test_pollard_p1<cpp_int, long long>(cpp_int(n), atoi(count));
		else if (strncmp(argv[i], "--zv=", 5) == 0)
			test_zn_var(atoi(argv[i] + 5)) ;
		else if (strcmp(argv[i], "--power") == 0)
			test_power(static_cast<int>(a), static_cast<int>(b), static_cast<int>(m)) ;
		else if (strcmp(argv[i], "--qr") == 0)
			test_quadratic_residue<long long>(a, m, ps);
#ifdef HAVE_BOOST
		else if (strcmp(argv[i], "--qrc") == 0)
			test_quadratic_residue<cpp_int>(cpp_int(a1), cpp_int(m1), cpp_int(ps1));
#endif
#ifdef HAVE_BOOST
		else if (strcmp(argv[i], "--boost") == 0)
			test_zn<cpp_int>(a, b, m) ;
#endif
		else if (strncmp(argv[i], "--a=", 4) == 0)
			a = atoll(a1 = argv[i] + 4) ;
		else if (strncmp(argv[i], "--b=", 4) == 0)
			b = atoll(b1 = argv[i] + 4) ;
		else if (strncmp(argv[i], "--m=", 4) == 0)
			m = atoll(m1 = argv[i] + 4);
		else if (strncmp(argv[i], "--ps=", 5) == 0)
			ps = atoll(ps1 = argv[i] + 5);
	}
    catch (std::exception &exc)
    {
        std::cout << "Exception: " << exc.what() << std::endl ;
    }
    return 0 ;
}
