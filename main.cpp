#include <iostream>
#include <string.h>
#ifdef HAVE_BOOST
#include <boost/multiprecision/cpp_int.hpp>
#endif
#include "znbasic.h"
#include "zngroup.h"
#include "zneratosthenes_sieve.h"
#include <sstream>

namespace zn
{
    
    
    template <class T>
    struct bit_handling_t
    {
        static int highest(const T &t)
        {
            return sizeof(T) * 8 - 1;
        }
        static bool isset(const T &t, int i)
        {
            return t & (1 << i) ;
        }
    } ;
    
    template <class N, class E>
    N power(N base, const E &exp)
    {
     /*   if (exp < 0)
            return power(base, -exp) ;*/
        N result = 1 ;
        size_t b = bit_handling_t<E>::highest(exp) ;
        for (size_t i = 0 ; i < b ; i++)
        {
            if (bit_handling_t<E>::isset(exp, i))
                result *= base ;
            base *= base ;
        }
        return result ;
    }
}


using namespace zn ;

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

void test_power(int a, int b, int m)
{
    module_var_t<int, 0>::set(m) ;
    zn_var_t a1(a) ;
    auto pwr = power(a1, b) ;
    std::cout << a1.value() << "^" << b << " = " << pwr.value() 
              << "(mod " << m << ")" << std::endl ;    
}

template <class T>
int test_zn(T a, T b, T m)
{
    typedef zn_t<T, module_var_t<T, 0>> zn ;
    zn za(a), zb(b), zc ;
    module_var_t<T, 0>::set(m) ;

    
    zn s = za + zb ;
    zn d = za - zb ;
    if (s + d != za + za)
        return -1 ;
    std::cout << "za = " << za.value() 
              << ", zb = " << zb.value() 
              << ", s = " << s.value() << std::endl ;
    return 0 ;
}

#ifdef HAVE_BOOST
using namespace boost::multiprecision ;
template <> cpp_int module_var_t<cpp_int, 0>::value_ = cpp_int(0) ;
void test_boost(void)
{
    cpp_int a(2), b(3) ;
    zn_t<cpp_int, module_var_t<cpp_int, 0>> x ;
    module_var_t<cpp_int, 0>::set(b) ;
    std::cout << a * b << std::endl ;
}
#endif

void test_zn_var(int m)
{
    module_var_t<int, 0>::set(m) ;
    std::cout << "Testing for m = " << m << std::endl ;
    zn_var_t a(3), b(5) ;
    std::cout << "a / b = " << (a / b).value() << "(mod " << m << ")" << std::endl ;
}

int main(int argc, char *argv[])
{
    int a = 2 ;
    int b = 7 ; 
    int m = 19 ;
    try
    {
    for (int i = 0 ; i < argc ; i++)
        if (strncmp(argv[i], "--eratosthenes=", 15) == 0)
            test_eratosthenes(atoi(argv[i] + 15)) ;
        else if (strcmp(argv[i], "--zn") == 0)
            test_zn<int>(a, b, m) ;
        else if (strncmp(argv[i], "--zv=", 5) == 0)
            test_zn_var(atoi(argv[i] + 5)) ;
        else if (strcmp(argv[i], "--power") == 0)
            test_power(a, b, m) ;
#ifdef HAVE_BOOST
        else if (strcmp(argv[i], "--boost") == 0)
            test_zn<cpp_int>(a, b, m) ;
#endif
        else if (strncmp(argv[i], "--a=", 4) == 0)
            a = atoi(argv[i] + 4) ;
        else if (strncmp(argv[i], "--b=", 4) == 0)
            b = atoi(argv[i] + 4) ;
        else if (strncmp(argv[i], "--m=", 4) == 0)
            m = atoi(argv[i] + 4) ;
    }
    catch (std::exception &exc)
    {
        std::cout << "Exception: " << exc.what() << std::endl ;
    }
    return 0 ;
}
