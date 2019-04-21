#include <iostream>
#include <string.h>
#include "znbasic.h"
#include "zneratosthenes_sieve.h"

namespace zn
{
    template <class T, int N>
    struct module_t
    {
        T module(void) const { return N ; }
    } ;
    
    template <class T, int Tag = 0>
    struct module_var_t
    {
    public:
        const T &module(void) { return value_ ;}
        void set(T t) { value_ = t ;}
    private:
        T value_ ;
    } ;
    //
    // class for modular arithmetic
    //
    template <class T, class M>
    class zn_t : public M // this way can use empty base optimization
    {
    public:
        typedef M module_type ;
        zn_t(const T &value = T()) : value_(value) {}
        T value(void) const { return value_ ;}
        zn_t<T, M> inverse(void) const
        {
            auto ext = extended_euclidean_algorithm(this->module(), value_) ;
            if (std::get<0>(ext) != 1)
                throw std::runtime_error("") ;
            T result = std::get<2>(ext) ;
            if (result < 0)
                result += this->module() ;
            return zn_t<T, M>(result) ;
        }
        zn_t &operator+=(const zn_t<T, M> &rhs)
        {
            value_ = (value + rhs.value_) % this->module() ;
            return *this ;
        }
        zn_t &operator-=(const zn_t<T, M> &rhs)
        {
            value_ = (value - rhs.value_) ;
            if (value_ < 0)
                value += this->module() ;
            return *this ;
        }
        zn_t &operator*=(const zn_t<T, M> &rhs)
        {
            value_ = (value * rhs.value_) % this->module() ;
            return *this ;
        }
    private:
        T  value_ ; 
     } ;
     
    template <class T, class M>
    zn_t<T, M> operator+(const zn_t<T, M> &lhs, const zn_t<T, M> &rhs)
    {
        return zn_t<T, M>((lhs.value() + rhs.value()) % rhs.module()) ;
    }

    template <class T, class M>
    zn_t<T, M> operator-(const zn_t<T, M> &lhs, const zn_t<T, M> &rhs)
    {
        T result = lhs.value() - rhs.value() ;
        if (result < 0)
            result += lhs.module() ;
        return zn_t<T, M>(result) ;
    }
    
    template <class T, class M>
    zn_t<T, M> operator*(const zn_t<T, M> &lhs, const zn_t<T, M> &rhs)
    {
        return zn_t<T, M>((lhs.value() * rhs.value()) % lhs.module()) ;
    }

    template <class T, class M>
    zn_t<T, M> operator/(const zn_t<T, M> &lhs, const zn_t<T, M> &rhs)
    {
        return lhs * rhs.inverse() ;
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

void test_zn(void)
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

int main(int argc, char *argv[])
{

    try
    {
    for (int i = 0 ; i < argc ; i++)
        if (strncmp(argv[i], "--eratosthenes=", 15) == 0)
            test_eratosthenes(atoi(argv[i] + 15)) ;
        else if (strncmp(argv[i], "--zn=", 5) == 0)
            test_zn() ;
    }
    catch (std::exception &exc)
    {
        std::cout << "Exception: " << exc.what() << std::endl ;
    }
    return 0 ;
}
