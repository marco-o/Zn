#include <iostream>
#include <string.h>
#include "znbasic.h"
#include "zneratosthenes_sieve.h"
#include <sstream>

namespace zn
{
    class exception
    {
    public:
        
    } ;
    
    template <class T>
    class missing_inverse_t : public std::exception
    {
    public:
        missing_inverse_t(const T &lhs, const T &rhs) : lhs_(lhs), rhs_(rhs) {}
        virtual const char *what() const noexcept override
        {
            std::ostringstream ost ;
            ost << lhs_ << " - " << rhs_ ;
            temp_ = ost.str() ;
            return temp_.c_str() ;
        }
    private:
        T lhs_ ;
        T rhs_ ;
        mutable std::string temp_ ;
    } ;
    
    
    template <class T, int N>
    struct module_t
    {
        T module(void) const { return N ; }
    } ;
    
    template <class T, int Tag = 0>
    struct module_var_t
    {
    public:
        static const T &module(void) { return value_ ;}
        static void set(T t) { value_ = t ;}
    private:
        static T value_ ;
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
            return zn_t<T, M>(inverse_value()) ;
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
        zn_t &operator/=(const zn_t<T, M> &rhs)
        {
            value_ = (value * rhs.inverse_value()) % this->module() ;
            return *this ;
        }
    private:
        T inverse_value(void) const
        {
            auto ext = extended_euclidean_algorithm(this->module(), value_) ;
            if (std::get<0>(ext) != 1)
                throw missing_inverse_t<T>(this->module(), value_) ;
            T result = std::get<2>(ext) ;
            if (result < 0)
                result += this->module() ;
            return result ;
        }
        
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


template<>
int module_var_t<int, 0>::value_ = 0 ;

typedef zn_t<int, module_var_t<int, 0> > mod_var_t ;

void test_zn_var(int m)
{
    module_var_t<int, 0>::set(m) ;
    std::cout << "Testing for m = " << m << std::endl ;
    mod_var_t a(3), b(5) ;
    std::cout << "a / b = " << (a / b).value() << "(mod " << m << ")" << std::endl ;
}

int main(int argc, char *argv[])
{

    try
    {
    for (int i = 0 ; i < argc ; i++)
        if (strncmp(argv[i], "--eratosthenes=", 15) == 0)
            test_eratosthenes(atoi(argv[i] + 15)) ;
        else if (strcmp(argv[i], "--zn") == 0)
            test_zn() ;
        else if (strncmp(argv[i], "--zv=", 5) == 0)
            test_zn_var(atoi(argv[i] + 5)) ;
    }
    catch (std::exception &exc)
    {
        std::cout << "Exception: " << exc.what() << std::endl ;
    }
    return 0 ;
}
