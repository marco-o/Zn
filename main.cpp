#include <iostream>
#include <tuple>

namespace zn
{
    template <class T, int N>
    struct module_t
    {
        static T module(void) { return N ; }
    } ;
    
    template <class T>
    T gcd(T a, T b)
    {
        while (a)
        {
            T c = b % a ;
            b = a ;
            a = c ;
        }
        return b ;
    }
    
    template <class T>
    std::tuple<T, T, T> extended_euclidean_algorithm(const T &a, const T &b)
    {
        T q ;
        int i, j ;
        T r[2] = {a, b} ;
        T s[2] = {1, 0} ;
        T t[2] = {0, 1} ;
        
        for (i = 0, j = 1 ; r[j] ; i=j, j = (j + 1) & 1)
        {
            q = r[i] / r[j] ;
            r[i] = r[j] - q * r[i] ;
            s[i] = s[j] - q * s[i] ;
            t[i] = t[j] - q * t[i] ;
        }
        return std::tuple<T, T, T>(r[i], s[i], t[i]) ;
    }
    //
    // class for modular arithmetic
    //
    template <class T, class M>
    class zn_t
    {
    public:
        typedef M module_type ;
        zn_t(const T &value = T()) : value_(value) {}
        T value(void) const { return value_ ;}
        zn_t &operator+=(const zn_t<T, M> &rhs)
        {
            value_ = (value + rhs.value_) % module_type::module() ;
            return *this ;
        }
        zn_t &operator-=(const zn_t<T, M> &rhs)
        {
            value_ = (value - rhs.value_) ;
            if (value_ < 0)
                value += module_type::module() ;
            return *this ;
        }
        zn_t &operator*=(const zn_t<T, M> &rhs)
        {
            value_ = (value * rhs.value_) % module_type::module() ;
            return *this ;
        }
    private:
        T  value_ ; 
     } ;
     
    template <class T, class M>
    zn_t<T, M> operator+(const zn_t<T, M> &lhs, const zn_t<T, M> &rhs)
    {
        return zn_t<T, M>((lhs.value() + rhs.value()) % M::module()) ;
    }

    template <class T, class M>
    zn_t<T, M> operator-(const zn_t<T, M> &lhs, const zn_t<T, M> &rhs)
    {
        T result = lhs.value() - rhs.value() ;
        if (result < 0)
            result += M::module() ;
        return zn_t<T, M>(result) ;
    }
    
    template <class T, class M>
    zn_t<T, M> operator*(const zn_t<T, M> &lhs, const zn_t<T, M> &rhs)
    {
        return zn_t<T, M>((lhs.value() * rhs.value()) % M::module()) ;
    }

    template <class T, class M>
    zn_t<T, M> operator/(const zn_t<T, M> &lhs, const zn_t<T, M> &rhs)
    {
        auto result = extended_euclidean_algorithm(M::module(), rhs.value()) ;
        if (std::get<0>(result) != 1)
            throw std::runtime_error("") ;
        return zn_t<T, M>((lhs.value() * std::get<1>(result)) % M::module()) ;
    }
}

using namespace zn ;
int main(void)
{
    typedef module_t<int, 17> mod17_t ;
    typedef zn_t<int, module_t<int, 17> > z17_t ;
    
    z17_t a(15), b(7) ;
    std::cout << (a + b).value() << std::endl ;
    std::cout << (a - b).value() << std::endl ;
    std::cout << (a * b).value() << std::endl ;
    std::cout << (a / b).value() << std::endl ;
    std::cout << gcd(27, 84) << std::endl ;
    return 0 ;
}
