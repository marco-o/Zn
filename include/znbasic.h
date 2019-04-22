#ifndef znbasic_H
#define znbasic_H

#include <tuple>
#include <sstream>


namespace zn
{
    template <class T>
    class not_relatively_prime_t : public std::exception
    {
    public:
        not_relatively_prime_t(const T &lhs, const T &rhs) : lhs_(lhs), rhs_(rhs) {}
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
            r[i] = r[i] - q * r[j] ;
            s[i] = s[i] - q * s[j] ;
            t[i] = t[i] - q * t[j] ;
        }
        return std::tuple<T, T, T>(r[i], s[i], t[i]) ;
    }

}
#endif