#include <iostream>

namespace zn
{
    template <class T, int N>
    struct module_t
    {
        static T module(void) { return N ; }
    } ;
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
}

using namespace zn ;
int main(void)
{
    typedef module_t<int, 17> mod17_t ;
    typedef zn_t<int, module_t<int, 17> > z17_t ;
    
    z17_t a(15), b(7) ;
    std::cout << (a + b).value() << std::endl ;
}
