//---------------------------------------------------------------------------------
//
//  Zn 
//  Copyright Marco Oman 2019
//
// Distributed under the Boost Software License, Version 1.0. 
// (See accompanying file LICENSE_1_0.txt or copy at 
// http://www.boost.org/LICENSE_1_0.txt)
//
#ifndef zngroup_H
#define zngroup_H

#include "znbasic.h"

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
        zn_t(const T &value = T())
        {
            value_ = value % this->module() ;
            if (value_ < 0)
                value_ += this->module() ;
        }
        T value(void) const { return value_ ;}
        explicit operator T (void) const { return value_ ;}
        zn_t<T, M> inverse(void) const
        {
            return zn_t<T, M>(inverse_value()) ;
        }
        zn_t operator-(void) const 
        {
            return zn_t<T, M>(this->module() - value_) ;
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
            value_ = (value_ * rhs.value_) % this->module() ;
            return *this ;
        }
        zn_t &operator/=(const zn_t<T, M> &rhs)
        {
            value_ = (value_ * rhs.inverse_value()) % this->module() ;
            return *this ;
        }
    private:
        T inverse_value(void) const
        {
            auto ext = extended_euclidean_algorithm(this->module(), value_) ;
            if (std::get<0>(ext) != 1)
                throw not_relatively_prime_t<T>(this->module(), value_) ;
            T result = std::get<2>(ext) ;
            if (result < 0)
                result += this->module() ;
            return result ;
        }
        
        T  value_ ; // assumed in the range [0; module()[
     } ;
     
    template <class T, class M>
    zn_t<T, M> operator+(const zn_t<T, M> &lhs, const zn_t<T, M> &rhs)
    {
        return zn_t<T, M>((lhs.value() + rhs.value()) % rhs.module()) ;
    }

    template <class T, class M, class U>
    zn_t<T, M> operator+(const zn_t<T, M> &lhs, const U &rhs)
    {
        return zn_t<T, M>((lhs.value() + rhs) % lhs.module()) ;
    }
    
    template <class T, class M, class U>
    zn_t<T, M> operator+(const U &lhs, const zn_t<T, M> &rhs)
    {
        return zn_t<T, M>((lhs + rhs.value()) % rhs.module()) ;
    }

    
    template <class T, class M>
    zn_t<T, M> operator-(const zn_t<T, M> &lhs, const zn_t<T, M> &rhs)
    {
        T result = (lhs.value() - rhs.value()) % lhs.module() ;
        if (result < 0)
            result += lhs.module() ;
        return zn_t<T, M>(result) ;
    }
    
    template <class T, class M, class U>
    zn_t<T, M> operator-(const zn_t<T, M> &lhs, const U &rhs)
    {
        T result = (lhs.value() - rhs) % lhs.module() ;
        if (result < 0)
            result += lhs.module() ;
        return zn_t<T, M>(result) ;
    }
    
    template <class T, class M, class U>
    zn_t<T, M> operator-(const U &lhs, const zn_t<T, M> &rhs)
    {
        T result = (lhs - rhs.value()) % rhs.module() ;
        if (result < 0)
            result += rhs.module() ;
        return zn_t<T, M>(result) ;
    }
    
    template <class T, class M>
    zn_t<T, M> operator*(const zn_t<T, M> &lhs, const zn_t<T, M> &rhs)
    {
        return zn_t<T, M>((lhs.value() * rhs.value()) % lhs.module()) ;
    }

    template <class T, class M, class U>
    zn_t<T, M> operator*(const zn_t<T, M> &lhs, const U &rhs)
    {
        return zn_t<T, M>((lhs.value() * rhs) % lhs.module()) ;
    }
    
    template <class T, class M, class U>
    zn_t<T, M> operator*(const U &lhs, const zn_t<T, M> &rhs)
    {
        return zn_t<T, M>((lhs * rhs.value()) % rhs.module()) ;
    }

    template <class T, class M>
    zn_t<T, M> operator/(const zn_t<T, M> &lhs, const zn_t<T, M> &rhs)
    {
        return lhs * rhs.inverse() ;
    }
    
    template <class T, class M, class U>
    zn_t<T, M> operator/(const U &lhs, const zn_t<T, M> &rhs)
    {
        return lhs * rhs.inverse() ;
    }

    template <class T, class M, class U>
    zn_t<T, M> operator/(const zn_t<T, M> &lhs, const U &rhs)
    {
        return lhs * zn_t<T, M>(rhs).inverse() ;
    }

    template <class T, class M>
    bool operator==(const zn_t<T, M> &lhs, const zn_t<T, M> &rhs)
    {
        return lhs.value() == rhs.value();
    }

    template <class T, class M, class U>
    bool operator==(const zn_t<T, M> &lhs, const U &rhs)
    {
        return (lhs.value() - rhs) % lhs.module() == 0 ;
    }

    template <class T, class M, class U>
    bool operator==(const U &lhs, const zn_t<T, M> &rhs)
    {
        return (lhs - rhs.value()) % rhs.module() == 0 ;
    }

    template <class T, class M>
    bool operator!=(const zn_t<T, M> &lhs, const zn_t<T, M> &rhs)
    {
        return lhs.value() != rhs.value();
    }

    template <class T, class M>
    bool operator!=(const zn_t<T, M> &lhs, const T &rhs)
    {
        return lhs.value() != rhs;
    }

    template <class T, class M>
    bool operator!=(const T &lhs, const zn_t<T, M> &rhs)
    {
        return lhs != rhs.value();
    }
}

#endif