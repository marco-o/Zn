//---------------------------------------------------------------------------------
//
//  Zn 
//  Copyright Marco Oman 2019
//
// Distributed under the Boost Software License, Version 1.0. 
// (See accompanying file LICENSE_1_0.txt or copy at 
// http://www.boost.org/LICENSE_1_0.txt)
//
#ifndef znbasic_H
#define znbasic_H

#include <tuple>
#include <sstream>


char *may_break(void)
{
	return "!";
}


#define ZNASSERT(x) if (!(x)) std::cerr << "Assertion failed: " << #x << may_break() << std::endl ;
#define DBG_SIEVE_ERROR		1
#define DBG_SIEVE_WARNINIG  2
#define DGB_SIEVE_INFO		3
#define DBG_SIEVE_TRACE		4
#define DBG_SIEVE_DEBUG		5
#ifdef _DEBUG
#define DBG_SIEVE			DBG_SIEVE_TRACE
#else
#define DBG_SIEVE			DBG_SIEVE_INFO
#endif

namespace zn
{
	template <class D, class S, bool>
	struct safe_cast_imp {};

	template <class D, class S>
	struct safe_cast_imp<D, S, false>
	{
		static D exec(const S &s)
		{	// For some reason the const_cast is required...
			return const_cast<S &>(s).convert_to<D>();
		}
	};

	template <class D, class S>
	struct safe_cast_imp<D, S, true>
	{
		static D exec(const S &s)
		{
			return static_cast<D>(s);
		}
	};

	template <class D, class S>
	D safe_cast(const S &s)
	{
		return safe_cast_imp<D, S, std::is_arithmetic<S>::value>::exec(s);
	}

	template <class N>
	bool divide_qr1(N &n, const N &d)
	{
		N r, q;
		divide_qr(n, d, q, r);
		if (r == 0)
		{
			n = q;
			return true;
		}
		return false;
	}

	struct system_info_t
	{
		static const size_t memory(void)
		{
#if !defined(_MSC_VER) || defined(_M_X64)
			const size_t max_mem = 4 * 1024 * 1048576LL; // 4GB
#else		// that's just for 32 bits on Windows
			const size_t max_mem = 1 * 512 * 1048576LL; // 0.5GB
#endif	
			return max_mem;
		}
	};


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
	//
	// Some functions define in boost::multipfrecision are adapted
	// here for predefined types
	//
	template <class N>
	typename std::enable_if<std::is_integral<N>::value, int>::type  bit_test(N n, unsigned int bit)
	{
		return (n >> bit) & 1;
	}


	template <class N>
	typename std::enable_if<std::is_integral<N>::value, int>::type  signbit(N n)
	{
		return static_cast<int>(n);
	}

	template <class N>
	typename std::enable_if<std::is_integral<N>::value, unsigned int>::type msb(N n)
	{
		unsigned int i = 0;
		if (signbit(n) < 0)
			n = -n;
		for (n /= 2; n; i++)
			n >>= 1;
		return i;
	}

	template <class N>
	typename std::enable_if<std::is_integral<N>::value>::type divide_qr(N n, N d, N &q, N &r)
	{
		q = n / d;
		r = n % d;
	}


	template <class N, class E>
	N power(N base, const E &exp)
	{
		//   if (signbit(exp) < 0)
		//      return power(base, -exp) ;
		N result(1);
		unsigned int b = msb(exp);
		for (unsigned int i = 0; i <= b; i++)
		{
			if (bit_test(exp, i))
				result *= base;
			base *= base;
		}
		return result;
	}

	template <class N>
	typename std::enable_if<std::is_integral<N>::value, N>::type powm(N b, N e, N p)
	{
		N r = 1;
		for (; e; e >>= 1)
		{
			if (e & 1)
				r = (r * b) % p;
			b = (b * b) % p;
		}
		return r;
	}
}
#endif