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
	//
	// Some stuff that may be useful 
	//
	inline bool bit_test(int n, int bit)
	{
		return (n >> bit) & 1;
	}

	inline int msb(unsigned int n)
	{
		int i = 0;
		for (n /= 2; n; i++)
			n >>= 1;
		return i;
	}

	inline void divide_qr(int n, int d, int &q, int &r)
	{
		q = n / d;
		r = n % d;
	}

	inline int signbit(int n)
	{
		return n;
	}

	template <class N, class E>
	N power(N base, const E &exp)
	{
		//   if (signbit(exp) < 0)
		//      return power(base, -exp) ;
		N result(1);
		size_t b = msb(exp);
		for (size_t i = 0; i <= b; i++)
		{
			if (bit_test(exp, i))
				result *= base;
			base *= base;
		}
		return result;
	}

	int powm(int b, int e, int p)
	{
		int r = 1;
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