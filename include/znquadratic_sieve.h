//---------------------------------------------------------------------------------
//
//  Zn 
//  Copyright Marco Oman 2019
//
// Distributed under the Boost Software License, Version 1.0. 
// (See accompanying file LICENSE_1_0.txt or copy at 
// http://www.boost.org/LICENSE_1_0.txt)
//
#ifndef znquadratic_sieve_H
#define znquadratic_sieve_H

#include <tuple>

char *may_break(void)
{
	return "!";
}

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
#define ZNASSERT(x) if (!(x)) std::cerr << "Assertion failed: " << #x << may_break() << std::endl ;
#define DBG_SIEVE

namespace zn
{
	template <class large_int, class small_int, class real>
	class quadratic_sieve_t
	{
	public:
		typedef typename std::make_unsigned<small_int>::type small_uint;
		enum { small_uint_size = sizeof(small_uint) * 8 };
		struct base_t
		{
			small_int	prime;
			small_int	residue; // quadratic residue
			real		logp;
#ifdef DBG_SIEVE
			small_int   prime0;
#endif // DBG_SIEVE

			base_t(small_int p, small_int r) : 
				prime(p), residue(r), logp(static_cast<real>(-std::log(p)))
			{
#ifdef DBG_SIEVE
				prime0 = prime;
#endif // DBG_SIEVE			
			}
			base_t operator-(void) const
			{
				base_t result = *this;
				result.residue = prime - residue;
				return result;
			}
			void compose(const base_t &rhs, const large_int &n)
			{
				prime *= rhs.prime;
				logp += rhs.logp;
				small_int n1 = safe_cast<small_int>(n % prime);
				residue = quadratic_residue(n1, prime); // actually a power of prime
			}
		};
		struct smooth_t
		{
			large_int				n ;
			std::vector<small_uint>	e ;
			smooth_t(const large_int &n1, int base_size) : 
				n(n1), e((base_size + small_uint_size - 1) / small_uint_size, 0) {}
			void toggle(int pos)
			{
				int slot = pos / small_uint_size;
				int bit  = pos % small_uint_size;
				e[slot] ^= 1 << bit;
			}
		};
		quadratic_sieve_t(const large_int &n, small_int base_size) : n_(n)
		{
			// double the range; half of them won't be a quadratic residue
			small_int range = primes_range(base_size * 2);
			auto primes = eratosthenes_sieve<small_int>(range);
			small_int r;
			for (auto p : primes)
			{
				small_int n1 = safe_cast<small_int>(n % p);
				if ((r = quadratic_residue(n1, p)) != 0)
					base_.push_back(base_t(p, r));
			}
		}
		void sieve(void)
		{
			large_int n1 = safe_cast<large_int>(sqrt(n_) + 1);
			small_int range = static_cast<small_int>(std::pow(base_.size(), 2.8));
			auto values = build_sieving_range(n1, range);
			sieve_range(values, n1);
			small_int base_size = base_.size();
			large_int q, r;
			std::vector<smooth_t> smooths;
			for (small_int i = 0; i < range; i++)
				if (values[i] < 1)
				{
					large_int n = n1 + i;
					n = n * n - n_;
					smooth_t s(n, base_.size());
					for (small_int j = 0 ; j < base_size ; j++)
					{
						const auto &base = base_[j];
						for (; ;)
						{
							divide_qr<large_int>(n, base.prime, q, r);
							if (signbit(r) == 0)
							{
								s.toggle(j);
								n = q;
							}
							else
								break;
						}
					}
					if (n == 1)
						smooths.push_back(s);
					else
						std::cout << (n1 + i) << ", " << values[i] << std::endl;
				}
			std::sort(values.begin(), values.end());
		}
	private:
		void sieve_range(std::vector<real> &values, const large_int &begin)
		{
			auto size = values.size();
			for (const auto &base : base_)
			{
				sieve_range(values, begin, base);
				if (base.prime != 2)
					sieve_range(values, begin, -base);
				base_t powers = base;
				for (int i = 0; i < 10 && (powers.prime < static_cast<small_int>(size / base.prime)); i++)
				{
					powers.compose(base, n_);
					if (powers.residue == 0)
						continue;
					sieve_range(values, begin, powers);
					if (powers.prime != 2)
						sieve_range(values, begin, -powers);
				}
			}
		}
		void sieve_range(std::vector<real> &values, const large_int &begin, const base_t &base)
		{
			large_int n = (begin / base.prime) * base.prime + base.residue;
			if (n < begin)
				n += base.prime;
			auto size = values.size();
			std::vector<real>::size_type pos = safe_cast<small_int, large_int>(n - begin);
			ZNASSERT(((n *n - n_) % base.prime) == 0);
#ifdef DBG_SIEVE
			int errs = 0;
#endif
			for (; pos < size; pos += base.prime)
			{
				values[pos] += base.logp;
#ifdef DBG_SIEVE
				if (ns_[pos] % base.prime0 != 0)
					errs++;
				else
					ns_[pos] /= base.prime0;
#endif
			}

		}
		std::vector<real> build_sieving_range(large_int n1, small_int range)
		{
			std::vector<real> data(range);
			large_int n2 = n1 * n1 - n_;
			for (small_int i = 0; i < range; i++)
			{
				data[i] = std::log(safe_cast<real>(n2));
#ifdef DBG_SIEVE
				ns_.push_back(n2);
#endif
				n2 += n1 * 2 + 1;
				n1++;
			}
			return data;
		}
		//
		// this function does an approximate reverse of estimation
		// of prime numbers pi(n) = n / log(n)
		//
		small_int primes_range(small_int base_size)
		{
			double result = base_size;
			for (int i = 0; i < 10; i++)
				result = base_size * std::log(result);
			return static_cast<int>(result);
		}

		large_int			n_;
		std::vector<base_t> base_;
#ifdef DBG_SIEVE
		std::vector<large_int> ns_;
#endif
	};



	template <class large_int, class small_int = int, class real = float>
	large_int quadratic_sieve(const large_int &n, small_int base_size)
	{
		std::cout << "Factorization of " << n << std::endl;
		quadratic_sieve_t<large_int, small_int, real> qs(n, base_size);
		qs.sieve();
		return n;
	}

}

#endif
