//---------------------------------------------------------------------------------
//
//  Zn 
//  Copyright Marco Oman 2019
//
// Distributed under the Boost Software License, Version 1.0. 
// (See accompanying file LICENSE_1_0.txt or copy at 
// http://www.boost.org/LICENSE_1_0.txt)
//
#ifndef znmpqs_H
#define znmpqs_H

#include <tuple>
#include <thread>
#include <vector>
#include <list>
#include <map>
#include "znqueue.h"


namespace zn
{
	template <class large_int, class small_int, class real>
	class quadratic_sieve_base_t
	{
	public:
		struct base_t
		{
			std::vector<small_int>	prime;
			std::vector<small_int>	residue; // quadratic residue
			real					logp;
			base_t(small_int p) : logp(static_cast<real>(-std::log(p)))
			{
			}
			bool valid_for_polynomial(void) const { return prime.size() > 1; }
			small_int prime1(void) const
			{
				return prime[0];
			}
			small_int prime2(void) const
			{
				return prime[1];
			}
			small_int residue2(void) const
			{
				return residue[1];
			}
			base_t operator-(void) const
			{
				base_t result = *this;
				result.residue = prime - residue;
				return result;
			}
			bool eval_residue(const large_int &n)
			{
				large_int n1 = n % prime;
				residue = static_cast<small_int>(quadratic_residue<large_int>(n1, prime, prime / prime0)); // actually a power of prime
				return residue != 0;
			}
			static base_t build(small_int prime, const large_int &n)
			{
				base_t result(prime);
				small_int prime_pwr = prime;
				small_int prime1 = 1;
				small_int prime_power_end = static_cast<small_int>(std::numeric_limits<int>::max() / prime);
				for (int i = 0; (i < 10) && (prime_pwr < prime_power_end); i++)
				{
					small_int n1 = safe_cast<small_int>(n % prime_pwr);
					small_int residue = quadratic_residue<small_int>(n1, prime_pwr, prime1); // actually a power of prime
					if (residue == 0)
						break;
					result.prime.push_back(prime_pwr);
					if (residue > prime_pwr / 2)
						residue = prime_pwr - residue;
					result.residue.push_back(residue);
					// what about employ some kind of lifting?
					prime1 = prime_pwr;
					prime_pwr *= prime;
				}
				return result;
			}
		};
		//
		// this function does an approximate reverse of estimation
		// of prime numbers pi(n) = n / log(n)
		//
		static small_int primes_range(small_int base_size)
		{
			double result = static_cast<double>(base_size);
			for (int i = 0; i < 10; i++)
				result = base_size * std::log(result);
			return static_cast<small_int>(result);
		}

	};


	template <class large_int, class small_int, class real>
	class multiple_polynomial_quadratic_sieve_t : public quadratic_sieve_base_t<large_int, small_int, real>
	{
	public:
		typedef long long long_t;
		typedef unsigned int slot_t;
		enum {bits_per_slot = 8 * sizeof(slot_t) };
		struct base_ref_t : public base_t
		{
			base_ref_t(const base_t &base) : base_t(base) {}

			std::vector<int> smooths; // indexes

		};
		struct polynomial_t
		{
			std::vector<int> index; // base used for the polinomial
			large_int a;
			large_int b;
			large_int c;
			large_int n;
			polynomial_t(const std::vector<int> &idx, 
				         const std::vector<base_ref_t> &base,
						 const large_int &n1) : index(idx), n(n1)
			{
				large_int r = 0;
				const base_ref_t &bp0 = base[idx[0]];
				a = bp0.prime2();
				b = bp0.residue2();
				size_t index_count = index.size();
				for (size_t i = 1 ; i < index_count ; i++)
				{
					const base_ref_t &bp = base[idx[i]];
					large_int g = bp.prime2();
					large_int db = (bp.residue2()- b);
					auto ext = extended_euclidean_algorithm(a, g);
#if DBG_SIEVE >= DBG_SIEVE_TRACE
					if (std::get<0>(ext) != 1)
						throw not_relatively_prime_t<large_int>(db, a);
#endif
					large_int h = (std::get<1>(ext) * db) % g ;
					b = b + h * a;
					a *= g;
					if (b > a / 2)
						b = a - b;
#if DBG_SIEVE >= DBG_SIEVE_TRACE
					if (((b * b) % a) != (n % a))
						throw std::runtime_error("Quadratic residue composition error");
					if (((n - b * b) % a) != 0)
						throw std::runtime_error("Quadratic residue internal error");
#endif
				}
				c = (b * b - n) / a;
			}

			large_int eval(const large_int &x) const
			{
				large_int x1 = a *x + 2 * b;
				return x1 * x + c;
			}
			// polynomial is negative in the range of the roots, including bounds
			std::pair<large_int, large_int> zeros(void) const
			{
				large_int d = safe_cast<large_int>(sqrt(n));
				large_int x1 = (-b - d) / a ;
				large_int x2 = (-b + d) / a;
#if DBG_SIEVE >= DBG_SIEVE_TRACE
				large_int y1 = eval(x1 - 1);
				large_int y0 = eval(x1 );
				if (y1 < 0 || y0 > 0)
					throw std::runtime_error("Error in finding x1");
				y1 = eval(x2 + 1);
				y0 = eval(x2);
				if (y1 < 0 || y0 > 0)
					throw std::runtime_error("Error in finding x2");
#endif
				return std::make_pair(x1, x2);
			}
		};
		enum smooth_status_e{ smooth_idle_e, smooth_valid_e, smooth_candidate_e };

		multiple_polynomial_quadratic_sieve_t(const large_int &n, const large_int &m, small_int base_size) : n_(n), m_(m)
		{
			// double the range; half of them won't be a quadratic residue
			small_int range;
			if (base_size != 0)
				range = primes_range(base_size * 2);
			else
			{
				double n1 = safe_cast<double>(n);
				double e1 = std::log(n1);
				double e2 = std::sqrt(e1 * std::log(e1)) / 2;
				range = static_cast<small_int>(std::exp(e2));
			}
			auto primes = eratosthenes_sieve<small_int>(static_cast<int>(range));
			small_int r = 1;
			for (auto p : primes)
			{
				base_t base = base_t::build(p, n);
				if (base.prime.size() > 0)
					base_.push_back(base);
			}
			small_int largest_sieving_prime = base_.rbegin()->prime[0];
#if DBG_SIEVE >= DBG_SIEVE_INFO
			std::cout << "Actual base size: " << base_.size() << ", largest = " << largest_sieving_prime << std::endl;
#endif // DBG_SIEVE	
			sieve_thrs_ = safe_cast<real>(std::log(largest_sieving_prime) * 2);
			smooth_thrs_ = largest_sieving_prime;
			smooth_thrs_ *= smooth_thrs_;
			if (m_ == 0)
				m_ = large_int(largest_sieving_prime) * 5;
		}
		large_int process(void)
		{
			large_int n1 = 2 * n_;
			large_int a2 = safe_cast<large_int>(sqrt(n1)) / m_;
#if DBG_SIEVE >= DBG_SIEVE_INFO
			std::cout << "a2 = " << a2 << std::endl;
#endif // DBG_SIEVE	
			std::vector<int> idx;
			idx.push_back(1);
			idx.push_back(2);
			polynomial_t p(idx, base_, n_);
			p.zeros();
			sieve(p);
			return 1;
		}
	private:
		void sieve(const polynomial_t &poly)
		{
			// build vector for sieving
			std::vector<real> values(safe_cast<size_t>(2 * m_));
			auto zeros = poly.zeros();
			fill_range(poly, values, -m_, zeros.first);
			fill_range(poly, values, zeros.first, zeros.second);
			fill_range(poly, values, zeros.second, m_);

			// use the base for sieving

			// collect smooth numbers
		}
		void fill_range(const polynomial_t &poly, 
			            std::vector<real> &values, 
			            const large_int &begin, 
			            const large_int &end)
		{
			large_int mid = (begin + end) / 2;
			large_int y_1 = poly.eval(begin);
			large_int y0  = poly.eval(mid);
			large_int y1  = poly.eval(end - 1);
			real t_1 = std::abs(safe_cast<real>(y_1));
			real t0 = std::abs(safe_cast<real>(y0));
			real t1 = std::abs(safe_cast<real>(y1));
			real q1 = t_1 / t0 ;
			real q0 = t0 / t1 ;
			real t = q1 / q0 + q0 / q1 - 2;
			if (t < 0.1)
			{
				fill_linear(values, begin, y_1, mid, y0);
				fill_linear(values, mid, y0, end, y1);
			}
			else 
				if (end - begin < 32)
					fill_exact(poly, values, begin, end);
				else
				{
					fill_range(poly, values, begin, mid);
					fill_range(poly, values, mid, end);
				}
		}
		void fill_linear(std::vector<real> &values, 
			             const large_int &x1, const large_int &y1, 
			             const large_int &x2, const large_int &y2)
		{
			size_t size = safe_cast<size_t>(x2 - x1);
			size_t offset = safe_cast<size_t>(x1 + m_);
			real t1 = std::log(std::abs(safe_cast<real>(y1)));
			real t2 = std::log(std::abs(safe_cast<real>(y2)));
			real m = (t2 - t1) / size;
			for (size_t i = 0; i < size; i++)
				values[i + offset] = t1 + m * i;
		}
		void fill_exact(const polynomial_t &poly, std::vector<real> &values, large_int begin, large_int end)
		{
			size_t size = safe_cast<size_t>(end - begin);
			size_t offset = safe_cast<size_t>(begin + m_);
			for (size_t i = 0; i < size; i++)
			{
				auto y = poly.eval(begin + i);
				values[offset + i] = std::log(std::abs(safe_cast<real>(y)));
			}
		}
		real						 sieve_thrs_;
		large_int					 smooth_thrs_; // square of last element of base
		large_int					 n_; // number to factor
		large_int					 m_; // size of sieving interval
		std::vector<base_ref_t>		 base_;
	};



	template <class large_int, class small_int = int, class real = float>
	large_int multiple_polynomial_quadratic_sieve(const large_int &n, const large_int &m, small_int base_size)
	{
#if DBG_SIEVE >= DBG_SIEVE_INFO
		std::cout << "Factorization of " << n << std::endl;
#endif
		multiple_polynomial_quadratic_sieve_t<large_int, small_int, real> qs(n, m, base_size);
		return qs.process();
	}

}

#endif
