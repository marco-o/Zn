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
#include <thread>
#include <vector>
#include <list>
#include <map>
#include "znqueue.h"

namespace zn
{
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
#ifdef _M_X64
			const size_t max_mem = 4 * 1024 * 1048576LL; // 4GB
#else
			const size_t max_mem = 1 * 512 * 1048576LL; // 0.5GB
#endif	
			return max_mem;
		}
	};

	template <class large_int, class small_int, class real>
	class quadratic_sieve_t
	{
	public:
		typedef long long long_t;
		typedef unsigned int slot_t;
		enum {bits_per_slot = 8 * sizeof(slot_t) };
		struct sieve_range_t
		{
			large_int first; // beginning
			small_int second; // size
			// additional info on polynomial, ecc
			int		dir_sign; // negative grows backwards
			large_int final_remainder; // value of polinomial at end of range
			sieve_range_t(const large_int &f = large_int(), const small_int &s = small_int()) : first(f), second(s) {}
			sieve_range_t &operator++(void)
			{
				first += dir_sign * second;
				return *this;
			}
		};
		struct base_t
		{
			small_int	prime;
			small_int	residue; // quadratic residue
			real		logp;
			small_int   prime0;
			base_t(small_int p, small_int r) :
				prime(p), residue(r), logp(static_cast<real>(-std::log(p)))
			{
				prime0 = prime;
			}
			base_t operator-(void) const
			{
				base_t result = *this;
				result.residue = prime - residue;
				return result;
			}
			// what about substitute composition with lifting?
			void compose(const base_t &rhs, const large_int &n)
			{
				prime *= rhs.prime;
				large_int n1 = n % prime;
				residue = static_cast<small_int>(quadratic_residue<large_int>(n1, prime, prime / rhs.prime0)); // actually a power of prime
			}
		};
		struct base_ref_t : public base_t
		{
			base_ref_t(small_int p, small_int r) : base_t(p, r) {}

			std::vector<int> smooths; // indexes

		};
		enum smooth_status_e{ smooth_idle_e, smooth_valid_e, smooth_candidate_e };
		struct smooth_t
		{
			large_int				n; // number squared
			std::vector<int>		factors;
			large_int				f; // remainder after trial division
			large_int				sqr;
			bool					sign_bit; // true if negative
			smooth_status_e			s;
			smooth_t(void) : n(1), f(1), sqr(1), sign_bit(false), s(smooth_valid_e) {}

			smooth_t(const large_int &n1, 
				     const large_int &r1, 
				     const large_int &m, 
				     const large_int &thrs, 
				     const std::vector<base_ref_t> &base) : n(n1), f(r1), sqr(1), sign_bit(false) 
			{
				if (f < 0)
				{
					sign_bit = true;
					f = -f;
				}
				size_t base_size = base.size();
				for (size_t j = 0; j < base_size; j++)
				{
					const auto &b = base[j];
					int rexp = 0;
					large_int p = b.prime;
					int power = 0;
					while (divide_qr1(f, p))
						if (++power % 2 == 0)
							sqr = (sqr * b.prime) % m;
					if (power & 1)
						factors.push_back(static_cast<int>(j));
				}
				if (f == 1)
					s = smooth_valid_e;
				else if (f < thrs)
					s = smooth_candidate_e;
				else
					s = smooth_idle_e;
			}
			bool valid(void) const { return s != smooth_idle_e; }
			large_int result(const large_int &m)
			{
				large_int a = gcd(m, n + sqr);
				if (a != 1 && a != m)
					return a;
				large_int b = (n > sqr ? n - sqr: sqr - n);
				return gcd(b, m);
			}
			void compose(const smooth_t &rhs, const large_int &m, const std::vector<base_ref_t> &base)
			{
				n = (n * rhs.n) % m;
				sqr = (sqr * rhs.sqr) % m;
				if (rhs.f == f)
				{
					sqr = (sqr * f) % m;
					f = 1;
				}

				sign_bit ^= rhs.sign_bit;
				std::vector<int>		fact;
				auto it1  = factors.begin();
				auto end1 = factors.end();
				auto it2 = rhs.factors.begin();
				auto end2 = rhs.factors.end();
				while (it1 != end1 && it2 != end2)
				{
					if (*it1 < *it2)
						fact.push_back(*it1++);
					else if (*it1 > *it2)
						fact.push_back(*it2++);
					else // same factor; gets skipped!
					{
						sqr = (sqr * base[*it1].prime) % m;
						++it1;
						++it2;
					}
				}
				if (it1 != end1)
					fact.insert(fact.end(), it1, end1);
				if (it2 != end2)
					fact.insert(fact.end(), it2, end2);
				factors = fact;
				s = smooth_valid_e;
			}
			bool invariant(const large_int &m, const std::vector<base_ref_t> &base) const
			{
				large_int s1 = (sqr * sqr * f) % m ;
				for (auto idx : factors)
					s1 = (s1 * base[idx].prime) % m;
				s1 = (s1 - n * n) % m;
				return s1 == 0;
			}
		};
		typedef std::vector<smooth_t> smooth_vect_t;
		typedef std::map<large_int, smooth_t> candidates_map_t;

		quadratic_sieve_t(const large_int &n, small_int base_size) : n_(n)
		{
			// double the range; half of them won't be a quadratic residue
			small_int range = primes_range(base_size * 2);
			auto primes = eratosthenes_sieve<small_int>(static_cast<int>(range));
			small_int r;
			for (auto p : primes)
			{
				small_int n1 = safe_cast<small_int>(n % p);
				if ((r = quadratic_residue(n1, p)) != 0)
					base_.push_back(base_ref_t(p, r));
			}
#if DBG_SIEVE >= DBG_SIEVE_INFO
			std::cout << "Actual base size: " << base_.size() << ", largest = " << base_.rbegin()->prime << std::endl;
#endif // DBG_SIEVE	
			sieve_thrs_ = safe_cast<real>(std::log(*primes.rbegin()) * 1.8);
			smooth_thrs_ = base_.rbegin()->prime;
			smooth_thrs_ *= smooth_thrs_;
		}
		large_int sieve(void)
		{
			sieve_range_t range;
			range.first = safe_cast<large_int>(sqrt(n_) + 1);
			range.second = static_cast<small_int>(std::pow(base_.size(), 2.5));
			const size_t max_mem = system_info_t::memory();
			auto cores = 1; // std::thread::hardware_concurrency();
			range.second = std::min<small_int>(range.second, static_cast<small_int>(max_mem / (sizeof(real) * cores)));

			smooth_vect_t smooths;
			candidates_map_t candidates;
			for ( ; smooths.size() < base_.size() ; ++range)
			{
#if DBG_SIEVE >= DBG_SIEVE_INFO
				std::cout << "Sieving " << range.first 
					      << " , " << range.second << "(" << smooths.size() << " )\r" << std::flush;
#endif
				auto smooth = sieve_range(range);
				for (auto &s : smooth)
				{
					if (s.s == smooth_candidate_e)
					{
						auto it = candidates.find(s.f);
						if (it != candidates.end())
						{
							s.compose(it->second, n_, base_);
#if DBG_SIEVE >= DBG_SIEVE_DEBUG
							if (!s.invariant(n_, base_))
								std::cout << "Hmmm";
#endif // DBG_SIEVE	
						}
						else // just put it aside
						{
							candidates[s.f] = s;
							continue;
						}
					}
					int si = static_cast<int>(smooths.size());
					for (auto f : s.factors)
						base_[f].smooths.push_back(si);
					smooths.push_back(s);
				}
			}
			auto result = solve(smooths);
			for (auto &item : result)
			{
				smooth_t s;
				for (auto index : item)
					s.compose(smooths[index], n_, base_);
				large_int r = s.result(n_);
				if (r != 1 && r != n_)
					return r;
			}
			return 1;
		}
	private:
		// linear system: rows is #base_, cols is number of smooths
		std::vector<std::vector<int>> solve(std::vector<smooth_t> &smooths)
		{
			size_t smooth_size = smooths.size();
			size_t base_size = base_.size() ;
			size_t slots = (smooth_size + bits_per_slot - 1) / bits_per_slot;
			std::vector<std::vector<slot_t>> matrix(base_size, std::vector<slot_t>(slots, 0));
			for (size_t i = 0 ; i < smooth_size; i++)
			{
				const smooth_t &smooth = smooths[i];
				slot_t slot = static_cast<slot_t>(i / bits_per_slot);
				slot_t mask = 1 << static_cast<slot_t>(i % bits_per_slot);
				if (smooth.sign_bit)
					matrix[0][slot] |= mask;
				for (auto idx : smooth.factors)
					matrix[idx + 1][slot] |= mask;
			}
			std::vector<int> base_perm(base_size);
			for (size_t i = 0; i < base_size; i++)
				base_perm[i] = i;
			for (size_t i = 0; i < base_size; i++)
			{
				size_t j;
				slot_t slot = i / bits_per_slot;
				slot_t mask = 1 << static_cast<slot_t>(i % bits_per_slot);
				for (j = i ; j < base_size ; j++) // pivoting
					if (matrix[i][slot] & mask)
					{
						if (i != j)
							std::swap(base_perm[i], base_perm[j]);
						break;
					}
				if (j == base_size)
					continue; // nothing to do...
				auto &v = matrix[i];
				for (j = i + 1; j < base_size; j++)
					for (size_t k = slot; k < slots; k++)
						matrix[j][k] ^= v[k];

			}
			std::vector<std::vector<int>> result;
			for (size_t i = base_size; i < smooth_size; i++)
			{
				std::vector<int> idx;
				slot_t slot = i / bits_per_slot;
				slot_t mask = 1 << static_cast<slot_t>(i % bits_per_slot);

				for (size_t j = 0; j < base_size; j++)
					if (matrix[j][slot] & mask)
						idx.push_back(base_perm[j]);
				idx.push_back(i);
				result.push_back(idx);
				if (result.size() > 30)
					break;
			}
			return result;
		}
		bool bit_test(const std::vector<slot_t> &v, size_t index)
		{
			size_t slot = i / bits_per_slot;
			return (v[slot] & (1 << static_cast<slot_t>(i % bits_per_slot))) != 0;
		}
		std::vector<smooth_t> collect_smooth(const sieve_range_t &range,
											 const std::vector<real> &values)
		{
			std::vector<smooth_t> smooths;
			for (int i = 0; i < range.second; i++)
				if (values[i] < sieve_thrs_)
				{
					large_int n = range.first + i;
					large_int r = n * n - n_;
					smooth_t s(n, r, n_, smooth_thrs_, base_);
					if (s.valid())
					{
						smooths.push_back(s);
#if DBG_SIEVE >= DBG_SIEVE_INFO
						if (!s.invariant(n_, base_))
							std::cout << "Hmm";
#endif
					}
#if DBG_SIEVE >= DBG_SIEVE_DEBUG
					else
						std::cout << (range.first + i) << ", " << values[i] << std::endl;
#endif
				}
			return smooths;
		}
		smooth_vect_t sieve_range(const sieve_range_t &range)
		{
			std::vector<real> values;
			values.reserve(static_cast<size_t>(range.second));
			build_sieving_range(range, values);
			sieve_range(values, range.first);
			return collect_smooth(range, values);
		}
		void sieve_range(std::vector<real> &values, const large_int &begin)
		{
			auto size = values.size();
			for (const auto &base : base_)
			{
				sieve_range(values, begin, base);
				if (base.prime != 2)
					sieve_range(values, begin, -base);
				base_t powers = base;
				small_int prime_power_end = static_cast<small_int>(std::numeric_limits<small_int>::max() / base.prime);
				for (int i = 0; (i < 10) && (powers.prime < prime_power_end); i++)
				{
					powers.compose(base, n_);
					if (powers.residue == 0)
						break;
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
			size_t pos = safe_cast<size_t, large_int>(n - begin);
#if DBG_SIEVE >= DBG_SIEVE_DEBUG
			if (((n *n - n_) % base.prime) != 0)
			{
				std::cout << "Error: n = " << n 
					      << "\nn_ = " << n_ 
					      << "\np = " << base.prime 
					      << "\nr = " << base.residue 
					      << "\np0= " << base.prime0
					      << std::endl;
				throw std::runtime_error("Quadratic residue problem");
			}
			int errs = 0;
#endif
			for (; pos < size; pos += static_cast<size_t>(base.prime))
			{
#if DBG_SIEVE >= DBG_SIEVE_DEBUG
				if (abs(values[pos] - std::log(safe_cast<real>(ns_[pos]))) > 1e-4)
					std::cout << "Hey!\n";
#endif
				values[pos] += base.logp;
#if DBG_SIEVE >= DBG_SIEVE_DEBUG
				if (ns_[pos] % base.prime0 != 0)
					errs++;
				else
					ns_[pos] /= base.prime0;
				if (abs(values[pos] - std::log(safe_cast<real>(ns_[pos]))) > 1e-4)
					std::cout << "Hey!\n";
#endif
			}
		}
		void build_sieving_range_exact(const sieve_range_t &range, std::vector<real> &values)
		{
			large_int n1 = range.first;
			large_int n2 = n1 * n1 - n_;
			std::vector<real> data(static_cast<size_t>(range.second));
			for (small_int i = 0; i < range.second; i++)
			{
				values.push_back(std::log(safe_cast<real>(n2)));
#if DBG_SIEVE >= DBG_SIEVE_DEBUG
				ns_.push_back(n2);
#endif
				n2 += n1 * 2 + 1;
				n1++;
			}
		}
		void  build_sieving_range(const sieve_range_t &range, std::vector<real> &values)
		{
			large_int n1 = range.first;
			large_int n2 = n1 * n1 - n_;
			large_int m1 = range.first + range.second;
			large_int m2 = m1 * m1 - n_;
			real rn = std::abs(safe_cast<real>(n2));
			real rm = std::abs(safe_cast<real>(m2));
			real t = rn / rm + rm / rn - 2;
			if (t < 3e-3)
			{
				std::vector<real> data(static_cast<size_t>(range.second));
				rn = std::log(rn);
				rm = std::log(rm);
				real delta = (rm - rn) / range.second;
				for (small_int i = 0; i < range.second; i++)
					values.push_back(rn + i * delta);
			}
			else if (range.second < 32)
				build_sieving_range_exact(range, values);
			else
			{
				sieve_range_t r1(range.first, range.second / 2);
				build_sieving_range(r1, values);
				sieve_range_t r2(r1.first + r1.second, range.second - r1.second);
				build_sieving_range(r2, values);
			}
		}
		// removes element that appears only once
		void sieving_thread(void)
		{
			try
			{
				auto range = ranges_to_sieve_.pop();
				for (; range.second != 0; range = ranges_to_sieve_.pop())

					smooths_found_.push(sieve_range(range));
			}
			catch (std::exception &exc)
			{
				std::cout << "Exception in thread: " << exc.what() << std::endl;
			}
		}
		//
		// this function does an approximate reverse of estimation
		// of prime numbers pi(n) = n / log(n)
		//
		small_int primes_range(small_int base_size)
		{
			double result = static_cast<double>(base_size);
			for (int i = 0; i < 10; i++)
				result = base_size * std::log(result);
			return static_cast<small_int>(result);
		}
		real				sieve_thrs_;
		large_int			smooth_thrs_; // square of last element of base
		large_int			n_;
		std::vector<base_ref_t>		 base_;
		shared_list_t<sieve_range_t> ranges_to_sieve_;
		shared_list_t<smooth_vect_t> smooths_found_;
		std::vector<std::thread>	 threads_;

	};



	template <class large_int, class small_int = int, class real = float>
	large_int quadratic_sieve(const large_int &n, small_int base_size)
	{
#if DBG_SIEVE >= DBG_SIEVE_INFO
		std::cout << "Factorization of " << n << std::endl;
#endif
		quadratic_sieve_t<large_int, small_int, real> qs(n, base_size);
		return qs.sieve();
	}

}

#endif
