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
#define DBG_SIEVE_ERROR		1
#define DBG_SIEVE_WARNINIG  2
#define DGB_SIEVE_INFO		3
#define DBG_SIEVE_TRACE		4
#define DBG_SIEVE_DEBUG		5
#define DBG_SIEVE			DBG_SIEVE_TRACE


namespace zn
{
	template <class large_int, class small_int, class real>
	class quadratic_sieve_t
	{
	public:
		typedef typename std::make_unsigned<small_int>::type small_uint;
		enum { small_uint_bits = sizeof(small_uint) * 8, max_base_count = 4 };
		struct base_t
		{
			small_int	prime;
			small_int	residue; // quadratic residue
			int			count ;
			int			index[max_base_count];
			real		logp;
#if DBG_SIEVE >= DBG_SIEVE_DEBUG
			small_int   prime0;
#endif // DBG_SIEVE

			base_t(small_int p, small_int r) : 
				prime(p), residue(r), logp(static_cast<real>(-std::log(p)))
			{
#if DBG_SIEVE >= DBG_SIEVE_DEBUG
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
				large_int n1 = n % prime;
				residue = static_cast<int>(quadratic_residue<large_int>(n1, prime, prime / rhs.prime)); // actually a power of prime
			}
		};
		enum smooth_status_e { smooth_delete_e, smooth_touch_e, smooth_valid_e };
		struct smooth_t
		{
			smooth_status_e			s;
			large_int				n;
			large_int				r;
			std::vector<small_uint>	e ;
			smooth_t(large_int r1, const large_int &n1, const large_int &m, const std::vector<base_t> &base) :
				r(1), n(n1), e(base.size() / small_uint_bits + 1, 0), s(smooth_valid_e) 
			{
				large_int r2, q;
				small_int base_size = base.size();
				for (small_int j = 0; j < base_size; j++)
				{
					const auto &b = base[j];
					int rexp = 0;
					large_int p = b.prime;
					for (; ;)
					{
						divide_qr(r1, p, q, r2);
						if (r2 == 0)
						{
							r1 = q;
							toggle(j);
							if (++rexp == 2)
							{
								r = (r * b.prime) % m;
								rexp = 0;
							}
						}
						else
							break;
					}
				}
				if (r1 != 1)
					s = smooth_delete_e;
			}
			bool valid(void) const { return s != smooth_delete_e; }
			void toggle(int pos)
			{
				int slot = pos / small_uint_bits;
				int bit  = pos % small_uint_bits;
				e[slot] ^= 1 << bit;
			}
			void combine(smooth_t &rhs, const large_int &m, const std::vector<base_t> &base)
			{
				s = smooth_touch_e;
				n = (n * rhs.n) % m;
				r = (r * rhs.r) % m;
				size_t size = e.size();
				for (size_t i = 0; i < size; i++)
				{
					small_uint v = e[i] & rhs.e[i];
					e[i] ^= rhs.e[i];
					if (v)
						for (int b = 0; b < small_uint_bits; b++)
							if (bit_test(v, b))
								r = (r * base[i * small_uint_bits + b].prime) % m;
				}
				rhs.s = smooth_delete_e;
			}
			// n * n = r * r * odd_factors_of_base
			bool verify(const large_int &m, const std::vector<base_t> &base)
			{
				small_uint v;
				large_int r1 = r * r % m; 
				size_t size = e.size();
				for (size_t i = 0; i < size; i++)
				{
					if ((v = e[i]) != 0)
						for (int b = 0; b < small_uint_bits; b++)
							if (bit_test(v, b))
								r1 = (r1 * base[i * small_uint_bits + b].prime) % m;
				}
				large_int n1 = (n * n - r1) % m;
				return n1 == 0;
			}
			void remap(const std::vector<int> &index, size_t base_slots)
			{
				size_t size = e.size();
				for (size_t i = 0; i < size; i++)
					if (e[i])
					{
						auto value = e[i];
						e[i] = 0; // a bit tricky but ok
						for (size_t j = 0; j < small_uint_bits; j++)
							if (bit_test(value, j))
							{
								size_t old_pos = i * small_uint_bits + j;
								size_t new_pos = index[old_pos];
								int slot = new_pos / small_uint_bits;
								int bit = new_pos % small_uint_bits;
								e[slot] ^= 1 << bit;
							}
					}
				e.erase(e.begin() + base_slots, e.end());
				s = smooth_valid_e;
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
		large_int sieve(void)
		{
			std::pair<large_int, small_int> range;
			range.first = safe_cast<large_int>(sqrt(n_) + 1);
			range.second = static_cast<small_int>(std::pow(base_.size(), 2.5));
			auto values = build_sieving_range(range);
			sieve_range(values, range.first);
#if DBG_SIEVE >= DBG_SIEVE_TRACE
			std::cout << "Sieved  " << values.size() << " values" << std::endl;
#endif
			std::vector<smooth_t> smooths;
			collect_smooth(smooths, range, values);
#if DBG_SIEVE >= DBG_SIEVE_INFO
			std::cout << "Found  " << smooths.size() << " values" << std::endl;
#endif
#if DBG_SIEVE >= DBG_SIEVE_DEBUG
			std::sort(values.begin(), values.end());
#endif
			for (; ;)
			{
				int eraseable = compact_base(smooths);
				for (int count = 2; count <= max_base_count; count++)
					eraseable += simplyfy_cycles(smooths, count);
#if DBG_SIEVE >= DBG_SIEVE_INFO
				std::cout << "Reduced base size = " << base_.size() - eraseable << std::endl;
#endif
				if (eraseable > 0)
					compact(smooths);
				else
					break;
			}
			return solve(smooths);
		}
	private:
		large_int solve(std::vector<smooth_t> &smooths)
		{
			size_t base_size = base_.size() ;
			size_t smooth_size = smooths.size();
			size_t key_col = base_size + 1;
			int slot = key_col / small_uint_bits;
			int bit = key_col % small_uint_bits;
	/*		small_uint bitp = 1 << bit;
			for (auto &smooth : smooths)
				smooth.e[slot] |= bitp;*/
			// std gaussian elimination
			for (size_t i = 0; i < base_size; i++)
			{
				// search pivot
				slot = i / small_uint_bits;
				bit  = i % small_uint_bits;
				for (size_t j = i ; j < smooth_size ; j++)
					if (bit_test(smooths[j].e[slot], bit))
					{
						std::swap(smooths[j], smooths[i]);
						break;
					}
				if (!bit_test<small_uint>(smooths[i].e[slot], bit)) // not found
				{
					base_[i].count = 0;
					continue;
				}
				smooth_t &smooth = smooths[i];
				for (size_t j = i + 1; j < smooth_size; j++)
					if (bit_test<small_uint>(smooths[j].e[slot], bit))
					{
						smooths[j].combine(smooth, n_, base_);
#if DBG_SIEVE >= DBG_SIEVE_DEBUG
						smooths[j].verify(n_, base_);
#endif
					}
			}
			for (size_t i = base_size; i < smooth_size; i++)
			{
				large_int a = smooths[i].n;
				large_int b = smooths[i].r;
#if DBG_SIEVE >= DBG_SIEVE_TRACE
				if (a < b)
					std::swap(a, b);
				std::cout << "Testing a = " << a << ", b = " << b 
					      << " ( " << ((a * a - b * b) % n_) << ")" << std::endl;
#endif
				large_int p1 = gcd(a + b, n_);
				if (p1 != 1 && p1 != n_)
					return p1;
				if (a < b)
					p1 = gcd(b - a, n_);
				else
					p1 = gcd(a - b, n_);
				if (p1 != 1 && p1 != n_)
					return p1;
			}
			return 1;
		}
		void collect_smooth(std::vector<smooth_t> &smooths, 
			                const std::pair<large_int, small_int> &range, 
			                const std::vector<real> &values)
		{
			small_int base_size = base_.size();
			for (small_int i = 0; i < range.second; i++)
				if (values[i] < 3)
				{
					large_int n = range.first + i;
					large_int r = n * n - n_;
					smooth_t s(r, n, n_, base_);
					if (s.valid())
					{
						smooths.push_back(s);
#if DBG_SIEVE >= DBG_SIEVE_DEBUG
						s.verify(n_, base_);
#endif
					}
#if DBG_SIEVE >= DBG_SIEVE_DEBUG
					else
						std::cout << (range.first + i) << ", " << values[i] << std::endl;
#endif
				}
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
			std::vector<real>::size_type pos = safe_cast<small_int, large_int>(n - begin);
			ZNASSERT(((n *n - n_) % base.prime) == 0);
#if DBG_SIEVE >= DBG_SIEVE_DEBUG
			int errs = 0;
#endif
			for (; pos < size; pos += base.prime)
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
		std::vector<real> build_sieving_range(const std::pair<large_int, small_int> &range)
		{
			large_int n1 = range.first;
			large_int n2 = n1 * n1 - n_;
			std::vector<real> data(range.second);
			for (small_int i = 0; i < range.second; i++)
			{
				data[i] = std::log(safe_cast<real>(n2));
#if DBG_SIEVE >= DBG_SIEVE_DEBUG
				ns_.push_back(n2);
#endif
				n2 += n1 * 2 + 1;
				n1++;
			}
			return data;
		}
		// removes element that appears only once
		size_t compact_base(std::vector<smooth_t> &smooths)
		{
			for (auto &base : base_)
				base.count = 0;

			size_t smooth_count = smooths.size();
			for (size_t k = 0 ; k < smooth_count ; k++)
			{
				const auto &smooth = smooths[k];
				size_t slots = smooth.e.size();
				for (size_t i = 0 ; i < slots ; i++)
					if (smooth.e[i])
					{
						auto slot = smooth.e[i];
						for (size_t j = 0; j < small_uint_bits; j++)
							if (bit_test(slot, j))
							{
								auto &base = base_[i * small_uint_bits + j];
								if (base.count < max_base_count)
									base.index[base.count] = k;
								base.count++;
							}
					}
			}
			// mark smooth numbers to remove
			size_t result = 0;
			for (auto &base : base_)
				if (base.count < 2)
				{
					if (base.count == 1)
						smooths[base.index[0]].s = smooth_delete_e;
					base.count = 0; // mark removed
					result++;
				}
			return result;
		}
		size_t simplyfy_cycles(std::vector<smooth_t> &smooths, int cycle_size)
		{
			size_t result = 0;
			auto rend = base_.rend();
			for (auto rit = base_.rbegin(); rit != rend; ++rit)
				if (rit->count == cycle_size)
				{
					int unusable = rit->count;
					for (int i = 0; i < rit->count; i++)
						if (smooths[rit->index[i]].s == smooth_valid_e)
							unusable--;

					if (unusable > 0)
						continue;

					result++;
					smooth_t &pivot = smooths[rit->index[0]];
					for (int i = 1; i < rit->count; i++)
						smooths[rit->index[i]].combine(pivot, n_, base_);
					rit->count = 0; // base element removeable
				}
			return result;
		}
		void compact(std::vector<smooth_t> &smooths)
		{
			// TODO smooth indexes could be remapped, not rebuilt
			// but the base.index has to be rebuilt anyway 
			// since new base element could drop in max_base_count range
#if 0
			std::vector<int> smooth_indexes;
			size_t smooths_size = smooths_.size();
			int smooth_index = 0;
			for (size_t i = 0; i < smooths_size; i++)
			{
				smooth_indexes.push_back(smooth_index);
				smooth_index += (smooths_[i].s != smooth_delete_e);
			}
#endif
			auto smooth_end = std::remove_if(smooths.begin(), smooths.end(), [](const smooth_t &smooth)
			{ return smooth.s == smooth_delete_e; });
			smooths.erase(smooth_end, smooths.end());

			std::vector<int> base_indexes;
			size_t base_size = base_.size();
			int base_index = 0;
			for (size_t i = 0; i < base_size; i++)
			{
				base_indexes.push_back(base_index);
				base_index += (base_[i].count > 1);
			}
			auto new_end = std::remove_if(base_.begin(), base_.end(), [](const base_t &base)
			{ return base.count < 2; });
			base_.erase(new_end, base_.end());

			// now rebuild smooth indexes
			size_t base_slots = base_.size() / small_uint_bits + 1; // leave 1 extrta bit
			for (auto &smooth : smooths)
				smooth.remap(base_indexes, base_slots);
#if 0
			for (auto &base : base_)
				base.remap(smooth_indexes);
#endif
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
#if DBG_SIEVE >= DBG_SIEVE_DEBUG
		std::vector<large_int> ns_;
#endif
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
