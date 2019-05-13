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

#define HAVE_MULTIPLE

namespace zn
{
	template <class large_int, class small_int, class real>
	class quadratic_sieve_t
	{
	public:
		typedef long long long_t;
		typedef unsigned int slot_t;
		enum {bits_per_slot = 8 * sizeof(slot_t) };
		class sieve_range_t
		{
		public:
			sieve_range_t(void) : begin_(0), size_(0), a_(0), b_(0), c_(0), m_(0), dir_sign_(0) {}
			sieve_range_t(const small_int &size, 
						  const large_int &m, 
				          const large_int &c,
						  const large_int &b = large_int(0),
						  const large_int &a = large_int(1)) : size_(size), m_(m), a_(a), b_(b), c_(c), dir_sign_(1)
			{
				large_int d = b * b - 4 * a *c;
				order_ = safe_cast<small_int>(-c_ / m_);
				begin_ = safe_cast<large_int>((-b + sqrt(d) / (2 * a)) + 1);
				step_ = 0;
			}
			large_int first(void) const { return begin_; }
			small_int second(void) const { return size_; }
			bool      sign(void) const { return dir_sign_ < 0; }
			int		  step(void) const { return step_; }
			large_int module(void) const { return m_;}
			small_int order(void) const { return order_; }
			sieve_range_t negate(void) const
			{
				sieve_range_t result = *this;
				result.dir_sign_ = -dir_sign_;
				++result;
				result.step_ = 0;
				return result;
			}
			sieve_range_t &operator++(void)
			{
				begin_ += dir_sign_ * size_;
				step_++;
				return *this;
			}
			large_int mid_value(void) const
			{
				return eval(size_ / 2);
			}
			sieve_range_t first_half(void) const
			{
				sieve_range_t result = *this;
				result.size_ /= 2;
				return result;
			}
			sieve_range_t second_half(void) const
			{
				sieve_range_t result = *this;
				small_int s1 = size_ / 2;
				result.begin_ += s1;
				result.size_ -= s1;
				return result;
			}
			large_int eval(const small_int dx) const
			{
				return eval_base(begin_ + dx);
			}
			std::pair<large_int, large_int> eval_pair(const small_int &dx) const
			{
				std::pair<large_int, large_int> result;
				result.second = begin_ + dx;
				result.first = eval_base(result.second);
				return result;
			}
			large_int eval_base(const large_int &x) const
			{
				large_int y = (a_ * x + b_) % m_;
				y = (y * x + c_) % m_;
				if (dir_sign_ > 0)
					return y;
				else
					return -y;
			}
			large_int eval_plain(const large_int &x) const
			{
				return (a_ * x + b_) * x + c_;
			}
		private:
			large_int a_; // polynomial
			large_int b_;
			large_int c_;

			large_int m_; // module

			large_int begin_; // beginning of sieving range
			small_int size_; // size
			small_int order_;
							  // additional info on polynomial, ecc
			int		dir_sign_; // negative grows backwards
			int     step_;
		};
		class range_handler_t
		{
			bool is_square1(small_int n)
			{
				if (n < 4)
					return false; // 1 is OK
				small_int n1 = static_cast<small_int>(std::sqrt(n) + 0.25);
				return n1 * n1 == n;
			}
		public:
			range_handler_t(const large_int &m, small_int base_size) : m_(m) 
			{
				const auto max_mem = system_info_t::memory();
				const auto cores = system_info_t::cores();
				small_int size = std::min<small_int>(static_cast<small_int>(std::pow(base_size, 2.6)), 
					                                 static_cast<small_int>(max_mem / (sizeof(real) * cores)));
#ifdef HAVE_MULTIPLE
				for (small_int i = 1; i < 10; i++)
					if (!is_square1(i))
#else
				for (small_int i = 1; i < 2; i++)
#endif
				{
					sieve_range_t rangep(size, m, -m * i);
					ranges_[rangep.mid_value()] = rangep;
					sieve_range_t rangen = rangep.negate();
					ranges_[rangen.mid_value()] = rangen;
					/*
					 * Small check
					 */
					large_int y0 = rangep.eval(0);
					large_int y1 = rangep.eval(-1);
					if (y0 < 0 || y1 > 0)
					{
						std::cout << "Wrong range at i = " << i << ":\n" 
								  << "y0 = " << y0 << "\n" 
								  << "y1 = " << y1 << std::endl;
						throw std::runtime_error("Polynomial root problem");
					}
				}
			}
			sieve_range_t next(void)
			{
				auto it = ranges_.begin();
				auto result = it->second;
				ranges_.erase(it);
				auto next = result;
				++next;
				ranges_[next.mid_value()] = next;
				return result;
			}
			sieve_range_t null(void) const
			{
				return sieve_range_t();
			}
		public:
			const large_int m_;
			std::map<large_int, sieve_range_t> ranges_;
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
			bool eval_residue(const large_int &n)
			{
				large_int n1 = n % prime;
				residue = static_cast<small_int>(quadratic_residue<large_int>(n1, prime, prime / prime0)); // actually a power of prime
				return residue != 0 ;
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
		class smooth_t
		{
		public:
			smooth_t(void) : n(1), f(1), sqr(1), sign_bit(false), s(smooth_valid_e) {}

			smooth_t(const sieve_range_t &range,
					 size_t offset, 
				     const large_int &thrs, 
				     const std::vector<base_ref_t> &base) : sqr(1), sign_bit(range.sign()) 
			{
				auto r = range.eval_pair(offset);
				auto m = range.module();
				f = r.first;
				n = r.second;
				//temp_r = r1;
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
						factors_.push_back(static_cast<int>(j));
				}
				if (f == 1)
					s = smooth_valid_e;
				else if (f < thrs)
					s = smooth_candidate_e;
				else
					s = smooth_idle_e;
			}
			bool valid(void) const { return s != smooth_idle_e; }
			bool candidate(void) const { return s == smooth_candidate_e ; }
			bool square(void) const { return factors_.empty(); }
			bool sign_neg(void) const { return sign_bit; }
			large_int reminder(void) const { return f;}
			const std::vector<int>	&factors(void) const { return factors_; }
			void invalidate(void) { s = smooth_idle_e; }
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
				//temp_r = (temp_r * rhs.temp_r) % m;
				if (rhs.f == f)
				{
					sqr = (sqr * f) % m;
					f = 1;
				}

				sign_bit ^= rhs.sign_bit;
				std::vector<int>		fact;
				auto it1  = factors_.begin();
				auto end1 = factors_.end();
				auto it2 = rhs.factors_.begin();
				auto end2 = rhs.factors_.end();
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
				factors_ = fact;
				s = smooth_valid_e;
			}
			void remap_factors(const std::vector<int> &fmap)
			{
				for (auto &idx : factors_)
				{
					idx = fmap[idx];
#if DBG_SIEVE >= DBG_SIEVE_TRACE
					if (idx < 0)
						std::cout << "Factor has been removed...\n";
#endif
				}
			}
			bool invariant(const sieve_range_t &range, const std::vector<base_ref_t> &base) const
			{
				large_int m = range.module();
				large_int s1 = (sqr * sqr * f) % m ;
				for (auto idx : factors_)
					s1 = (s1 * base[idx].prime) % m;
				if (sign_bit)
					s1 = -s1;
				s1 = (s1 - n * n) % m;
				return s1 == 0;
			}
		private:
			large_int				n; // number squared
			std::vector<int>		factors_;
			large_int				f; // remainder after trial division
			//large_int				temp_r; // not necessary
			large_int				sqr;
			bool					sign_bit; // true if negative
			smooth_status_e			s;
		};
		typedef std::vector<smooth_t> smooth_vect_t;
		typedef std::map<large_int, smooth_t> candidates_map_t;

		quadratic_sieve_t(const large_int &n, small_int base_size) : n_(n)
		{
			// double the range; half of them won't be a quadratic residue
			small_int range;
			if (base_size != 0)
#ifdef HAVE_MULTIPLE
				range = primes_range(base_size);
#else
				range = primes_range(base_size * 2);
#endif
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
				small_int n1 = safe_cast<small_int>(n % p);
#ifndef HAVE_MULTIPLE
				if ((r = quadratic_residue(n1, p)) != 0)
#endif
					base_.push_back(base_ref_t(p, r));
			}
#if DBG_SIEVE >= DBG_SIEVE_INFO
			std::cout << "Actual base size: " << base_.size() << ", largest = " << base_.rbegin()->prime << std::endl;
#endif // DBG_SIEVE	
			sieve_thrs_ = safe_cast<real>(std::log(*primes.rbegin()) * 2);
			smooth_thrs_ = base_.rbegin()->prime;
			smooth_thrs_ *= smooth_thrs_;
		}
		large_int sieve(void)
		{
			range_handler_t range_handler(n_, base_.size());

			int count = 0;
			smooth_vect_t smooths;
			candidates_map_t candidates;
#ifdef HAVE_THREADING
			int cores = system_info_t::cores();
			for (int i = 0; i < cores; i++)
			{
				threads_.emplace_back(&quadratic_sieve_t::sieving_thread, this);
				ranges_to_sieve_.push(range_handler.next());
			}
#endif
			for ( ; smooths.size() < base_.size() ; )
			{
#ifdef HAVE_THREADING
				auto chunk = smooths_found_.pop();
#else
				auto chunk = sieve_range(range_handler.next());
#endif
				process_candidate_chunk(candidates, chunk, smooths);
#if DBG_SIEVE >= DBG_SIEVE_INFO
				std::cout << "Found = " << smooths.size() << " smooths (" << count++ << ")\r" << std::flush;
#endif
#ifdef HAVE_THREADING
				ranges_to_sieve_.push(range_handler.next());
#endif
			}
#ifdef HAVE_THREADING
			ranges_to_sieve_.clear();
			for (int i = 0; i < cores; i++)
				ranges_to_sieve_.push(range_handler.null());
			for (auto &thread : threads_)
				thread.join();
			// pick remaining smooths
			while (!smooths_found_.empty())
			{
				auto chunk = smooths_found_.pop();
				process_candidate_chunk(candidates, chunk, smooths);
			}
#endif
#if DBG_SIEVE >= DBG_SIEVE_INFO
			int failed = 0;
			std::cout << "\nFound " << smooths.size() << std::endl;
#endif
			std::cout << "Start erasing" << std::endl;
			erase_base(smooths);
			std::cout << "Erasure done" << std::endl;
			auto result = solve(smooths);
			large_int r = 1;
			for (auto &item : result)
			{
				smooth_t s;
				for (auto index : item)
					s.compose(smooths[index], n_, base_);
#if DBG_SIEVE >= DBG_SIEVE_WARNING
				if (!s.square())
					std::cout << "Non null factors!\n";
#endif // DBG_SIEVE	
				r = s.result(n_);
				if (r != 1 && r != n_)
					break;
#if DBG_SIEVE >= DBG_SIEVE_INFO
				else
					failed++;
#endif // DBG_SIEVE	
			}
#if DBG_SIEVE >= DBG_SIEVE_INFO
			if (failed > 0)
				std::cout << "Failed " << failed << " attempts\n";
#endif // DBG_SIEVE	
					
			return r;
		}
	private:
		void process_candidate_chunk(candidates_map_t &candidates,
									 std::vector<smooth_t> &chunk,
									 std::vector<smooth_t> &smooths)
		{
			for (auto &s : chunk)
			{
				if (s.candidate())
				{
					auto it = candidates.find(s.reminder());
					if (it != candidates.end())
					{
						s.compose(it->second, n_, base_);
#if 0 // DBG_SIEVE >= DBG_SIEVE_ERROR
						if (!s.invariant(n_, base_))
							std::cout << "Hmmm";
#endif // DBG_SIEVE	
					}
					else // just put it aside
					{
						candidates[s.reminder()] = s;
						continue;
					}
				}
				int si = static_cast<int>(smooths.size());
				for (auto f : s.factors())
					base_[f].smooths.push_back(si);
				smooths.push_back(s);
			}
		}
		void erase_base(base_ref_t &base, std::vector<smooth_t> &smooths)
		{
			size_t count = base.smooths.size();
			smooth_t &ref = smooths[base.smooths[0]];
			for (size_t i = 1; i < count; i++)
			{
				smooth_t &value = smooths[base.smooths[i]];
				value.compose(ref);
			}
			ref.invalidate();
		}
		void erase_base(std::vector<smooth_t> &smooths)
		{
			std::vector<int> base_remapping;
			int base_size = static_cast<int>(base_.size());
			int base_map = 0;
#if DBG_SIEVE >= DBG_SIEVE_INFO
			std::cout << "Removing unused bases: start from " << base_size << " and " << smooths.size() << std::endl;
#endif
			for (int i = 0; i < base_size; i++)
			{
				auto &base = base_[i];
				if (base.smooths.size() < 2)
				{
					base_remapping.push_back(-1);
					if (!base.smooths.empty())
						smooths[base.smooths[0]].invalidate();
					continue; // remove it
				}
				base_remapping.push_back(base_map);
				base.smooths.clear();
				if (base_map != i)
					base_[base_map] = base;
				base_map++;
			}
			base_.erase(base_.begin() + base_map, base_.end());
			// compact smooths
			int smooths_size = static_cast<int>(smooths.size());
			int smooths_map = 0;
			for (int i = 0; i < smooths_size; i++)
			{
				auto &smooth = smooths[i];
				if (!smooth.valid())
					continue;
				smooth.remap_factors(base_remapping);
				if (i != smooths_map)
					smooths[smooths_map] = smooth;
				smooths_map++;
			}
			smooths.erase(smooths.begin() + smooths_map, smooths.end());
#if DBG_SIEVE >= DBG_SIEVE_INFO
			std::cout << "Base reduced to " << base_.size() << ", smooth to " << smooths.size()  << std::endl;
#endif
		}
		int find_nonzero(const std::vector<slot_t> &v, int size, int i)
		{
			int slot = static_cast<int>(i / bits_per_slot);
			int j = static_cast<int>(i % bits_per_slot);
			slot_t mask = 1 << j;
			for (i = slot * bits_per_slot ; i < size; i += bits_per_slot, slot++)
			{
				slot_t s = v[slot];
				if (s)
					for (; j < bits_per_slot; j++, mask <<= 1)
						if (s & mask)
							return i + j;
				j = 0;
				mask = 1;
			}
			return -1;
		}
		void swap_bit(std::vector<std::vector<slot_t>> &matrix, int i, int j) // assume j > i
		{
			const int sloti = i / bits_per_slot;
			const int slotj = j / bits_per_slot;
			const size_t size = matrix.size();
			const slot_t di = static_cast<slot_t>(i % bits_per_slot);
			const slot_t dj = static_cast<slot_t>(j % bits_per_slot);
			const slot_t maski = 1 << di ;
			const slot_t maskj = 1 << dj;
			if (sloti == slotj)
			{
				const slot_t delta = dj - di;
				const slot_t mask = maski | maskj;
				const slot_t maskn = ~mask;
				for (size_t k = 0; k < size; k++)
				{
					slot_t x = matrix[k][sloti];
					slot_t mi = (x & maski) << delta ;
					slot_t mj = (x & maskj) >> delta ;
					matrix[k][sloti] = (x &maskn) | mi | mj;
				}
			}
			else
			{
				int delta = dj - di;
				const slot_t maskni = ~maski;
				const slot_t masknj = ~maskj;
				if (delta > 0)
					for (size_t k = 0; k < size; k++)
					{
						slot_t &xi = matrix[k][sloti];
						slot_t &xj = matrix[k][slotj];
						slot_t mi = (xi & maski) << delta;
						slot_t mj = (xj & maskj) >> delta;
						xi = (xi &maskni) | mj;
						xj = (xj &masknj) | mi;
					}
				else
				{
					delta = -delta;
					for (size_t k = 0; k < size; k++)
					{
						slot_t &xi = matrix[k][sloti];
						slot_t &xj = matrix[k][slotj];
						slot_t mi = (xi & maski) >> delta;
						slot_t mj = (xj & maskj) << delta;
						xi = (xi &maskni) | mj;
						xj = (xj &masknj) | mi;
					}
				}
			}
		}
		void trace_matrix(const std::vector<std::vector<slot_t>> &m, size_t hsize)
		{
			int lcount = 0;
			if (hsize > 80)
				return;
			for (auto &v : m)
			{
				size_t scount = 0;
				for (auto s : v)
					for (size_t i = 0; i < bits_per_slot && scount < hsize ; i++, scount++)
						std::cout << (s & (1 << i) ? '1' : '0') << (i % 8 == 7 ? " " : "") ;
				std::cout << "\n" << (++lcount % 8 == 0 ? "\n": "") ;
			}
			std::cout << std::flush;
		}
		// linear system: rows is #base_, cols is number of smooths
		std::vector<std::vector<int>> solve(std::vector<smooth_t> &smooths)
		{
			size_t smooth_size = smooths.size();
			size_t base_size = base_.size() + 1 ;
			size_t slots = (smooth_size + bits_per_slot - 1) / bits_per_slot;
			std::vector<std::vector<slot_t>> matrix(base_size, std::vector<slot_t>(slots, 0));
			for (size_t i = 0 ; i < smooth_size; i++)
			{
				const smooth_t &smooth = smooths[i];
				int slot = static_cast<int>(i / bits_per_slot);
				slot_t mask = 1 << static_cast<slot_t>(i % bits_per_slot);
				if (smooth.sign_neg())
					matrix[0][slot] |= mask;
				for (const auto idx : smooth.factors())
					matrix[idx + 1][slot] |= mask;
			}
#if DBG_SIEVE >= DBG_SIEVE_TRACE
			std::cout << "Initial status\n";
			trace_matrix(matrix, smooth_size);
#endif
			std::vector<int> smooth_perm(smooth_size);
			for (size_t i = 0; i < smooth_size; i++)
				smooth_perm[i] = static_cast<int>(i);
			int i = 0; // actual row used as pivot
			std::vector<int> rows;
			std::vector<int> rows_idx;
			for (size_t k = 0; k < base_size; k++)
			{
				rows_idx.push_back(static_cast<int>(i));
#if DBG_SIEVE >= DBG_SIEVE_TRACE
				std::cout << "\nPass " << i << " out of " << base_size << "\n";
				trace_matrix(matrix, smooth_size);
#endif
				int	   slot = i / bits_per_slot;
				slot_t mask = 1 << static_cast<slot_t>(i % bits_per_slot);
				auto &v = matrix[k];
				int j = find_nonzero(v, static_cast<int>(smooth_size), i);
				if (j < 0)
					continue;
				if (j != i)
				{
					swap_bit(matrix, i, j);
					std::swap(smooth_perm[i], smooth_perm[j]);
#if DBG_SIEVE >= DBG_SIEVE_TRACE
					std::cout << "Swap columns " << i << ", " << j << std::endl;
#endif
				}
				for (j = static_cast<int>(k + 1); j < static_cast<int>(base_size); j++)
					if (matrix[j][slot] & mask)
						for (size_t h = slot; h < slots; h++)
							matrix[j][h] ^= v[h];
				rows.push_back(static_cast<int>(k));
				i++;
			}
			// backsubstitution
			size_t rows_size = rows.size();
			for (size_t c = rows_size - 1; c > 0; c--) // indexing on column, which is != from row, i.e. 'i'
			{
				int i = rows[c];
#if DBG_SIEVE >= DBG_SIEVE_TRACE
				std::cout << "\nBack " << i << " out of " << base_size << "\n";
				trace_matrix(matrix, smooth_size);
#endif
				auto &v = matrix[i];
				int	   slot = static_cast<int>(c / bits_per_slot);
				slot_t mask = 1 << static_cast<slot_t>(c % bits_per_slot);
				if (v[slot] & mask) // may be zero...
					for (int j = i - 1 ; j >= 0 ; j--)
						if (matrix[j][slot] & mask)
						{
							auto &w = matrix[j];
							for (size_t k = slot; k < slots; k++)
								w[k] ^= v[k];
						}
			}
			
			std::vector<std::vector<int>> result;
#if DBG_SIEVE >= DBG_SIEVE_TRACE
			std::cout << "Result:\n";
			trace_matrix(matrix, smooth_size);
#endif
			for (size_t i = rows_size; i < smooth_size; i++)
			{
				std::vector<int> idx;
				slot_t slot = static_cast<slot_t>(i / bits_per_slot);
				slot_t mask = 1 << static_cast<slot_t>(i % bits_per_slot);

				for (size_t j = 0; j < base_size; j++)
					if (matrix[j][slot] & mask)
						idx.push_back(smooth_perm[rows_idx[j]]);  // make use of smooths_perm
				idx.push_back(smooth_perm[i]);
				result.push_back(idx);
			}
			return result;
		}
		bool bit_test(const std::vector<slot_t> &v, size_t index)
		{
			size_t slot = index / bits_per_slot;
			return (v[slot] & (1 << static_cast<slot_t>(index % bits_per_slot))) != 0;
		}
		std::vector<smooth_t> collect_smooth(const sieve_range_t &range,
											 const std::vector<real> &values)
		{
			std::vector<smooth_t> smooths;
			size_t size = static_cast<size_t>(range.second());
			for (size_t i = 0; i < size; i++)
				if (values[i] < sieve_thrs_)
				{
					smooth_t s(range, i, smooth_thrs_, base_);
					if (s.valid())
					{
						smooths.push_back(s);
#if DBG_SIEVE >= DBG_SIEVE_INFO
						if (!s.invariant(range, base_))
							std::cout << "Hmm";
#endif
					}
				}
#if DBG_SIEVE >= DBG_SIEVE_INFO
			std::cout << "Range order = " << range.order() 
				      << ", sign = " << range.sign() 
					  << ", step = " << range.step() 
				      << ", smooths = " << smooths.size() << std::endl;
#endif
#if DBG_SIEVE >= DBG_SIEVE_DEBUG
			if (smooths.size() == 0)
				std::cout << "\nNo smooths on range " << range.order() << "\n"
				          << "Prev value " << range.eval(-1) << "\n"
				          << "First value " << range.eval(0) << "\n"
						  << "Last value " << range.eval(range.second() - 1) << std::endl ;
#endif				
			return smooths;
		}
		smooth_vect_t sieve_range(const sieve_range_t &range)
		{
			std::vector<real> values;
			values.reserve(static_cast<size_t>(range.second()));
			build_sieving_range(range, values);
			sieve_range(values, range, range.first());
			return collect_smooth(range, values);
		}
		void sieve_range(std::vector<real> &values, const sieve_range_t &range, const large_int &begin)
		{
			auto size = values.size();
			for (const auto &base1 : base_)
			{
				base_t base = base1;
				if (!base.eval_residue(-range.eval_plain(0)))
					continue;
				sieve_range(values, range, begin, base);
				if (base.prime != 2)
					sieve_range(values, range, begin, -base);
				base_t powers = base;
				small_int prime_power_end = static_cast<small_int>(std::numeric_limits<int>::max() / base.prime);
				for (int i = 0; (i < 10) && (powers.prime < prime_power_end); i++)
				{
					powers.compose(base, -range.eval_plain(0));
					if (powers.residue == 0)
						break;
					sieve_range(values, range, begin, powers);
					if (powers.prime != 2)
						sieve_range(values, range, begin, -powers);
				}
			}
		}
		void sieve_range(std::vector<real> &values, const sieve_range_t &range, const large_int &begin, const base_t &base)
		{
			large_int n = (begin / base.prime) * base.prime + base.residue;
			if (n < begin)
				n += base.prime;
			auto size = values.size();
			size_t pos = safe_cast<size_t, large_int>(n - begin);
#if DBG_SIEVE >= DBG_SIEVE_INFO
			if ((range.eval_plain(n) % base.prime) != 0)
			{
				std::cout << "Error: n = " << n 
					      << "\np(n) = " << range.eval_plain(n)
					      << "\np = " << base.prime 
					      << "\nr = " << base.residue 
					      << "\np0= " << base.prime0
					      << "\npos= " << pos << ", n-begin=" << (n - begin) 
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
			size_t size = static_cast<size_t>(range.second());
			std::vector<real> data(size);
			for (size_t i = 0; i < size; i++)
			{
				large_int n2 = range.eval(i);
				values.push_back(std::log(safe_cast<real>(n2)));
			}
		}
		void  build_sieving_range(const sieve_range_t &range, std::vector<real> &values)
		{
			large_int n2 = range.eval(0);
			large_int m2 = range.eval(range.second() - 1);
			real rn = std::abs(safe_cast<real>(n2));
			real rm = std::abs(safe_cast<real>(m2));
			real t = rn / rm + rm / rn - 2;
			size_t size = static_cast<size_t>(range.second());
			if (t < 3e-3)
			{
				std::vector<real> data(size);
				rn = std::log(rn);
				rm = std::log(rm);
				real delta = (rm - rn) / size;
				for (size_t i = 0; i < size; i++)
					values.push_back(rn + i * delta);
			}
			else if (size < 32)
				build_sieving_range_exact(range, values);
			else
			{
				build_sieving_range(range.first_half(), values);
				build_sieving_range(range.second_half(), values);
			}
		}
		// removes element that appears only once
		void sieving_thread(void)
		{
			try
			{
				auto range = ranges_to_sieve_.pop();
				for (; range.second() != 0; range = ranges_to_sieve_.pop())
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
