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
#include "znquadratic_sieve_base.h"

namespace zn
{

	template <class large_int, class small_int, class real>
	class multiple_polynomial_quadratic_sieve_t : public quadratic_sieve_base_t<large_int, small_int, real>
	{
	public:
                typedef quadratic_sieve_base_t<large_int, small_int, real> inherit_t ;
                typedef typename inherit_t::base_ref_t base_ref_t ;
                typedef typename inherit_t::smooth_status_e smooth_status_e;

		struct polynomial_seed_t
		{
			std::vector<int> index; // base used for the polynomial
			real target_log;
			polynomial_seed_t(void) : target_log(-1) {}
			bool is_null(void) const { return target_log < -0.5; }
		};
		struct polynomial_t
		{
			std::vector<int> index; // base used for the polynomial
			large_int a0;
			large_int a;  // a = a0 * a0
			large_int b;
			large_int c;
			small_int x1;
			small_int x2;
			bool valid;
			large_int residue(const base_ref_t &base, const large_int &a, const large_int &a2, const large_int &n)
			{
				if (base.powers() > 1)
					return base.residue(1);
				else
					return quadratic_residue<large_int>(n, a2, a); // actually a power of prime
			}
			polynomial_t(const std::vector<int> &idx,
						 const std::vector<base_ref_t> &base,
						 const large_int &n) : index(idx), valid(true)
			{
				large_int r = 0;
				const base_ref_t &bp0 = base[idx[0]];
				a0 = bp0.prime(0);
				a = a0 * a0;
				b = residue(bp0, a0, a, n);
				valid = (b != 0);
				size_t index_count = index.size();
				for (size_t i = 1; valid && (i < index_count); i++)
				{
					const base_ref_t &bp = base[idx[i]];
					large_int p1 = bp.prime(0);
					a0 *= p1;
					large_int g = p1 * p1;
					large_int r1 = residue(bp, p1, g, n);
					if (r1 == 0)
						valid = false;
					large_int db = (r1 - b);
					auto ext = extended_euclidean_algorithm(a, g);
#if DBG_SIEVE >= DBG_SIEVE_TRACE
					if (std::get<0>(ext) != 1)
						throw not_relatively_prime_t<large_int>(db, a);
#endif
					large_int h = (std::get<1>(ext) * db) % g;
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
				compute_zeros(n);
			}

			large_int eval(const large_int &x) const
			{
				large_int x1 = a * x + 2 * b;
				return x1 * x + c;
			}
			real eval_log(small_int x) const
			{
#if 1
				large_int y = abs(eval(x));
				return real_op_t<real>::log1(y);
#else
				return loga + log(abs((x - z1) * (x - z2)) + 1);
#endif
			}
			std::pair<small_int, small_int> zeros(void) const
			{
				return std::make_pair(x1, x2);
			}
			// polynomial is negative in the range of the roots, including bounds
			void compute_zeros(const large_int &n)
			{
				large_int d = safe_cast<large_int>(sqrt(n));
				x1 = safe_cast<small_int>((-b - d) / a);
				x2 = safe_cast<small_int>((-b + d) / a);
#if DBG_SIEVE >= DBG_SIEVE_TRACE
				large_int y1 = eval(x1 - 1);
				large_int y0 = eval(x1);
				if (y1 < 0 || y0 > 0)
					throw std::runtime_error("Error in finding x1");
				y1 = eval(x2 + 1);
				y0 = eval(x2);
				if (y1 < 0 || y0 > 0)
					throw std::runtime_error("Error in finding x2");
#endif
			}
		};

		class polynomial_generator_t
		{
			enum { first_base_e = 2 };
		public:
			polynomial_generator_t(const large_int &n, small_int m, const std::vector<base_ref_t> &base)
			{
				large_int n1 = 2 * n;
				large_int a2 = safe_cast<large_int>(sqrt(n1)) / m;
				polynomial_seed_t result;
				target_ = real_op_t<real>::log1(a2) / 2;
				index_.target_log = 0;
				index_.index.push_back(first_base_e);
				index_.index.push_back(base.size() - 1);
				coarse_init(base);
			}
			polynomial_seed_t operator()(const std::vector<base_ref_t> &base)
			{
				increment(base);
				while (!within_target(index_))
					increment(base);
				return index_;
			}
		private:
			bool within_target(const polynomial_seed_t &seed)
			{
				return abs(seed.target_log - target_) < 8;
			}
			void increment(const std::vector<base_ref_t> &base)
			{
				size_t order = index_.index.size();
				if (index_.target_log < target_) // short of target, incfrease first element
				{
					size_t idx = order - 2;
					index_.target_log += base[index_.index[idx]].logp();
					index_.index[idx]++;
					index_.target_log -= base[index_.index[idx]].logp();
				}
				else
				{
					size_t idx = order - 1;
					index_.target_log += base[index_.index[idx]].logp();
					index_.index[idx]--;
					index_.target_log -= base[index_.index[idx]].logp();
				}
				auto rit = index_.index.rbegin();
				if (rit[0] != rit[1])
					return;
				
				if (order > 2)
				{
					size_t i = order - 2;
					for (; i > 0; i--)
						if (index_.index[i] + 1 < index_.index[i + 1])
							break;
					index_.index[i]++;
					for (i++; i < order - 1; i++)
						index_.index[i] = index_.index[i - 1] + 1;
					index_.index[order - 1] = base.size() - 1;
				}
				// promote to a larger order
				if (index_.index[order - 2] >= index_.index[order - 1])
				{
					for (size_t i = 0; i < order; i++)
						index_.index[i] = i + first_base_e;
					index_.index.push_back(base.size() - 1);
				}
				coarse_init(base);
			}
			void coarse_init(const std::vector<base_ref_t> &base)
			{
				init_seed(index_, base);
				auto it = base.begin();
				for (size_t idx = index_.index.size() - 1; idx > 0 && !within_target(index_) ; idx--)
					if (index_.target_log < target_) // too light; increase first
					{
						index_.target_log += base[index_.index[idx - 1]].logp();
						it = std::lower_bound(base.begin() + index_.index[idx - 1] + 1, 
							base.begin() + index_.index[idx] - 1,
							target_ - index_.target_log,
							[](const base_ref_t &base, real p)		{
							return -base.logp() < p;
						});
						index_.index[idx - 1] = it - base.begin();
						index_.target_log -= base[index_.index[idx - 1]].logp();
					}
					else // too heavy, lower the high
					{
						index_.target_log += base[index_.index[idx]].logp();
						it = std::upper_bound(base.begin() + index_.index[idx - 1] + 1,
							base.begin() + index_.index[idx] - 1,
							target_ - index_.target_log,
							[](real p, const base_ref_t &base) {
							return p < -base.logp() ;
						});
						index_.index[idx]  = it - base.begin();
						index_.target_log -= base[index_.index[idx]].logp();
					}
			}
			void init_seed(polynomial_seed_t &seed, const std::vector<base_ref_t> &base)
			{
				seed.target_log = 0;
				size_t size = seed.index.size();
				for (size_t i = 0; i < size; i++)
					seed.target_log -= base[seed.index[i]].logp();
			}
			polynomial_seed_t index_;
			real target_;
		};

		class smooth_t
		{
		public:
			smooth_t(void) : axb(1), f(1), sqr(1), sign_bit(false), s(inherit_t::smooth_valid_e) {}

			smooth_t(const polynomial_t &poly,
				small_int x,
				const large_int &thrs,
				const large_int &n,
				const std::vector<base_ref_t> &base) : sqr(poly.a0), sign_bit(false)
			{
				axb = poly.a * x + poly.b;
				f = poly.eval(x);
				if (f < 0)
				{
					f = -f;
					sign_bit = true;
				}
				size_t base_size = base.size();
				for (size_t j = 0; j < base_size; j++)
				{
					const auto &b = base[j];
					int rexp = 0;
					large_int p = b.prime(0);
					int power = 0;
					while (divide_qr1(f, p))
						if (++power % 2 == 0)
							sqr = (sqr * p) % n;
					if (power & 1)
						factors_.push_back(static_cast<int>(j));
				}
				if (f == 1)
					s = inherit_t::smooth_valid_e;
				else if (f < thrs)
					s = inherit_t::smooth_candidate_e;
				else
					s = inherit_t::smooth_idle_e;
			}
			smooth_status_e type(void) const { return s; }
			bool square(void) const { return factors_.empty(); }
			bool sign_neg(void) const { return sign_bit; }
			large_int reminder(void) const { return f; }
			const std::vector<int>	&factors(void) const { return factors_; }
			void invalidate(void) { s = inherit_t::smooth_idle_e; }
			large_int result(const large_int &n)
			{
				large_int a = gcd(n, axb + sqr);
				if (a != 1 && a != n)
					return a;
				large_int b = (axb > sqr ? axb - sqr : sqr - axb);
				return gcd(b, n);
			}
			void compose(const smooth_t &rhs, const large_int &n, const std::vector<base_ref_t> &base)
			{
				axb = (axb * rhs.axb) % n;
				sqr = (sqr * rhs.sqr) % n;
				if (rhs.f == f)
				{
					sqr = (sqr * f) % n;
					f = 1;
				}

				sign_bit ^= rhs.sign_bit;
				std::vector<int>		fact;
				auto it1 = factors_.begin();
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
						sqr = (sqr * base[*it1].prime(0)) % n;
						++it1;
						++it2;
					}
				}
				if (it1 != end1)
					fact.insert(fact.end(), it1, end1);
				if (it2 != end2)
					fact.insert(fact.end(), it2, end2);
				factors_ = fact;
				s = inherit_t::smooth_valid_e;
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
			bool invariant(const large_int &n, const std::vector<base_ref_t> &base) const
			{
				large_int s1 = (sqr * sqr * f) % n;
				for (auto idx : factors_)
					s1 = (s1 * base[idx].prime(0)) % n;
				if (sign_bit)
					s1 = -s1;
				s1 = (s1 - axb * axb) % n;
				return s1 == 0;
			}
		private:
			std::vector<int>		factors_;
			large_int				f; // remainder after trial division
			large_int				sqr; // product of prime with even exponents (/ 2)
			large_int				axb;// a *x + b, thenumber squared
			bool					sign_bit; // true if negative
			smooth_status_e			s;
		};
		typedef std::map<large_int, smooth_t> candidates_map_t;
		typedef std::vector<smooth_t> smooth_vector_t;

		multiple_polynomial_quadratic_sieve_t(const large_int &n, const large_int &m, small_int base_size) : n_(n), m_(m)
		{
			// double the range; half of them won't be a quadratic residue
			small_int range;
			if (base_size != 0)
				range = inherit_t::primes_range(base_size * 2);
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
				typename inherit_t::base_t base = inherit_t::base_t::build(p, n);
				if (base.prime_.size() > 0)
					base_.push_back(base);
			}
			small_int largest_sieving_prime = base_.rbegin()->prime(0);
#if DBG_SIEVE >= DBG_SIEVE_INFO
			std::cout << "Actual base size: " << base_.size() << ", largest = " << largest_sieving_prime << std::endl;
#endif // DBG_SIEVE	
			sieve_thrs_ = real_op_t<real>::log1(largest_sieving_prime * 2);
			smooth_thrs_ = largest_sieving_prime;
			smooth_thrs_ *= smooth_thrs_;
			if (m_ < sqrt(largest_sieving_prime))
				m_ = largest_sieving_prime * 2;
		}
		large_int process(void)
		{
			candidates_map_t candidates;
			smooth_vector_t smooths;
			polynomial_generator_t generator(n_, m_, base_);

//			for (int i = 0; i < 10000; i++)
//				generator(base_);

#ifdef HAVE_THREADING
			int cores = system_info_t::cores();
			for (int i = 0; i < cores; i++)
			{
				threads_.emplace_back(&multiple_polynomial_quadratic_sieve_t::sieving_thread, this);
				polynomials_.push(generator(base_));
				polynomials_.push(generator(base_));
			}
#else
				polynomials_.push_back(generator(base_));
#endif

			int count = 0;
			size_t actual_bsize = 0;
			while (smooths.size() < actual_bsize + 5 + actual_bsize / 100)
			{
#ifdef HAVE_THREADING
				auto chunk = smooths_found_.pop();
				polynomials_.push(generator(base_));
#else
				auto item = *polynomials_.begin();
				polynomials_.pop_front();
				polynomials_.push_back(generator(base_));
				polynomial_t p(item.index, base_, n_);
				if (!p.valid)
					continue;
				auto chunk = sieve(p);
#endif
				inherit_t::process_candidates_chunk(base_, candidates, chunk, smooths, n_);
				count++;
				actual_bsize = inherit_t::actual_base_size(base_);
#if DBG_SIEVE >= DBG_SIEVE_INFO
				std::cout << count << ". Sieving.. " << smooths.size() << ", candidates " << candidates.size() << ", base = " << actual_bsize << "\r" << std::flush;
#endif
			}
#ifdef HAVE_THREADING
			polynomials_.clear();
			for (int i = 0; i < cores ; i++)
				polynomials_.push(polynomial_seed_t());
			for (auto &thread : threads_)
				thread.join();
			// pick remaining smooths
			while (!smooths_found_.empty())
			{
				auto chunk = smooths_found_.pop();
				inherit_t::process_candidates_chunk(base_, candidates, chunk, smooths, n_);
			}
#endif
#if DBG_SIEVE >= DBG_SIEVE_INFO
			std::cout << std::endl;
#endif
			if (smooths.size() < actual_bsize)
			{
#if 1 // DBG_SIEVE >= DBG_SIEVE_ERROR
				std::cout << "Found only " << smooths.size() << std::endl;
#endif
				return 1;
			}
#ifdef HAVE_CANDIDATE_ANALYSYS
			inherit_t::print_analysis(base_.rbegin()->prime(0));
#endif
			inherit_t::erase_base(base_, smooths);
			linear_solver_t solver;
			auto basemix = solver.solve(smooths, base_.size() + 1);
			return build_solution(smooths, basemix);
		}
	private:
		large_int  build_solution(const smooth_vector_t &smooths,
								  const std::vector<std::vector<int>> &basemix)
		{
			int failed = 0;
			large_int r = 1;
			for (auto &item : basemix)
			{
				smooth_t s;
				for (auto index : item)
					s.compose(smooths[index], n_, base_);
#if DBG_SIEVE >= DBG_SIEVE_WARNING
				if (!s.square())
					std::cout << "Non null factors!\n";
				if (!s.invariant(n_, base_))
					std::cout << "Hmmmm";
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
		smooth_vector_t sieve(const polynomial_t &poly)
		{
			// build vector for sieving
			std::vector<real> values(safe_cast<size_t>(2 * m_));
			auto zeros = poly.zeros();
			if (-m_ < zeros.first)
			{
				real u = poly.eval_log(zeros.first);
				real v = poly.eval_log(zeros.second);
				fill_range(poly, values, -m_, poly.eval_log(-m_), zeros.first, u);
				fill_range(poly, values, zeros.first, u, zeros.second, v);
				fill_range(poly, values, zeros.second, v, m_, poly.eval_log(m_));
			}
			else // enters here only during tests
				fill_range(poly, values, -m_, poly.eval_log(-m_), m_, poly.eval_log(m_));

			// use the base for sieving
			sieve_values(poly, values);
			return collect_smooth(poly, values);
		}
		smooth_vector_t  collect_smooth(const polynomial_t &poly,
										const std::vector<real> &values)
		{
			std::vector<smooth_t> result;
			real sieve_thrs = -2 * base_.rbegin()->logp_ - real_op_t<real>::unit(); // small prime variation
			large_int largest_prime = base_.rbegin()->prime(0);
			large_int candidate_thrs = largest_prime * largest_prime ;
			size_t size = values.size();
			for (size_t i = 0; i < size; i++)
				if (values[i] < sieve_thrs)
				{
					smooth_t s(poly, i - m_, candidate_thrs, n_, base_);
					if (s.type() != inherit_t::smooth_idle_e)
					{
#if DBG_SIEVE >= DBG_SIEVE_INFO
						if (!s.invariant(n_, base_))
							std::cout << "Hm\n";
#endif // DBG_SIEVE	
						result.push_back(s);
					}
				}
			return result;
		}
		void sieve_values(const polynomial_t &poly, std::vector<real> &values)
		{
			small_int size = values.size();
			for (const auto &base : base_)
			{
				if (base.prime(0) < 20) // small prime variation
					continue;
				size_t powers = base.powers();
				real logp = base.logp();
				for (size_t i = 0; i < powers; i++)
				{
					small_int prime = base.prime(i);
					small_int residue = base.residue(i);
					small_int a = safe_cast<small_int>(poly.a % prime);
					if (a == 0)
						break; // we are hitting a divisor of a
					small_int b = safe_cast<small_int>(poly.b % prime);
					small_int a1 = std::get<1>(extended_euclidean_algorithm<small_int>(a, prime));
					small_int x = ((prime - b + residue) * a1) % prime;
#if 0 //DBG_SIEVE >= DBG_SIEVE_TRACE
					large_int y = poly.eval(x);
					small_int x0 = safe_cast<small_int>(y % prime);
					if (x0 != 0)
						throw std::runtime_error("Sieving with wrong offset");
#endif // DBG_SIEVE	
					// here x = 0 maps to values[m_], so x -= (m_ % prime) * prime
					size_t p1 = static_cast<size_t>(prime);
					size_t index = static_cast<size_t>(m_ + x) % p1;
#if 0 // DBG_SIEVE >= DBG_SIEVE_TRACE
					large_int y1 = poly.eval(index - m_);
					real t1 = std::log(std::abs(safe_cast<real>(y1)));
					if (abs(t1 - values[index]) > 1e-5)
						throw std::runtime_error("Something wrong in sieving");
#endif // DBG_SIEVE	
					for (size_t i = index; i < size; i += p1)
						values[i] += logp;
					if (prime != 2)
					{
						x = ((2 * prime - b - residue) * a1) % prime;
						index = static_cast<size_t>(m_ + x) % p1;
						for (size_t i = index; i < size; i += p1)
							values[i] += logp;
					}
				}
			}
		}
		void fill_range(const polynomial_t &poly,
						std::vector<real> &values,
						small_int begin, real t_1,
						small_int end,   real t1)
		{
			small_int mid = (begin + end) / 2;
			real t0 = poly.eval_log(mid);
			auto q1 = t_1 / static_cast<double>(t0);
			auto q0 = static_cast<double>(t0) / t1;
			auto t = q1 / q0 + q0 / q1 - 2;
			if (t < 0.002)
			{
				fill_linear(values, begin, t_1, mid, t0);
				fill_linear(values, mid, t0, end, t1);
			}
			else
				if (end - begin < 16)
					fill_exact(poly, values, begin, end);
				else
				{
					fill_range(poly, values, begin, t_1, mid, t0);
					fill_range(poly, values, mid, t0, end, t1);
				}
		}
		void fill_linear(std::vector<real> &values,
						const small_int &x1, real t1,
						const small_int &x2, real t2)
		{
			size_t size	  = static_cast<size_t>(x2 - x1);
			size_t offset = static_cast<size_t>(x1 + m_);
			double m = static_cast<double>(t2 - t1) / size;
			for (size_t i = 0; i < size; i++)
				values[i + offset] = static_cast<real>(t1 + m * i);
		}
		void fill_exact(const polynomial_t &poly, std::vector<real> &values, small_int begin, small_int end)
		{
			size_t size = safe_cast<size_t>(end - begin);
			size_t offset = safe_cast<size_t>(begin + m_);
			for (size_t i = 0; i < size; i++)
				values[offset + i] = poly.eval_log(begin + i);
		}
#ifdef HAVE_THREADING
		void sieving_thread(void)
		{
			try
			{
				for (auto seed = polynomials_.pop(); !seed.is_null(); seed = polynomials_.pop())
				{
					polynomial_t p(seed.index, base_, n_);
					if (!p.valid)
						continue;
					auto chunk = sieve(p);
					smooths_found_.push(chunk);
				}
			}
			catch (std::exception &exc)
			{
				std::cout << "Exception in thread: " << exc.what() << std::endl;
			}
		}
#endif
		real						 sieve_thrs_;
		large_int					 smooth_thrs_; // square of last element of base
		large_int					 n_; // number to factor
		small_int					 m_; // size of sieving interval
		std::vector<base_ref_t>		 base_;
#ifdef HAVE_THREADING
		shared_list_t<polynomial_seed_t> polynomials_;
		shared_list_t<smooth_vector_t>   smooths_found_;
		std::vector<std::thread>	     threads_;
#else
		std::list<polynomial_seed_t> polynomials_;
#endif
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

};

#endif
