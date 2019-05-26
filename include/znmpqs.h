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
		struct polynomial_seed_t
		{
			std::vector<int> index; // base used for the polynomial
			real quality;
			polynomial_seed_t(void) : quality(-1) {}
			bool is_null(void) const { return quality < -0.5; }
			bool operator<(const polynomial_seed_t &rhs) const
			{
				return quality < rhs.quality;
			}
		};
		std::vector<polynomial_seed_t> generate(const std::vector<int> &index, real log1, real log2)
		{
			auto first = base_.begin();
			auto last = base_.begin() + base_.size() - 1;
			std::vector<polynomial_seed_t> result;
			size_t size = index.size() + 2;
			do
			{
				real loga = log1 - first->logp() - last->logp();
				polynomial_seed_t seed;
				seed.index = index;
				seed.index.push_back(first - base_.begin());
				seed.index.push_back(last - base_.begin());
				std::sort(seed.index.begin(), seed.index.end());
				for (size_t i = 1; i < size; i++)
					if (seed.index[i - 1] == seed.index[i])
						seed.quality = -1; // invalid
				seed.quality = abs(loga - log2);
				if (seed.quality < 0.5)
					result.push_back(seed);
				if (loga > log2) // select one smaller than last
					last--;
				else
					first++;
			} while (first != last);
			return result;
		}
		struct polynomial_t
		{
			std::vector<int> index; // base used for the polynomial
			large_int a0; 
			large_int a;  // a = a0 * a0
			large_int b;
			large_int c;
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
				valid = (b != 0) ;
				size_t index_count = index.size();
				for (size_t i = 1 ; valid && (i < index_count); i++)
				{
					const base_ref_t &bp = base[idx[i]];
					large_int p1 = bp.prime(0);
					a0 *= p1;
					large_int g = p1 * p1;
					large_int r1 = residue(bp, p1, g, n);
					if (r1 == 0)
						valid = false;
					large_int db = (r1- b);
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
				large_int x1 = a * x + 2 * b;
				return x1 * x + c;
			}
			// polynomial is negative in the range of the roots, including bounds
			std::pair<large_int, large_int> zeros(const large_int &n) const
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
		class smooth_t
		{
		public:
			smooth_t(void) : axb(1), f(1), sqr(1), sign_bit(false), s(smooth_valid_e) {}

			smooth_t(const polynomial_t &poly,
					small_int x,
					const large_int &thrs,
					const large_int &n,
					const std::vector<base_ref_t> &base) : sqr(poly.a0), sign_bit(false)
			{
				axb = poly.a * x + poly.b ;
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
					s = smooth_valid_e;
				else if (f < thrs)
					s = smooth_candidate_e;
				else
					s = smooth_idle_e;
			}
			smooth_status_e type(void) const { return s; }
			bool square(void) const { return factors_.empty(); }
			bool sign_neg(void) const { return sign_bit; }
			large_int reminder(void) const { return f; }
			const std::vector<int>	&factors(void) const { return factors_; }
			void invalidate(void) { s = smooth_idle_e; }
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
				//temp_r = (temp_r * rhs.temp_r) % m;
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
				if (base.prime_.size() > 0)
					base_.push_back(base);
			}
			small_int largest_sieving_prime = base_.rbegin()->prime(0);
#if DBG_SIEVE >= DBG_SIEVE_INFO
			std::cout << "Actual base size: " << base_.size() << ", largest = " << largest_sieving_prime << std::endl;
#endif // DBG_SIEVE	
			sieve_thrs_ = safe_cast<real>(std::log(largest_sieving_prime) * 2);
			smooth_thrs_ = largest_sieving_prime;
			smooth_thrs_ *= smooth_thrs_;
			if (m_ == 0)
				m_ = largest_sieving_prime * 5;
		}
		large_int process(void)
		{
			candidates_map_t candidates;
			smooth_vector_t smooths;
			build_polynomials();
			size_t smooths_required = base_.size() + 5 + base_.size() / 100;
			int count = 0;
			while (smooths.size() < smooths_required)
			{
#ifdef HAVE_THREADING
				auto chunk = smooths_found_.pop();
#else
				if (polynomials_.empty())
					break;
				auto seed = *polynomials_.begin();
				polynomials_.pop_front();
				polynomial_t p(seed.index, base_, n_);
				if (!p.valid)
					continue;
				auto chunk = sieve(p);
#endif
				process_candidates_chunk(base_, candidates, chunk, smooths, n_);
				count++;
#if DBG_SIEVE >= DBG_SIEVE_INFO
				std::cout << count << ". Sieving.. " << smooths.size() << ", candidates " << candidates.size() << "\r" << std::flush;
#endif
			}
#ifdef HAVE_THREADING
			int cores = system_info_t::cores();
			polynomials_.clear();
			for (int i = 0; i < cores - 1; i++)
				polynomials_.push(polynomial_seed_t());
			for (auto &thread : threads_)
				thread.join();
			// pick remaining smooths
			while (!smooths_found_.empty())
			{
				auto chunk = smooths_found_.pop();
				process_candidates_chunk(base_, candidates, chunk, smooths, n_);
			}
#endif
#if DBG_SIEVE >= DBG_SIEVE_INFO
			std::cout << std::endl;
#endif
			if (smooths.size() < base_.size())
			{
#if DBG_SIEVE >= DBG_SIEVE_ERROR
				std::cout << "Found only " << smooths.size() << std::endl;
#endif
				return 1;
			}
			erase_base(base_, smooths);
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
			auto zeros = poly.zeros(n_);
			if (-m_ < zeros.first)
			{
				fill_range(poly, values, -m_, zeros.first);
				fill_range(poly, values, zeros.first, zeros.second);
				fill_range(poly, values, zeros.second, m_);
			}
			else // enters here only during tests
				fill_range(poly, values, -m_, m_);

			// use the base for sieving
			sieve_values(poly, values);
			//std::sort(values.begin(), values.end());
			return collect_smooth(poly, values);
		}
		smooth_vector_t  collect_smooth(const polynomial_t &poly,
											const std::vector<real> &values)
		{
			std::vector<smooth_t> result;
			real sieve_thrs = -2 * base_.rbegin()->logp_;
			large_int largest_prime = base_.rbegin()->prime(0);
			size_t size = values.size();
			for (size_t i = 0; i < size; i++)
				if (values[i] < sieve_thrs)
				{
					smooth_t s(poly, i - m_, largest_prime * largest_prime, n_, base_);
					if (s.type() != smooth_idle_e)
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
#if DBG_SIEVE >= DBG_SIEVE_TRACE
					large_int y = poly.eval(x);
					small_int x0 = safe_cast<small_int>(y % prime);
					if (x0 != 0)
						throw std::runtime_error("Sieving with wrong offset");
#endif // DBG_SIEVE	
					// here x = 0 maps to values[m_], so x -= (m_ % prime) * prime
					size_t p1 = static_cast<size_t>(prime);
					size_t index = static_cast<size_t>(m_ + x) % p1 ;
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
		void build_polynomials(void)
		{
			large_int n1 = 2 * n_;
			large_int a2 = safe_cast<large_int>(sqrt(n1)) / m_;
			real loga = log(safe_cast<real>(a2)) / 2;
			real factors = -loga / base_.rbegin()->logp();
			std::vector<polynomial_seed_t> seeds;
			if (factors < 1.8)
				seeds = generate(std::vector<int>(), 0, loga);
			else if (factors < 2.8)
			{
				real logp = loga / ceil(factors);
				auto it = std::lower_bound(base_.begin(), base_.end(), logp, [](const base_ref_t &base, real p) {
					return p > -base.logp();
				});
				for (int i = it - base_.begin(); seeds.size() < base_.size() * 3; i++)
				{
					std::vector<int> idx(1, i);
					auto seeds3 = generate(idx, 0, log(safe_cast<real>(a2)) / 2 + it->logp());
					seeds.insert(seeds.end(), seeds3.begin(), seeds3.end());
				}
			}
			else
				throw std::runtime_error("Three factors not enough");
		
			std::sort(seeds.begin(), seeds.end());
#ifdef HAVE_THREADING
			int cores = system_info_t::cores() - 1;
			for (int i = 0; i < cores; i++)
				threads_.emplace_back(&multiple_polynomial_quadratic_sieve_t::sieving_thread, this);
			for (auto &seed : seeds)
				polynomials_.push(seed);
#else
			for (auto &seed : seeds)
				polynomials_.push_back(seed);
#endif
#if DBG_SIEVE >= DBG_SIEVE_INFO
			std::cout << "Polynomial seeds: " << seeds.size() << std::endl;
#endif // DBG_SIEVE	
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
