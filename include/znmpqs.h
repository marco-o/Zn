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
#include "znlinear_solver.h"
#include "znpolynomial.h"


namespace zn
{

	template <class large_int, class small_int, class real>
	class multiple_polynomial_quadratic_sieve_t : public quadratic_sieve_base_t<large_int, small_int, real>
	{
	public:
		typedef quadratic_sieve_base_t<large_int, small_int, real> inherit_t;
		typedef typename inherit_t::base_ref_t base_ref_t;
		typedef typename inherit_t::smooth_status_e smooth_status_e;

		class smooth_t
		{
		public:
			smooth_t(void) : axb(1), f(1), sqr(1), sign_bit(false), s(inherit_t::smooth_valid_e) {}

			smooth_t(const polynomial_t<large_int, small_int> &poly,
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
			void compose(const smooth_t &rhs,
				const large_int &n,
				const std::vector<base_ref_t> &base,
				std::vector<int> *erased = nullptr,
				std::vector<int> *added = nullptr)
			{
				axb = (axb * rhs.axb) % n;
				sqr = (sqr * rhs.sqr) % n;
				if ((rhs.f == f) && (f != 1))
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
					{
						if (added)
							added->push_back(*it2);
						fact.push_back(*it2++);
					}
					else // same factor; gets skipped!
					{
						if (erased)
							erased->push_back(*it1);
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
					int i1 = idx;
					idx = fmap[idx];
#if DBG_SIEVE >= DBG_SIEVE_TRACE
					if (idx < 0)
						LOG_ERROR << "Factor " << i1 << " has been removed..." << log_base_t::newline_t();
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
		multiple_polynomial_quadratic_sieve_t(const large_int &n, const large_int &m, small_int base_size, int k = 0) : n_(n), m_(m)
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
#ifdef HAVE_MULTIPLIER
			if (k == 0)
				k_ = premultiplier(n, primes);
			else
				k_ = k;
			LOG_INFO << "Premultiplier = " << k_ << log_base_t::newline_t();
			n_ *= k_;
#endif
			int valid_for_a = 0;
			for (auto p : primes)
			{
				typename inherit_t::base_t base = inherit_t::base_t::build(p, n_);
				if (base.prime_.size() > 0)
				{
					base_.push_back(base);
					if (base.powers() > 1)
						valid_for_a++;
				}
			}
			small_int largest_sieving_prime = base_.rbegin()->prime(0);
			LOG_INFO << "Actual base size: " << base_.size()
				<< ", largest = " << largest_sieving_prime << " (" << valid_for_a << ")" << log_base_t::newline_t();
			smooth_thrs_ = largest_sieving_prime;
			smooth_thrs_ *= smooth_thrs_;
			if (m_ < sqrt(largest_sieving_prime))
				m_ = largest_sieving_prime * 2;
		}
		large_int process(void)
		{
			candidates_map_t candidates;
			smooth_vector_t smooths;
			std::vector<prime_info_t<small_int>> base_info;
			for (auto &item : base_)
				if (item.powers() > 1)
					base_info.push_back(prime_info_t<small_int>{item.prime(0), item.prime(1), item.residue(1)});
			polynomial_generator_t<large_int, small_int> generator(n_, m_, base_info);

#ifdef HAVE_TIMING
			time_estimator_t time_estimator(base_.size());
#endif

#ifdef HAVE_THREADING
			int cores = system_info_t::cores();
			for (int i = 0; i < cores; i++)
			{
				threads_.emplace_back(&multiple_polynomial_quadratic_sieve_t::sieving_thread, this, std::ref(base_info));
				polynomials_.push(generator());
				polynomials_.push(generator());
			}
#else
			polynomials_.push_back(generator());
#endif

			int count = 0;
			int promoted = 0;
			size_t actual_bsize = 0;
			while (smooths.size() < actual_bsize + 5 + actual_bsize / 100)
			{
#ifdef HAVE_THREADING
				auto chunk = smooths_found_.pop();
				polynomials_.push(generator());
#else
				auto item = *polynomials_.begin();
				polynomials_.pop_front();
				polynomials_.push_back(generator());
				polynomial_t<large_int, small_int> p(item.index, base_info, n_);
				if (!p.valid)
					continue;
				auto chunk = sieve(p);
#endif
				promoted += inherit_t::process_candidates_chunk(base_, candidates, chunk, smooths, n_);
				count++;
				actual_bsize = inherit_t::actual_base_size(base_);
#ifdef HAVE_TIMING
				time_estimator.update(static_cast<int>(smooths.size()), promoted);
#endif
				LOG_INFO << count << ". Sieving.. " << smooths.size()
					<< ", candidates " << candidates.size()
					<< ", base = " << actual_bsize
#ifdef HAVE_TIMING
					<< " (" << time_estimator.elapsed()
					<< " - " << time_estimator.estimated() << ")"
#endif					     
					<< "   \r" << log_base_t::flush_t();
			}
#ifdef HAVE_THREADING
			polynomials_.clear();
			for (int i = 0; i < cores; i++)
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
			LOG_INFO << log_base_t::newline_t();
			if (smooths.size() < actual_bsize)
			{
				LOG_ERROR << "Found only " << smooths.size() << log_base_t::newline_t();
				return 1;
			}
#ifdef HAVE_CANDIDATE_ANALYSYS
			inherit_t::print_analysis(base_.rbegin()->prime(0));
#endif
			inherit_t::erase_base(base_, smooths, n_);
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
#ifdef HAVE_MULTIPLIER
			large_int n1 = n_ / k_;
#endif
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
#ifdef HAVE_MULTIPLIER
				r = s.result(n1);
				if (r != 1 && r != n1)
					break;
#else
				r = s.result(n_);
				if (r != 1 && r != n_)
					break;
#endif
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
		smooth_vector_t sieve(const polynomial_t<large_int, small_int> &poly)
		{
			// build vector for sieving
			std::vector<real> values(safe_cast<size_t>(2 * m_));
			auto zeros = poly.zeros();
			if (-m_ < zeros.first)
			{
				real u = poly.eval_log<real>(zeros.first);
				real v = poly.eval_log<real>(zeros.second);
				fill_range(poly, values, -m_, poly.eval_log<real>(-m_), zeros.first, u);
				fill_range(poly, values, zeros.first, u, zeros.second, v);
				fill_range(poly, values, zeros.second, v, m_, poly.eval_log<real>(m_));
			}
			else // enters here only during tests
				fill_range(poly, values, -m_, poly.eval_log<real>(-m_), m_, poly.eval_log<real>(m_));

			// use the base for sieving
			sieve_values(poly, values);
			return collect_smooth(poly, values);
		}
		smooth_vector_t  collect_smooth(const polynomial_t<large_int, small_int> &poly,
			const std::vector<real> &values)
		{
			std::vector<smooth_t> result;
			real sieve_thrs = -2 * base_.rbegin()->logp_ - real_op_t<real>::unit(); // small prime variation
			large_int largest_prime = base_.rbegin()->prime(0);
			large_int candidate_thrs = largest_prime * largest_prime;
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

		void sieve_values(const polynomial_t<large_int, small_int> &poly, std::vector<real> &values)
		{
			size_t size = values.size();
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
		void fill_range(const polynomial_t<large_int, small_int> &poly,
						std::vector<real> &values,
						small_int begin, real t_1,
						small_int end,   real t1)
		{
			small_int mid = (begin + end) / 2;
			real t0 = poly.eval_log<real>(mid);
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
		void fill_exact(const polynomial_t<large_int, small_int> &poly, std::vector<real> &values, small_int begin, small_int end)
		{
			size_t size = safe_cast<size_t>(end - begin);
			size_t offset = safe_cast<size_t>(begin + m_);
			for (size_t i = 0; i < size; i++)
				values[offset + i] = poly.eval_log<real>(begin + i);
		}
#ifdef HAVE_THREADING
		void sieving_thread(const std::vector<prime_info_t<small_int>> &base_info)
		{
			try
			{
				for (auto seed = polynomials_.pop(); !seed.is_null(); seed = polynomials_.pop())
				{
					polynomial_t<large_int, small_int> p(seed.index, base_info, n_);
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
		large_int					smooth_thrs_; // square of last element of base
		large_int					n_; // number to factor
		small_int					m_; // size of sieving interval
#ifdef HAVE_MULTIPLIER
		small_int					k_;
#endif
		std::vector<base_ref_t>		base_;
#ifdef HAVE_THREADING
		shared_list_t<polynomial_seed_t> polynomials_;
		shared_list_t<smooth_vector_t>   smooths_found_;
		std::vector<std::thread>	     threads_;
#else
		std::list<polynomial_seed_t> polynomials_;
#endif
	};



	template <class large_int, class small_int = int, class real = float>
	large_int multiple_polynomial_quadratic_sieve(const large_int &n, const large_int &m, small_int base_size, int k = 0)
	{
#if DBG_SIEVE >= DBG_SIEVE_INFO
		std::cout << "Factorization of " << n << std::endl;
#endif
		multiple_polynomial_quadratic_sieve_t<large_int, small_int, real> qs(n, m, base_size, k);
		return qs.process();
	}

};

#endif
