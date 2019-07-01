//---------------------------------------------------------------------------------
//
//  Zn 
//  Copyright Marco Oman 2019
//
// Distributed under the Boost Software License, Version 1.0. 
// (See accompanying file LICENSE_1_0.txt or copy at 
// http://www.boost.org/LICENSE_1_0.txt)
//
#ifndef zsipqs_H
#define zsipqs_H

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
	class self_initializing_quadratic_sieve_t : public quadratic_sieve_base_t<large_int, small_int, real>
	{
	public:
        typedef quadratic_sieve_base_t<large_int, small_int, real> inherit_t ;
        typedef typename inherit_t::base_ref_t base_ref_t ;
        typedef typename inherit_t::smooth_status_e smooth_status_e;
		typedef polynomial_siqs_t<large_int, small_int> poly_t;
		typedef sieve_range_t<poly_t, real>				sieve_t;

		class smooth_t
		{
		public:
			smooth_t(void) : axb(1), f(1), sqr(1), sign_bit(false), s(inherit_t::smooth_valid_e) {}

			smooth_t(const polynomial_siqs_t<large_int, small_int> &poly,
					small_int idx,
					small_int x,
					const large_int &thrs,
					const large_int &n,
					const std::vector<typename sieve_t::sieve_run_t> &runs) : sqr(poly.a0), sign_bit(false)
			{
				axb = poly.a * x + poly.b;
				f = poly.eval(x);
				if (f < 0)
				{
					f = -f;
					sign_bit = true;
				}
				size_t runs_size = runs.size();
				for (size_t j = 0; j < runs_size; j++)
				{
					const auto &run = runs[j];
#if 1
					if (run.pwr || ((x - run.x) % run.p != 0))
						continue;
#endif
					int rexp = 0;
					large_int p = run.p;
					int power = 0;
					while (divide_qr1(f, p))
						if (++power % 2 == 0)
							sqr = (sqr * p) % n;
					if (power & 1)
						factors_.push_back(static_cast<int>(run.bix));
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
						LOG_ERROR <<  "Factor " << i1 << " has been removed..." << log_base_t::newline_t();
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
			large_int				axb;// a *x + b, the number squared
			bool					sign_bit; // true if negative
			smooth_status_e			s;
		};
		typedef std::map<large_int, smooth_t> candidates_map_t;
		typedef std::vector<smooth_t> smooth_vector_t;
		self_initializing_quadratic_sieve_t(const large_int &n, const large_int &m, small_int base_size, int k = 0) : n_(n), m_(m)
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
				     << ", largest = " << largest_sieving_prime << " ("<< valid_for_a << ")" << log_base_t::newline_t() ;
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
				threads_.emplace_back(&self_initializing_quadratic_sieve_t::sieving_thread, this, std::ref(base_info));
				polynomials_.push(generator());
				polynomials_.push(generator());
			}
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
				auto seed = generator();
				polynomial_siqs_t<large_int, small_int> poly(seed.index, base_info, n_);
				auto chunk = sieve(poly);
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
			LOG_INFO << log_base_t::newline_t();
			LOG_INFO << "Attempted " << smooth_attempts_ << ", failed " << smooth_failures_ << "\n";
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
		smooth_vector_t sieve(polynomial_siqs_t<large_int, small_int> &poly)
		{
			// build vector for sieving
			smooth_vector_t result;
			poly.select(0);
			size_t count = poly.count();
			std::vector<real> values_init = sieve_t(poly, static_cast<size_t>(m_)).fill();
			std::vector<sieve_t::sieve_run_t> runs;
			sieve_t::build_run(poly, base_, runs);
			size_t size = values_init.size();
			for (size_t c = 1; c <= count; c++)
			{
				std::vector<real> values(values_init);
				for (auto &run : runs)
				{
					size_t index = static_cast<size_t>(m_ + run.x) % run.p;
					for (size_t i = index; i < size; i += run.p)
						values[i] += run.lg;
				}
				collect_smooth(poly, values, runs, result);
				if (c < count)
				{
					poly.select(c);
					sieve_t::update_run(poly, runs);
				}
			}
			return result;
		}
		void collect_smooth(const polynomial_siqs_t<large_int, small_int> &poly,
							const std::vector<real> &values, 
							const std::vector<typename sieve_t::sieve_run_t> &runs, 
							smooth_vector_t &result)
		{
			real sieve_thrs = -2 * base_.rbegin()->logp_ - real_op_t<real>::unit(); // small prime variation
			large_int largest_prime = base_.rbegin()->prime(0);
			large_int candidate_thrs = largest_prime * largest_prime ;
			size_t size = values.size();
			for (size_t i = 0; i < size; i++)
				if (values[i] < sieve_thrs)
				{
					smooth_t s(poly, i, i - m_, candidate_thrs, n_, runs);
					smooth_attempts_++;
					
					if (s.type() != inherit_t::smooth_idle_e)
					{
#if DBG_SIEVE >= DBG_SIEVE_INFO
						if (!s.invariant(n_, base_))
							std::cout << "Hm\n";
#endif // DBG_SIEVE	
						result.push_back(s);
					}
					else
						smooth_failures_++;
				}
		}
#ifdef HAVE_THREADING
		void sieving_thread(const std::vector<prime_info_t<small_int>> &base_info)
		{
			try
			{
				for (auto seed = polynomials_.pop(); !seed.is_null(); seed = polynomials_.pop())
				{
					polynomial_siqs_t<large_int, small_int> poly(seed.index, base_info, n_);
					auto chunk = sieve(poly);
					smooths_found_.push(chunk);
				}
			}
			catch (std::exception &exc)
			{
				std::cout << "Exception in thread: " << exc.what() << std::endl;
			}
		}
#endif
		real						sieve_thrs_;
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
#endif
		int smooth_attempts_ = 0;
		int smooth_failures_ = 0;
	};



	template <class large_int, class small_int = int, class real = float>
	large_int self_initializing_quadratic_sieve(const large_int &n, const large_int &m, small_int base_size, int k = 0)
	{
#if DBG_SIEVE >= DBG_SIEVE_INFO
		std::cout << "Factorization of " << n << std::endl;
#endif
		self_initializing_quadratic_sieve_t<large_int, small_int, real> qs(n, m, base_size, k);
		return qs.process();
	}

};

#endif
