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
#include "znelliptic_curve_fact.h"


namespace zn
{
#ifdef HAVE_DOUBLE_LARGE_PRIME
	template <class smooth_t>
	class relations_graph_t
	{
	public:
		typedef size_t id_t;
		typedef int color_t;
		typedef typename smooth_t::small_int small_int ;
		struct vertex_t
		{
			color_t				color = 0; // of connected components
			std::list<size_t>	edges; // index into edges array
		};
		typedef std::map<small_int, vertex_t> vertex_map_t;
		typedef typename vertex_map_t::iterator vertex_descriptor_t;

		relations_graph_t(void)
		{
			// insert the '1' vertex
			vertex_t v;
			v.color = 1;
			vertexes_[1] = v;
		}
		size_t size(void) { return edges_.size(); }
		void add_smooth(const smooth_t &smooth)
		{
			auto edge_id = edges_.size();
			edges_.push_back(smooth);
			small_int f0 = smooth.factor(0);
			small_int f1 = (smooth.count() == 1 ? 1 : smooth.factor(1));

			auto it0 = vertexes_.insert(std::make_pair(f0, vertex_t())).first;
			it0->second.edges.push_back(edge_id);
			auto it1 = vertexes_.insert(std::make_pair(f1, vertex_t())).first;
			it1->second.edges.push_back(edge_id);

			if (it0->second.color == 0) // newly inserted
			{
				if (it1->second.color != 0) // extend connected component
					add_connected(it0->second, it1->second.color);
				else
					create_connected(it0->second, it1->second);
			}
			else
			{
				if (it1->second.color == 0)
					add_connected(it1->second, it0->second.color);
				else
					if (it0->second.color != it1->second.color) 
					{ // merge two connected components
						auto cit0 = connected_.find(it0->second.color);
						auto cit1 = connected_.find(it1->second.color);
						// both should be valid
						if (cit0->second < cit1->second)
						{
							connected_.erase(cit0);
							merge_connected_as(it0, cit1->second);
						}
						else
						{
							connected_.erase(cit0);
							merge_connected_as(it1, cit0->second);
						}
					}
					else
					{// found a cycle!!
					}
			}
		}
	private:
		void merge_connected_as(vertex_descriptor_t &v, int color)
		{
			if (v->second.color == color) // already ok
				return;

			v->second.color = color;
			vertex_descriptor_t v1;
			for (auto idx : v->second.edges)
			{
				auto &smooth = edges_[idx];
				small_int f = smooth.factor(smooth.factor(0) == v->first);
				merge_connected_as(vertexes_.find(f), color);
			}
		}
		void add_connected(vertex_t &v, int color)
		{
			v.color = color;
			connected_.find(color)->second++;
		}
		void create_connected(vertex_t &v1, vertex_t &v2)
		{
			v1.color = v2.color = color_++;
			connected_[v1.color] = 2;
		}
		vertex_map_t			vertexes_;
		std::vector<smooth_t>	edges_;
		std::map<int, int>		connected_; // maps id of connected components to its size
		int						color_ = 2; // used to mark connected components
	};
#endif


	template <class large_int, class small_int, class real>
	class self_initializing_quadratic_sieve_t : public quadratic_sieve_base_t<large_int, small_int, real>
	{
	public:
        typedef quadratic_sieve_base_t<large_int, small_int, real> inherit_t ;
        typedef typename inherit_t::base_ref_t base_ref_t ;
        typedef typename inherit_t::smooth_status_e smooth_status_e;
		typedef polynomial_siqs_t<large_int, small_int> poly_t;
		typedef sieve_range_t<poly_t, real>				sieve_t;
		struct smooth_info_t
		{
			poly_t		poly;
			large_int	thrs;
			large_int	n;
			std::vector<typename sieve_t::sieve_run_t> runs;
		};
		class smooth_t
		{
		public:
			typedef typename large_int large_int;
			typedef typename small_int small_int;
			smooth_t(void) : axb(1), sqr(1), sign_bit(false), s(inherit_t::smooth_valid_e) {}

			smooth_t(const smooth_info_t &info, small_int x) : sqr(info.poly.a0), sign_bit(false)
			{
				axb = info.poly.a * x + info.poly.b;
				auto f = info.poly.eval(x);
				if (f < 0)
				{
					f = -f;
					sign_bit = true;
				}
				size_t runs_size = info.runs.size();
				for (size_t j = 0; j < runs_size; j++)
				{
					const auto &run = info.runs[j];
#if 1
					if ((run.bix < 0) || ((x - run.x) % run.p != 0))
					{
#ifdef _DEBUG1
						large_int f1 = f;
						if (divide_qr1(f1, large_int(run.p)))
							std::cout << "Hmm   " << run.p << "\n";
#endif
						continue;
					}
#endif
					int rexp = 0;
					large_int p = run.p;
					int power = 0;
					while (divide_qr1(f, p))
						if (++power % 2 == 0)
							sqr = (sqr * p) % info.n;
#ifdef _DEBUG
					if (power == 0)
						std::cout << "Boh.. " << run.p << "\n";
#endif
					if (power & 1)
						factors_.push_back(run.bix);
				}
				if (f == 1)
					s = inherit_t::smooth_valid_e;
				else if (f < info.thrs)
					s = inherit_t::smooth_candidate_e;
				else
				{
					s = inherit_t::smooth_idle_e;
				}
			}
			//
			//	Accessors
			//
			smooth_status_e type(void) const { return s; }
			int				count(void) const { return static_cast<int>(s);}
			small_int		factor(int index) const { return large_f[index]; }
			small_int		reminder(void) const { return large_f[0]; }

			bool square(void) const { return factors_.empty(); }
			bool sign_neg(void) const { return sign_bit; }
			const std::vector<int>	&factors(void) const { return factors_; }
			//
			//	Operations
			//
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
				for (int i = 0; i < static_cast<int>(s); i++)
					for (int j = 0; j < static_cast<int>(rhs.s); j++)
					{
						sqr = (sqr * large_f[i]) % n;
						large_f[i] = 1;
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
				large_int s1 = (sqr * sqr) % n;
				for (int i = 0; i < static_cast<int>(s); i++)
					s1 = (s1 * large_f[i]) % n;
				for (auto idx : factors_)
					s1 = (s1 * base[idx].prime(0)) % n;
				if (sign_bit)
					s1 = -s1;
				s1 = (s1 - axb * axb) % n;
				return s1 == 0;
			}
			bool large_composite(large_int f) const
			{
				int seeds[] = { 7, 43};
				for (auto s : seeds)
					if (powm<large_int>(s, f - 1, f) != 1)
						return true;
				return false;
			}
		private:
			std::vector<int>		factors_; // indexes into base array
			small_int			    large_f[2];
			large_int				sqr; // product of prime with even exponents (/ 2)
			large_int				axb;// a *x + b, the number squared
			bool					sign_bit; // true if negative

			smooth_status_e			s;
		};
		struct sieve_stuff_t
		{
			std::vector<real> values;
			std::vector<real> init;
			real			  offset;
		};
		typedef std::map<large_int, smooth_t> candidates_map_t;
		typedef std::vector<smooth_t> smooth_vector_t;
		self_initializing_quadratic_sieve_t(const large_int &n, const large_int &m, small_int base_size) : n_(n), m_(m)
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
			smooth_thrs_ = largest_sieving_prime;
			smooth_thrs_ *= smooth_thrs_;
			if (m_ < sqrt(largest_sieving_prime))
				m_ = largest_sieving_prime * 2;
		}
		large_int process(int order)
		{
#ifdef HAVE_DOUBLE_LARGE_PRIME
			relations_graph_t<smooth_t> candidates;
#else
			candidates_map_t candidates;
#endif
			smooth_vector_t smooths;
			std::vector<prime_info_t<small_int>> base_info;
			for (auto &item : base_)
				if (item.powers() > 1)
					base_info.push_back(prime_info_t<small_int>{item.prime(0), item.prime(1), item.residue(1)});
			polynomial_generator_t<large_int, small_int> generator(n_, m_, base_info);
			if (order > 0)
				generator.order_init(order);
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
#else
			polynomial_siqs_t<large_int, small_int> poly(generator().index, base_info, n_);
			sieve_stuff_t sieve_stuff;
			sieve_stuff.init = sieve_t(poly, static_cast<size_t>(m_)).fill(sieve_stuff.offset);
			sieve_stuff.values = sieve_stuff.init;
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
				auto chunk = sieve(poly, sieve_stuff);
#endif
				promoted += process_candidates_chunk(base_, candidates, chunk, smooths, n_);
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
			LOG_INFO << "Attempted " << smooth_attempts_ 
				     << ", failed " << smooth_failures_ 
#ifdef HAVE_DOUBLE_LARGE_PRIME
					  << ", of which "  << smooth_large_composite_ 
					  << " (" << smooth_large_composite_ * 100 / smooth_failures_ << "%) are composite."
#endif // HAVE_DOUBLE_LARGE_PRIME
				     << "\n";
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
		smooth_vector_t sieve(const polynomial_siqs_t<large_int, small_int> &poly,
							  sieve_stuff_t &sieve_stuff)
		{
			// build vector for sieving
			smooth_vector_t result;
			large_int largest_prime = base_.rbegin()->prime(0);
			smooth_info_t info = { poly, largest_prime * largest_prime , n_ };
			info.runs.reserve(base_.size() * 2);
			size_t count = info.poly.count();
			for (size_t c = 0; c < count; c++)
			{
				info.poly.select(c);
				if (c == 0)
					sieve_t::build_run(info.poly, base_, info.runs);
				else
					sieve_t::update_run(info.poly, info.runs);

				std::copy(sieve_stuff.init.begin(), sieve_stuff.init.end(), sieve_stuff.values.begin());
				size_t size = sieve_stuff.values.size();
				for (auto &run : info.runs)
				{
					size_t index = static_cast<size_t>(m_ + run.x) % run.p;
					for (size_t i = index; i < size; i += run.p)
						sieve_stuff.values[i] -= run.lg;
				}
				collect_smooth(info, sieve_stuff, result);
			}
			return result;
		}
		void collect_smooth(const smooth_info_t &info,
							const sieve_stuff_t &sieve,
							smooth_vector_t &result)
		{
#ifdef HAVE_DOUBLE_LARGE_PRIME
			real sieve_thrs = 3 * base_.rbegin()->logp_ + sieve.offset + 5 * real_op_t<real>::unit(); // small prime variation
#else
			real sieve_thrs = 2 * base_.rbegin()->logp_ + sieve.offset + 2 * real_op_t<real>::unit(); // small prime variation
#endif
			size_t size = sieve.values.size();
			for (size_t i = 0; i < size; i++)
				if (sieve.values[i] < sieve_thrs)
				{
					smooth_t s(info, i - m_);
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
					{
#ifdef HAVE_DOUBLE_LARGE_PRIME
						if (s.large_composite())
							smooth_large_composite_++;
#endif
						smooth_failures_++;
					}
				}
		}
#ifdef HAVE_THREADING
		void sieving_thread(const std::vector<prime_info_t<small_int>> &base_info)
		{
			sieve_stuff_t sieve_stuff;
			try
			{
				for (auto seed = polynomials_.pop(); !seed.is_null(); seed = polynomials_.pop())
				{
					polynomial_siqs_t<large_int, small_int> poly(seed.index, base_info, n_);
					if (sieve_stuff.init.empty())
					{
						sieve_stuff.init = sieve_t(poly, static_cast<size_t>(m_)).fill(sieve_stuff.offset);
						sieve_stuff.values = sieve_stuff.init;
					}
					auto chunk = sieve(poly, sieve_stuff);
					smooths_found_.push(chunk);
				}
			}
			catch (std::exception &exc)
			{
				std::cout << "Exception in thread: " << exc.what() << std::endl;
			}
		}
#endif
#ifdef HAVE_DOUBLE_LARGE_PRIME
		int  process_candidates_chunk(std::vector<base_ref_t> &base,
 									  relations_graph_t<smooth_t>  &relations,
									  std::vector<smooth_t> &chunk,
									  std::vector<smooth_t> &smooths,
									  const large_int &n)
		{
			for (auto &smooth : chunk)
			{
				relations.add_smooth(smooth);
			}
			return 0;
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
#endif
		int smooth_attempts_ = 0;
		int smooth_failures_ = 0;
#ifdef HAVE_DOUBLE_LARGE_PRIME
		int smooth_large_composite_ = 0;
#endif
	};



	template <class large_int, class small_int = int, class real = float>
	large_int self_initializing_quadratic_sieve(const large_int &n, const large_int &m, small_int base_size, int order = 0)
	{
#if DBG_SIEVE >= DBG_SIEVE_INFO
		std::cout << "Factorization of " << n << std::endl;
#endif
		self_initializing_quadratic_sieve_t<large_int, small_int, real> qs(n, m, base_size);
		return qs.process(order);
	}

};

#endif
