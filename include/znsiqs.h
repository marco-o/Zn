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
	bool rb_test(const boost::multiprecision::cpp_int &);
#ifdef HAVE_CANDIDATE_ANALYSYS
	struct analysis_t
	{
		int smooth_attempts = 0;
		int smooth_idle = 0;
#ifdef HAVE_DOUBLE_LARGE_PRIME
		int smooth_large_composite = 0;
		int smooth_unfactored = 0;
		int smooth_huge_prime = 0;
#endif
	};
#endif

	template <class smooth_t>
	class relations_graph_t
	{
	public:
		typedef size_t id_t;
		typedef int color_t;
		typedef typename smooth_t::small_int small_int;
		typedef typename smooth_t::large_int large_int;
		struct vertex_t
		{
			color_t				color = 0; // of connected components
			size_t				dist ;  // from '1'
			std::list<size_t>	edges; // index into edges array
			vertex_t(int c = 0, size_t d = std::numeric_limits<size_t>::max() / 2) : color(c), dist(d) {}
		};
		typedef std::map<small_int, vertex_t> vertex_map_t;
		typedef typename vertex_map_t::iterator vertex_descriptor_t;

		relations_graph_t(const large_int &n) : n_(n)
		{
			// insert the '1' vertex
			vertex_t v(1, 0);
			vertexes_[1] = v;
			connected_[1] = 1;
		}
		size_t size(void) { return edges_.size(); }
		template <class base_ref_t>
		bool add_smooth(smooth_t &smooth, const std::vector<base_ref_t> &base)
		{
			// for single large prime factor(1) is 1
			auto it0 = vertexes_.insert(std::make_pair(smooth.factor(0), vertex_t())).first;
			auto it1 = vertexes_.insert(std::make_pair(smooth.factor(1), vertex_t())).first;

			if (it0->second.color == 0) // newly inserted
			{
				if (it1->second.color != 0) // extend connected component
					add_connected(it0->second, it1->second);
				else
					create_connected(it0->second, it1->second);
			}
			else
			{
				if (it1->second.color == 0)
					add_connected(it1->second, it0->second);
				else
					if (it0->second.color != it1->second.color) 
					{ // merge two connected components
						auto cit0 = connected_.find(it0->second.color);
						auto cit1 = connected_.find(it1->second.color);
						// both should be valid
						if (cit0->second < cit1->second)
						{
							int expected = cit0->second + cit1->second;
							connected_.erase(cit0);
							merge_connected_as(it0, cit1, it1->second.dist + 1);
							if (cit1->second != expected)
								std::cout << "Hmm\n";
						}
						else
						{
							int expected = cit0->second + cit1->second;;
							connected_.erase(cit1);
							merge_connected_as(it1, cit0, it0->second.dist + 1);
							if (cit0->second != expected)
								std::cout << "Hmm\n";
						}
					}
					else // found a cycle!!
					{	 // the problem is which are the edges involved
						if (it0->second.color == 1)
						{
							path_to1(smooth, it0, base, true);
							path_to1(smooth, it1, base, false);
							if (smooth.type() != 0)
								std::cout << "Bad cycle\n";
						}
						else
						{
							std::cout << "TODO: cycle in Oort belt\n";
							return false;
						}
						return true;
					}
			}
			auto edge_id = edges_.size();
			it0->second.edges.push_back(edge_id);
			it1->second.edges.push_back(edge_id);
			edges_.push_back(smooth);
			return false;
		}
	private:
		template <class base_ref_t>
		void path_to1(smooth_t &smooth, 
			          vertex_descriptor_t v, 
			          const std::vector<base_ref_t> &base,
			          bool pre)
		{
			if (v->second.dist > 0)
			{
				for (auto idx : v->second.edges)
				{
					smooth_t &edge = edges_[idx];
					small_int f = edge.factor(edge.factor(0) == v->first);
					auto v1 = vertexes_.find(f);
					if (v1->second.dist < v->second.dist)
					{
						if (pre)
						{
							smooth.compose(edge, n_, base);
							path_to1(smooth, v1, base, pre);
						}
						else
						{
							path_to1(smooth, v1, base, pre);
							smooth.compose(edge, n_, base);
						}
						return;
					}
				}
			}
		}
		void merge_connected_as(vertex_descriptor_t &v, std::map<int, int>::iterator &iter, size_t dist)
		{
			if (v->second.color == iter->first) // already ok
			{
				if (v->second.dist > dist)
					v->second.dist = dist;
				return;
			}

			v->second.color = iter->first;
			v->second.dist = dist;
			iter->second++;
			vertex_descriptor_t v1;
			for (auto idx : v->second.edges)
			{
				auto &smooth = edges_[idx];
				small_int f = smooth.factor(smooth.factor(0) == v->first);
				merge_connected_as(vertexes_.find(f), iter, dist + 1);
			}
		}
		void add_connected(vertex_t &v, const vertex_t &w)
		{
			v.color = w.color;
			v.dist = w.dist + 1;
			connected_.find(w.color)->second++;
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
		large_int				n_;
	};

	struct qs_config_t
	{
		int multiplier = 0; 
		int order = 2; // numer of primes in 'a'
		int base_size = 0; // number of primes (approximate) in base
	};

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
			large_int	thrs3;
			large_int	n;
			bool		have_double;
			std::vector<typename sieve_t::sieve_run_t> runs;
		};
		class smooth_t
		{
		public:
			typedef typename large_int large_int;
			typedef typename small_int small_int;
			smooth_t(void) : axb(1), sqr(1), sign_bit(false), s(inherit_t::smooth_valid_e) 
			{
				large_f[0] = large_f[1] = 1;
			}

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
					bool plus = ((x - run.x[0]) % run.p == 0);
					bool minus = ((x - run.x[1]) % run.p == 0);
					//if (run.a != 0)
					if ((run.bix < 0) || !(plus || minus))
					{
#ifdef _DEBUG
						large_int f1 = f;
						if (divide_qr1(f1, large_int(run.p)))
							std::cout << "Missed " << run.p 
							          << ", b = " << run.b 
							          << ", c = " << run.r[0] 
							          << "\n";
#endif
						continue;
					}
#endif
#if 0
					if (run.p == 17)
					{
						std::cout << "p = " << run.p 
								  << ", a = " << safe_cast<small_int>(info.poly.a % run.p)
							      << ", b = " << safe_cast<small_int>(info.poly.b % run.p)
							      << ", c = " << safe_cast<small_int>(info.poly.c % run.p) << "\n";
						for (int i = 0; i < run.p; i++)
						{
							auto f1 = info.poly.eval(i);
							std::cout << "x = " << i << ", y = " << safe_cast<small_int>(f1 % run.p) << "\n";
						}
						std::cout << std::endl;
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
						if (run.a)
						{
							auto r2 = safe_cast<small_int>(info.n % run.p);
							auto r1 = run.r[0] * run.r[0] % run.p;
							std::cout << "Expected " << run.p << "\n";
						}
						else
						{
							std::cout << "p = " << run.p
								<< ", a = " << safe_cast<small_int>(info.poly.a % run.p)
								<< ", b = " << safe_cast<small_int>(info.poly.b % run.p)
								<< ", c = " << safe_cast<small_int>(info.poly.c % run.p) << "\n";
							std::cout << "f = 2bx + c?\n";
							std::cout << "p = " << run.p << ", b = " << run.b << " c = " << run.r[0] << "\n";
							for (int i = 0; i < run.p ; i++)
							{
								auto fi = info.poly.eval(i);
								std::cout << "x = " << i << ", f = " << safe_cast<small_int>(fi % run.p) << "\n";
							}
							std::cout << std::endl;
						}
//					if (run.a == 0 && power > 0)
//						std::cout << "Works\n";
#endif
					if (power & 1)
						factors_.push_back(run.bix);
				}
				s = inherit_t::smooth_idle_e;
				if (f == 1)
					s = inherit_t::smooth_valid_e;
				else if (f < info.thrs)
				{
					large_f[0] = safe_cast<small_int>(f);
					large_f[1] = 1;
					s = inherit_t::smooth_candidate_e;
				}
				else if (info.have_double)
				{
#ifdef HAVE_DOUBLE_LARGE_PRIME
					if (f > info.thrs3)
						s = inherit_t::smooth_idle_e;
					else if (!custom_prime_test(f))
					{
						int count = 400; // TODO: how to determine this?
						auto p1 = safe_cast<small_int>(pollards_rho(f, count, 2));
						if (p1 > 1)
						{
							auto p2 = safe_cast<small_int>(f / p1);
							if (p1 < p2)
							{
								large_f[0] = p1;
								large_f[1] = p2;
							}
							else
							{
								large_f[0] = p2;
								large_f[1] = p1;
							}
							s = inherit_t::smooth_double_e;
						}
						else
							s = inherit_t::smooth_unfactored_e;
					}
					else
						s = inherit_t::smooth_huge_prime_e;
#endif
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
			bool compose(const smooth_t &rhs, 
				         const large_int &n, 
				         const std::vector<base_ref_t> &base,
						 std::vector<int> *erased = nullptr,
						 std::vector<int> *added = nullptr)
			{
				small_int large[4];
				int index = 0;
				int i = 0;
				int i_rhs = 0;
				large_int sqr1 = sqr;
				while (i < static_cast<int>(s) && i_rhs < static_cast<int>(rhs.s))
					if (large_f[i] < rhs.large_f[i_rhs])
						large[index++] = large_f[i++];
					else if (large_f[i] > rhs.large_f[i_rhs])
						large[index++] = rhs.large_f[i_rhs++];
					else // equal!
					{
						sqr1 = (sqr1 * large_f[i]) % n;
						i++;
						i_rhs++;
					}
				while (i < static_cast<int>(s))
					large[index++] = large_f[i++];
				while (i_rhs < static_cast<int>(rhs.s))
					large[index++] = rhs.large_f[i_rhs++];
#ifdef HAVE_DOUBLE_LARGE_PRIME
				for (int i = index ; i < 2; i++)
					large[i] = 1;
				if (index > 2)
					return false;
#else
				if (index > 0)
					return false;
#endif
				for (int i = 0; i < index; i++)
					large_f[i] = large[i];
				for (int i = index; i < 2; i++)
					large_f[i] = 1;
				s = static_cast<smooth_status_e>(index);
				axb = (axb * rhs.axb) % n;
				sqr = (sqr1 * rhs.sqr) % n;
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
				return true;
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
			bool custom_prime_test(const large_int &n)
			{
				large_int nm1 = n - 1;
				//
				// Begin with a single Fermat test - it excludes a lot of candidates:
				//
				large_int q(228), x, y; // We know n is greater than this, as we've excluded small factors
				x = powm(q, nm1, n);
				if (x != 1u)
					return false;

				q = n - 1;
				unsigned k = lsb(q);
				q >>= k;

				//
				// Execute the trials:
				//
				int seeds[] = { 2, 5, 19, 41 };
				for (auto x : seeds)
				{
					y = powm<large_int>(x, q, n);
					unsigned j = 0;
					while (true)
					{
						if (y == nm1)
							break;
						if (y == 1)
						{
							if (j == 0)
								break;
							return false; // test failed
						}
						if (++j == k)
							return false; // failed
						y = (y * y) % n; // powm<large_int>(y, 2, n);
					}
				}
				return true;  // Yeheh! probably prime.
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
			enum { cached_size_e = 1 << 15 };
			std::vector<real> values;
			std::vector<real> init;

			size_t			  start;
			size_t			  steps;

			real			  offset;
			large_int         thrs2;
			large_int		  thrs3;
			bool			  have_double;
			sieve_stuff_t(large_int largest_prime, bool have_double)
			{
				have_double = have_double;
				thrs2 = largest_prime * largest_prime;
				thrs3 = thrs2 * largest_prime;
			}
			void begin(void) { start = 0; one_step(); }
			bool one_step(void)
			{
				if (start < steps)
				{
#ifdef HAVE_CACHED_SIEVE
					std::copy_n(init.begin() + start, cached_size_e, values.begin());
					start += cached_size_e;
#else
					values = init;
					start += steps;
#endif
					return true;
				}
				else
					return false;
			}
			small_int compute(polynomial_siqs_t<large_int, small_int> &poly, small_int m)
			{
				m = cached_size_e * static_cast<size_t>(m / cached_size_e);
				if (init.empty())
				{
#ifdef HAVE_CACHED_SIEVE
					steps = 2 * static_cast<size_t>(m);
					init = sieve_t(poly, static_cast<size_t>(steps)).fill(offset);
					values.resize(cached_size_e);
#else
					init = sieve_t(poly, static_cast<size_t>(m)).fill(offset);
					values = init ;
#endif
				}
				return m;
			}
		};
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
			k_ = premultiplier(n, primes);
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

			if (m_ < std::sqrt(largest_sieving_prime))
				m_ = largest_sieving_prime * 2;
		}
		large_int process(int order, bool have_double)
		{
			relations_graph_t<smooth_t> candidates(n_);
			smooth_vector_t smooths;
			std::vector<prime_info_t<small_int>> base_info;
			for (auto &item : base_)
				if (item.powers() > 1)
					base_info.push_back(prime_info_t<small_int>{item.prime(0), item.prime(1), item.residue(1)});
			polynomial_generator_t<large_int, small_int> generator(n_, m_, base_info);
			if (order > 0)
				generator.order_init(order);
			LOG_INFO << "Factors of 'a' = " << generator.order() << log_base_t::newline_t();
#ifdef HAVE_TIMING
			time_estimator_t time_estimator(base_.size());
#endif

#ifdef HAVE_THREADING
			int cores = system_info_t::cores();
			for (int i = 0; i < cores; i++)
			{
				threads_.emplace_back(&self_initializing_quadratic_sieve_t::sieving_thread, this, std::ref(base_info), have_double);
				polynomials_.push(generator());
				polynomials_.push(generator());
			}
#else
			polynomial_siqs_t<large_int, small_int> poly(generator().index, base_info, n_);

			sieve_stuff_t sieve_stuff(base_.rbegin()->prime(0), have_double);
			m_ = sieve_stuff.compute(poly, m_);
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
				process_candidates_chunk(base_, candidates, chunk, smooths, n_);
			}
#endif
			LOG_INFO << log_base_t::newline_t()
#ifdef HAVE_CANDIDATE_ANALYSYS
					 << "Attempted " << analysis_.smooth_attempts << "\n"
#ifdef HAVE_DOUBLE_LARGE_PRIME
					  << "Detail: "  << analysis_.smooth_huge_prime
					  << " (" << analysis_.smooth_huge_prime * 100 / analysis_.smooth_attempts << "%) likely prime "
					  << ", "  << analysis_.smooth_large_composite
					  << " (" << analysis_.smooth_large_composite * 100 / analysis_.smooth_attempts << "%) composite "
					  << " and " << analysis_.smooth_unfactored
					  << " (" << analysis_.smooth_unfactored * 100 / analysis_.smooth_attempts << "%) unfactored "
				  	  << " left out: " << analysis_.smooth_idle
#endif
#endif // HAVE_DOUBLE_LARGE_PRIME
				     << "\n";
			if (smooths.size() < actual_bsize)
			{
				LOG_ERROR << "Found only " << smooths.size() << log_base_t::newline_t();
				return 1;
			}
#ifdef HAVE_CANDIDATE_ANALYSYS1
			inherit_t::print_analysis(base_.rbegin()->prime(0));
#endif
#ifdef _DEBUG
			for (auto s : poly_stat_)
				LOG_INFO << "\t" << s / static_cast<double>(poly_count_) << log_base_t::newline_t();
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
			smooth_info_t info = { poly, sieve_stuff.thrs2, sieve_stuff.thrs3, n_, sieve_stuff.have_double};
			info.runs.reserve(base_.size() * 2);
			size_t count = info.poly.count();
#ifdef _DEBUG
			if (poly_stat_.empty())
				poly_stat_.resize(count);
			poly_count_++;
#endif
			for (size_t c = 0; c < count; c++)
			{
				info.poly.select(c);
				if (c == 0)
					sieve_t::build_run(info.poly, base_, info.runs);
				else
					sieve_t::update_run(info.poly, info.runs);

				sieve_stuff.begin(); 
				do 
				{
					size_t size = sieve_stuff.values.size();
					for (auto &run : info.runs)
					{
						if (run.p < 12)
							continue;
#ifdef HAVE_CACHED_SIEVE
						small_int m1 = m_ - sieve_stuff.start + sieve_stuff_t::cached_size_e + run.p * m_;
						int index0 = static_cast<int>((m1 + run.x[0]) % run.p);
						int index1 = static_cast<int>((m1 + run.x[1]) % run.p);
#else
						int index0 = static_cast<int>(m_ + run.x[0]) % run.p;
						int index1 = static_cast<int>(m_ + run.x[1]) % run.p;
#endif
						int index = std::min(index0, index1);
						int delta = std::max(index0, index1) - index;
						int size1 = static_cast<int>(size) - delta;
						int i = index;
						for (; i < size1; i += static_cast<int>(run.p))
						{
							sieve_stuff.values[i] -= run.lg;
							sieve_stuff.values[i + delta] -= run.lg;
						}
						if (i < static_cast<int>(size))
							sieve_stuff.values[i] -= run.lg;
					}
#ifdef _DEBUG
					size_t prev = result.size();
#endif
					collect_smooth(info, sieve_stuff, result);
#ifdef _DEBUG
					poly_stat_[c] += result.size() - prev;
#endif
				} while (sieve_stuff.one_step());
			}
			return result;
		}
		void collect_smooth(const smooth_info_t &info,
							const sieve_stuff_t &sieve,
							smooth_vector_t &result)
		{
			real sieve_thrs = 2 * base_.rbegin()->logp_ + sieve.offset + 1 * real_op_t<real>::unit(); // small prime variation
			if (info.have_double)
				sieve_thrs += base_.rbegin()->logp_ - 2 * real_op_t<real>::unit();

			size_t size = sieve.values.size();
			for (size_t i = 0; i < size; i++)
				if (sieve.values[i] < sieve_thrs)
				{
#ifdef HAVE_CACHED_SIEVE
					smooth_t s(info, i - m_ + sieve.start - sieve_stuff_t::cached_size_e);
#else
					smooth_t s(info, i - m_);
#endif
					analysis_.smooth_attempts++;
					
					if (s.type() <= inherit_t::smooth_double_e)
					{
#if DBG_SIEVE >= DBG_SIEVE_INFO
						if (!s.invariant(n_, base_))
							std::cout << "Hm\n";
#endif // DBG_SIEVE	
#ifdef HAVE_DOUBLE_LARGE_PRIME
						if (s.type() == smooth_double_e)
							analysis_.smooth_large_composite++;
#endif
						result.push_back(s);
					}
					else
					{
#ifdef HAVE_DOUBLE_LARGE_PRIME
						switch (s.type())
						{
						case smooth_unfactored_e:
							analysis_.smooth_unfactored++;
							break;
						case smooth_huge_prime_e:
							analysis_.smooth_huge_prime++;
							break;
						default:
							analysis_.smooth_idle++;
							break;
						}
#endif
					}
				}
		}
#ifdef HAVE_THREADING
		void sieving_thread(const std::vector<prime_info_t<small_int>> &base_info, bool have_double)
		{
			sieve_stuff_t sieve_stuff(base_.rbegin()->prime(0), have_double);
			try
			{
				for (auto seed = polynomials_.pop(); !seed.is_null(); seed = polynomials_.pop())
				{
					polynomial_siqs_t<large_int, small_int> poly(seed.index, base_info, n_);
					sieve_stuff.compute(poly, m_);
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
		int  process_candidates_chunk(std::vector<base_ref_t> &base,
 									  relations_graph_t<smooth_t>  &relations,
									  std::vector<smooth_t> &chunk,
									  std::vector<smooth_t> &smooths,
									  const large_int &n)
		{
			int promoted = 0;
			for (auto &smooth : chunk)
			{
				switch (smooth.type())
				{
				case smooth_valid_e:
					register_smooth(smooths, smooth);
					break;
				case smooth_candidate_e:
				case smooth_double_e:
					if (relations.add_smooth(smooth, base))
					{
						register_smooth(smooths, smooth);
						promoted++;
					}
					break;
				default: //, smooth_idle_e, smooth_unfactored_e
					break;
				}
			}
			return promoted;
		}
		void register_smooth(std::vector<smooth_t> &smooths, const smooth_t &smooth)
		{
			int si = static_cast<int>(smooths.size());
			for (auto f : smooth.factors())
				base_[f].smooths.push_back(si);
			smooths.push_back(smooth);
		}
		large_int					smooth_thrs_; // square of last element of base
		large_int					n_; // number to factor
		small_int					m_; // size of sieving interval
#ifdef HAVE_MULTIPLIER
		small_int					k_;
#endif
		std::vector<base_ref_t>		base_;
#ifdef _DEBUG
		std::vector<size_t> poly_stat_;
		int					poly_count_ = 0;
#endif
#ifdef HAVE_THREADING
		shared_list_t<polynomial_seed_t> polynomials_;
		shared_list_t<smooth_vector_t>   smooths_found_;
		std::vector<std::thread>	     threads_;
#endif
#ifdef HAVE_CANDIDATE_ANALYSYS
		analysis_t analysis_;
#endif
	};


	template <class large_int, class small_int = int, class real = float>
	large_int self_initializing_quadratic_sieve(const large_int &n, const large_int &m, small_int base_size, int order = 0, bool have_double = false)
	{
#if DBG_SIEVE >= DBG_SIEVE_INFO
		std::cout << "Factorization of " << n << std::endl;
#endif
		self_initializing_quadratic_sieve_t<large_int, small_int, real> qs(n, m, base_size);
		return qs.process(order, have_double);
	}

};

#endif
