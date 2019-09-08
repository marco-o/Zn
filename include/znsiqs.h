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
#include <atomic>
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

	template <class large_int>
	class trial_division_t
	{
	public:
		typedef unsigned int value_type;
		typedef long long temp_type;
		trial_division_t(const large_int &n)
		{
			boost::multiprecision::export_bits(n, std::back_inserter(value_), 32, false);
			temp_ = value_;
			order_ = static_cast<int>(value_.size()) - 1;
		}
		bool complete(void) const
		{
			return value_.size() == 1 && value_[0] == 1;
		}
		temp_type reminder(void) const
		{
			switch (value_.size())
			{
			case 1:
				return value_[0] ;
			case 2:
				return ((static_cast<temp_type>(value_[1]) << 32) + value_[0]) ;
			default:
				break;
			}
			return std::numeric_limits<temp_type>::max();
		}
		bool divide(int n)
		{
			temp_type r = 0;
			for (int i = order_; i >= 0; i--)
			{
				r = (r << 32) + value_[i];
				temp_[i] = static_cast<value_type>(r / n);
				r = r % n;
			}
			if (r != 0)
				return false;
			if (temp_[order_] == 0)
			{
				temp_.pop_back();
				order_--;
			}
			value_ = temp_;
			return true;
		}
	private:
		std::vector<value_type> value_;
		std::vector<value_type> temp_;
		int order_;
	};
	struct analysis_t
	{
		std::atomic<int> smooth_attempts = 0;
		std::atomic<int> smooth_prime_unused = 0;
		std::atomic<int> smooth_idle = 0;
		std::atomic<int> promoted = 0;
		std::atomic<int> partial_promoted = 0;
		std::atomic<int> direct = 0;
		std::atomic<int> smooth_large_composite = 0;
		std::atomic<int> smooth_unfactored = 0;
		std::atomic<int> smooth_huge_prime = 0;
		std::atomic<int> smooth_huge_composite = 0;
	};

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

	struct sieving_options_t
	{
		int order = 0; // number of factors osf 'a'
		int base_size = 0; // number of primes in the factoring base
		int multiplier = 0; // this forces the search for one, 1 avoids it
		bool have_double = false; // use partial-partial relations
		int m = 0 ; // sieving interval (half of it)
		int sieve_bias = 0; // to experiment with higher or lower threshold
#ifdef _DEBUG
		bool have_threading = false;
#else
		bool have_threading = true;
#endif
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
			small_int	thrs2;
			small_int	thrs20;
			small_int	thrs3;
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
				trial_division_t<large_int> f1(f);
				size_t runs_size = info.runs.size();
				for (size_t j = 0; j < runs_size; j++)
				{
					const auto &run = info.runs[j];
					bool plus = ((x - run.x[0]) % run.p == 0);
					bool minus = ((x - run.x[1]) % run.p == 0);
					//if (run.a != 0)
					if (!plus && !minus)
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
					int power = 0;
#if 1
					while (f1.divide(run.p))
						if (++power % 2 == 0)
							sqr = (sqr * run.p) % info.n;
#else
					large_int p = run.p;
					while (divide_qr(f, p))
						if (++power % 2 == 0)
							sqr = (sqr * p) % info.n;
#endif
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
				small_int reminder = f1.reminder();
				if (reminder == 1)
					s = inherit_t::smooth_valid_e;
				else if (reminder < info.thrs2)
				{
					large_f[0] = reminder;
					large_f[1] = 1;
					s = inherit_t::smooth_candidate_e;
				}
				else if (reminder < info.thrs20)
					s = inherit_t::smooth_prime_unused_e;
				else if (info.have_double)
				{
					if (reminder > info.thrs3)
						s = inherit_t::smooth_idle_e;
					else
					{
#ifdef HAVE_FACTORIZATION_TEST
						auto &file = test_file();
						file << reminder;
#endif
						if (!custom_prime_test(reminder))
						{
							int count = 800; // TODO: how to determine this?
#ifdef HAVE_FACTORIZATION_TEST
							file << " c ";
#endif
							large_f[0] = 1; // safe_cast<small_int>(pollards_rho(reminder, count, 2));
							if (large_f[0] > 1)
							{
								large_f[1] = safe_cast<small_int>(reminder / large_f[0]);
								if (large_f[0] > large_f[1])
									std::swap(large_f[0], large_f[1]);
#ifdef HAVE_FACTORIZATION_TEST
								test_file() << large_f[0] << " x " << large_f[1];
#endif
								if (large_f[1] > info.thrs2)
									s = inherit_t::smooth_huge_composite_e;
								else
									s = inherit_t::smooth_double_e;
							}
							else
							{
								s = inherit_t::smooth_unfactored_e;
							}
						}
#ifdef HAVE_FACTORIZATION_TEST
						else
							file << " p";
						file << "\n";
#endif
					}
				}
				else
					s = inherit_t::smooth_huge_prime_e;
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
				for (int i = index ; i < 2; i++)
					large[i] = 1;
				if (index > 2)
					return false;
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
#ifdef HAVE_FACTORIZATION_TEST
			std::ofstream &test_file(void)
			{
				static std::ofstream ofile("test.txt");
				return ofile;
			}
#endif
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

			size_t			  start;
			size_t			  steps;

			real			  offset;
			small_int         thrs2;  // reduced square; values bugger are not considered
			small_int		  thrs20; // exact square
			small_int		  thrs3;
			bool			  have_double;
			sieve_stuff_t(small_int largest_prime, bool have_dbl)
			{
				have_double = have_dbl;
				thrs20 = largest_prime * largest_prime;
				const int div = 8;
				thrs2 = thrs20 / div;
				thrs3 = thrs2 * largest_prime / div;
			}
			void compute(polynomial_siqs_t<large_int, small_int> &poly, small_int m)
			{
				if (init.empty())
					init = sieve_t(poly, static_cast<size_t>(m)).fill(offset);
			}
		};
		typedef std::vector<smooth_t> smooth_vector_t;
		self_initializing_quadratic_sieve_t(const large_int &n, const sieving_options_t &opt) : n_(n), options_(opt)
		{
			// double the range; half of them won't be a quadratic residue
#ifdef HAVE_TIMING
			start_ = std::chrono::steady_clock::now();
#endif
			small_int range;
			if (options_.base_size != 0)
				range = inherit_t::primes_range(options_.base_size * 2);
			else
			{
				double n1 = safe_cast<double>(n);
				double e1 = std::log(n1);
				double e2 = std::sqrt(e1 * std::log(e1)) / 2;
				range = static_cast<small_int>(std::exp(e2));
			}
			auto primes = eratosthenes_sieve<small_int>(static_cast<int>(range));
			if (options_.multiplier == 0)
				options_.multiplier = premultiplier(n, primes);
			if (options_.multiplier != 1)
			{
				LOG_INFO << "Premultiplier = " << options_.multiplier << log_base_t::newline_t();
				n_ *= options_.multiplier;
			}
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

			if (options_.m < std::sqrt(largest_sieving_prime))
				options_.m = static_cast<int>(largest_sieving_prime * 2);
		}
		large_int process(void)
		{
			relations_graph_t<smooth_t> candidates(n_);
			smooth_vector_t smooths;
			std::vector<prime_info_t<small_int>> base_info;
			for (auto &item : base_)
				if (item.powers() > 1)
					base_info.push_back(prime_info_t<small_int>{item.prime(0), item.prime(1), item.residue(1)});
			polynomial_generator_t<large_int, small_int> generator(n_, options_.m, base_info);
			if (options_.order > 0)
				generator.order_init(options_.order);
			LOG_INFO << "Factors of 'a' = " << generator.order() << log_base_t::newline_t();
#ifdef HAVE_TIMING
			time_estimator_t time_estimator(base_.size());
			log_time("Startup");
#endif
			sieve_stuff_t sieve_stuff(base_.rbegin()->prime(0), options_.have_double);
			int cores = system_info_t::cores();
			if (options_.have_threading)
				for (int i = 0; i < cores; i++)
				{
					threads_.emplace_back(&self_initializing_quadratic_sieve_t::sieving_thread, this, std::ref(base_info));
					polynomials_.push(generator());
					polynomials_.push(generator());
				}
			else
			{
				polynomial_siqs_t<large_int, small_int> poly(generator().index, base_info, n_);
				sieve_stuff.compute(poly, options_.m);
			}

			int count = 0;
			size_t actual_bsize = 0;
			while (smooths.size() < actual_bsize + 5 + actual_bsize / 100)
			{
				std::vector<smooth_t> chunk;
				if (options_.have_threading)
				{
					chunk = smooths_found_.pop();
					polynomials_.push(generator());
				}
				else
				{
					auto seed = generator();
					polynomial_siqs_t<large_int, small_int> poly(seed.index, base_info, n_);
					chunk = sieve(poly, sieve_stuff);
				}
				process_candidates_chunk(base_, candidates, chunk, smooths, n_);
				count++;
				actual_bsize = std::max(inherit_t::actual_base_size(base_), base_info.size() / 3) ;
#ifdef HAVE_TIMING
				time_estimator.update(static_cast<int>(smooths.size()), analysis_.promoted);
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
			if (options_.have_threading)
			{
				polynomials_.clear();
				for (int i = 0; i < cores; i++)
					polynomials_.push(polynomial_seed_t());
				for (auto &thread : threads_)
					thread.join();
				// pick remaining smooths
				while (!smooths_found_.empty())
				{
					auto chunk = smooths_found_.pop();
					process_candidates_chunk(base_, candidates, chunk, smooths, n_);
				}
			}
			LOG_INFO << log_base_t::newline_t()
			         << "Attempted  " << analysis_.smooth_attempts << "\n"
				     << "Direct     " << analysis_.direct << "\n"
					 << "Promoted   " << analysis_.promoted << "\n"
					 << "Large p.   " << analysis_.smooth_prime_unused << "\n"
					 << "Idle       " << analysis_.smooth_idle << "\n" ;
			if (options_.have_double)
				LOG_INFO
				     << "C.Partial  " << analysis_.smooth_large_composite << "\n"
					 << "P.Partial  " << analysis_.partial_promoted << "\n"
					 << "Huge p.    " << analysis_.smooth_huge_prime << "\n"
					 << "Huge c.    " << analysis_.smooth_huge_composite << "\n"
					 << "Unfactored " << analysis_.smooth_unfactored << "\n"
 				     << "\n";
			log_time("Sieving");
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

			large_int n1 = n_;
			if (options_.multiplier != 1)
				n1 /= options_.multiplier ;

			for (auto &item : basemix)
			{
				smooth_t s;
				for (auto index : item)
					s.compose(smooths[index], n_, base_);
#if DBG_SIEVE >= DBG_SIEVE_DEBUG
				if (!s.square())
					std::cout << "Non null factors!\n";
				if (!s.invariant(n_, base_))
					std::cout << "Hmmmm";
#endif // DBG_SIEVE	
				r = s.result(n1);
				if (r != 1 && r != n1)
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
			log_time("Final");
			return r;
		}
		smooth_vector_t sieve(const polynomial_siqs_t<large_int, small_int> &poly,
							  sieve_stuff_t &sieve_stuff)
		{
			// build vector for sieving
			smooth_vector_t result;
			smooth_info_t info = { poly, sieve_stuff.thrs2, sieve_stuff.thrs20, sieve_stuff.thrs3, n_, sieve_stuff.have_double};
			info.runs.reserve(base_.size() * 2);
			size_t count = info.poly.count();
#ifdef _DEBUG
			if (poly_stat_.empty())
				poly_stat_.resize(count);
			poly_count_++;
#endif
			std::vector<typename sieve_t::sieve_run_t>::iterator begin, end;
			for (size_t c = 0; c < count; c++)
			{
				info.poly.select(c);
				if (c == 0)
				{
					sieve_t::build_run(info.poly, base_, info.runs);
					end = info.runs.end();
					begin = info.runs.begin();
					for (; begin != end; ++begin)
						if (begin->p > 12)
							break;
				}
				else
					sieve_t::update_run(info.poly, info.runs);

				sieve_stuff.values = sieve_stuff.init;
				size_t size = sieve_stuff.values.size();
				auto it = begin;
#if 1
				const int loop_unroll = 4;
				for (; it != end; ++it)
				{
					auto &run = *it;
					int index0 = static_cast<int>((options_.m + run.x[0]) % run.p);
					int index1 = static_cast<int>((options_.m + run.x[1]) % run.p);
					int indexmin = std::min(index0, index1);
					int indexmax = std::max(index0, index1);
					int loops = static_cast<int>((size - indexmax) / run.p) - loop_unroll + 1;
					if (loops < loop_unroll * 2)
						break;
					int delta0 = indexmax - indexmin;
					int delta1 = static_cast<int>(run.p) - delta0;
					real *v = &sieve_stuff.values[0] + indexmin;
					for (int i = 0; i < loops; i += 4)
					{
						v[0] -= run.lg; v += delta0;
						v[0] -= run.lg;	v += delta1;
						v[0] -= run.lg; v += delta0;
						v[0] -= run.lg;	v += delta1;
						v[0] -= run.lg; v += delta0;
						v[0] -= run.lg;	v += delta1;
						v[0] -= run.lg; v += delta0;
						v[0] -= run.lg;	v += delta1;
					}
					int i = static_cast<int>(v - &sieve_stuff.values[0]);
					int size1 = static_cast<int>(size) - delta0;
					for (; i < size1; i += static_cast<int>(run.p))
					{
						sieve_stuff.values[i] -= run.lg;
						sieve_stuff.values[i + delta0] -= run.lg;
					}
					if (i < static_cast<int>(size))
						sieve_stuff.values[i] -= run.lg;
				}
#endif
				for (; it != end; ++it)
				{
					auto &run = *it;
					small_int m1 = options_.m + run.p;
					int index0 = static_cast<int>((m1 + run.x[0]) % run.p);
					int index1 = static_cast<int>((m1 + run.x[1]) % run.p);
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
			}
			return result;
		}
		void collect_smooth(const smooth_info_t &info,
							const sieve_stuff_t &sieve,
							smooth_vector_t &result)
		{
			real lg = base_.rbegin()->logp_;
			real un = real_op_t<real>::unit();
			real sieve_thrs = 2 * base_.rbegin()->logp_ + sieve.offset + 4 * real_op_t<real>::unit(); // small prime variation
			if (info.have_double)
				sieve_thrs += base_.rbegin()->logp_ - 2 * real_op_t<real>::unit();
			if (options_.sieve_bias)
				sieve_thrs = static_cast<real>(sieve_thrs + options_.sieve_bias);

			size_t size = sieve.values.size();
			for (size_t i = 0; i < size; i++)
				if (sieve.values[i] < sieve_thrs)
				{   // static cast required, otherwise operation is done on unsigend
					smooth_t s(info, static_cast<int>(i) - options_.m); 
					analysis_.smooth_attempts++;
					if (s.type() <= inherit_t::smooth_double_e)
					{
#if DBG_SIEVE >= DBG_SIEVE_DEBUG
						if (!s.invariant(n_, base_))
							std::cout << "Hm\n";
#endif // DBG_SIEVE	
						if (s.type() == smooth_double_e)
							analysis_.smooth_large_composite++;
						result.push_back(s);
					}
					else
					{
						switch (s.type())
						{
						case smooth_prime_unused_e:
							analysis_.smooth_prime_unused++;
							break;
						case smooth_unfactored_e:
							analysis_.smooth_unfactored++;
							break;
						case smooth_huge_prime_e:
							analysis_.smooth_huge_prime++;
							break;
						case smooth_huge_composite_e:
							analysis_.smooth_huge_composite++;
							break;
						default:
							analysis_.smooth_idle++;
							break;
						}
					}
				}
		}
		void sieving_thread(const std::vector<prime_info_t<small_int>> &base_info)
		{
			sieve_stuff_t sieve_stuff(base_.rbegin()->prime(0), options_.have_double);
			try
			{
				for (auto seed = polynomials_.pop(); !seed.is_null(); seed = polynomials_.pop())
				{
					polynomial_siqs_t<large_int, small_int> poly(seed.index, base_info, n_);
					sieve_stuff.compute(poly, options_.m);
					auto chunk = sieve(poly, sieve_stuff);
					smooths_found_.push(chunk);
				}
			}
			catch (std::exception &exc)
			{
				std::cout << "Exception in thread: " << exc.what() << std::endl;
			}
		}
		void process_candidates_chunk(std::vector<base_ref_t> &base,
 									  relations_graph_t<smooth_t>  &relations,
									  std::vector<smooth_t> &chunk,
									  std::vector<smooth_t> &smooths,
									  const large_int &n)
		{
			for (auto &smooth : chunk)
			{
				auto s = smooth.type();
				switch (s)
				{
				case smooth_valid_e:
					register_smooth(smooths, smooth);
					analysis_.direct++;
					break;
				case smooth_candidate_e:
				case smooth_double_e:
					if (relations.add_smooth(smooth, base))
					{
						register_smooth(smooths, smooth);
						if (s == smooth_candidate_e)
							analysis_.promoted++;
						else
							analysis_.partial_promoted++;
					}
					break;
				default: //, smooth_idle_e, smooth_unfactored_e
					break;
				}
			}
		}
		void register_smooth(std::vector<smooth_t> &smooths, const smooth_t &smooth)
		{
			int si = static_cast<int>(smooths.size());
			for (auto f : smooth.factors())
				base_[f].smooths.push_back(si);
			smooths.push_back(smooth);
		}
#ifdef HAVE_TIMING
		void log_time(const char *msg)
		{
			auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_).count();
			LOG_INFO << msg << " time " << elapsed / 1000.0<< log_base_t::newline_t();
		}
#else
		void log_time(const char *) {}
#endif
		large_int					n_; // number to factor
		std::vector<base_ref_t>		base_;
#ifdef _DEBUG
		std::vector<size_t> poly_stat_;
		int					poly_count_ = 0;
#endif

		shared_list_t<polynomial_seed_t> polynomials_;
		shared_list_t<smooth_vector_t>   smooths_found_;
		std::vector<std::thread>	     threads_;

		analysis_t analysis_;
#ifdef HAVE_TIMING
		std::chrono::steady_clock::time_point start_;
#endif
		sieving_options_t options_;
	};



	template <class large_int, class small_int = int, class real = float>
	large_int self_initializing_quadratic_sieve(const large_int &n, const sieving_options_t &opt)
	{
#if DBG_SIEVE >= DBG_SIEVE_INFO
		std::cout << "Factorization of " << n << std::endl;
#endif
		self_initializing_quadratic_sieve_t<large_int, small_int, real> qs(n, opt);
		return qs.process();
	}

};

#endif
