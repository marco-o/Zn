/*
	Test program to develop small number factoriation
	(for partial-partial relations)
*/

#include "znbasic.h"
#include "zneratosthenes_sieve.h"
#include "znelliptic_curve_fact.h"
#include "znquadratic_residue.h"
#include "boost/multiprecision/cpp_int.hpp"
#include <vector>
#include <fstream>
#include <string>
#include <numeric>
#include <map>

#define HAVE_LARGE_PRIME
//#define HAVE_MULTIPLIER_STATS

namespace zn
{
	class stats_t
	{
	public:
		void processed(void) { processed_++; }
		void direct(void) { direct_++; }
		void discarded(void) { discarded_++; }
		void promoted(void) { promoted_++; }
		void polynomials(void) { polynomials_++; }
		void large_primes(void) { large_primes_++; }
		void resume(void)
		{
			std::cout << "Processed = " << processed_ << "\n"
				      << "Discarded = " << discarded_ << "\n"
#ifdef HAVE_LARGE_PRIME
					  << "Direct = " << direct_ << "\n"
					  << "Large = " << large_primes_ << "\n"
					  << "Promoted = " << promoted_ << "\n"
#endif
					  << "Accepted = " << (processed_ - discarded_) << "\n"
					  << "Polynomials = " << polynomials_  << "\n"
				;
		}
	private:
		int processed_ = 0;
		int promoted_ = 0;
		int discarded_ = 0; // smooth detection failures
		int large_primes_ = 0;
		int direct_ = 0;
		int polynomials_ = 0;
	};
	struct stats_none_t
	{
	public:
		void processed(void) {}
		void direct(void) {}
		void discarded(void) {}
		void promoted(void) {}
		void polynomials(void) {}
		void large_primes(void) {}
		void resume(void) {}
	};
	struct config_qs_t
	{
		int		m2 = 1000 ;
		size_t	base_size = 32;
		double	sieve_offset = 0;
		double	pquality = 5.0;
		bool	more_polys = false;
	};
	// A quadratic sieve class for factoring numbers up to (about) 19 (2^63) digits
	template <class stat_t = stats_t, class large_int = long long, class small_int = int>
	class quadratic_sieve_cached_t
	{
	public:
		typedef unsigned char real_t;
		enum { log_unit_e = 8 };
		struct prime_t
		{
			small_int value;
			real_t	  logp;
			float	  logf;
			std::vector<small_int> inverse;
			std::vector<small_int> residue; // quadratic residues
			prime_t(small_int v = 2) : value(v), residue(v, 0), inverse(1, 0)
			{
				logf = static_cast<float>(std::log(v));
				logp = static_cast<real_t>(logf * log_unit_e);
				std::map<small_int, small_int> squares;
				for (small_int i = 1; i < v; i++)
				{
					inverse.push_back(std::get<2>(extended_euclidean_algorithm<small_int>(v, i)));
					small_int sq = i * i % v;
					squares.insert(std::make_pair(sq, i));
				}
				for (auto it : squares)
					residue[it.first] = it.second;
			}
		};
		class multiplier_t
		{
		public:
			multiplier_t(int km) : logk_(2, 0.0f)
			{
				for (int i = 2; i < km; i++)
					logk_.push_back(static_cast<float>(log(i)));
				for (int i = 2; i * i < km; i++)
					logk_[i*i] = 0;
			}
			template <class PrimeIt>
			std::pair<int, float> find_best(const large_int &n, const std::pair<PrimeIt, PrimeIt> &range)
			{
				float q;
				std::pair<int, float> result(3, quality(n, 3, range));
				int ks = static_cast<int>(logk_.size());
				for (int k = 5; k < ks; k += 2)
					if ((logk_[k] > 0) && ((q = quality(n, k, range))) > result.second)
						result = std::make_pair(k, q);
				return result;
			}
			template <class PrimeIt>
			float quality(const large_int &n, int k, const std::pair<PrimeIt, PrimeIt> &range)
			{
				large_int kn = k * n;
				float fk = -logk_[k] / 2;
				for (auto iter = range.first; iter != range.second ; ++iter)
					if (iter->value == 2)
					{
						switch (safe_cast<int>(kn % 8))
						{
						case 1:
							fk += 2 * iter->logp;
							break;
						case 2:
						case 3:
						case 6:
						case 7:
							fk += iter->logp / 2;
							break;
						case 5:
							fk += iter->logp;
							break;
						default:
							break;
						}
					}
					else
					{
						small_int knp = static_cast<small_int>(kn % iter->value);
						if (iter->residue[knp])
							if (k % iter->value == 0)
								fk += iter->logf / iter->value;
							else
								fk += 2 * iter->logf / (iter->value - 1);
					}
				return fk;
			}
		private:
			std::vector<float> logk_;
		};
		typedef long long factors_t;
		struct smooth_t
		{
			large_int sqr; // product of prime with even exponents (/ 2)
			large_int axb;// a *x + b, the number squared
			smooth_t(large_int ab = 1) : axb(ab), sqr(1) {}
		};
		struct poly_seed_t
		{
			double quality;
			int a1;
			int a2;
		};
		struct poly_t
		{
			small_int a;
			small_int a2;
			small_int b;
			large_int c;
			small_int m; // suggested value for sieve interval
			small_int m2; // m * 2
		};
		struct info_t
		{
			large_int					 n;
			int							 multi; // multiplier, most times 1
			std::vector<const prime_t *> base; 
			std::vector<real_t>			 values;
			std::vector<smooth_t>		 smooths;
			std::vector<factors_t>		 factors;
			int							 factors_index; // index of last smooth into factors
														// after first attempt is no longer smooths.size()
														// but restarts from 0
#ifdef HAVE_LARGE_PRIME
			std::map<small_int, smooth_t> large_primes;
#endif
			std::vector<poly_seed_t>	 poly;
			small_int					 sqr2n; // sqrt(n * 2)
			small_int					 sqrn; // sqrt(n)
			config_qs_t					 config;
		};
		quadratic_sieve_cached_t(small_int primes = 256) : multiplier_(12)
		{
			small_int bound = static_cast<small_int>(primes * std::pow(std::log(primes), 1.3));
			std::vector<small_int> prime = eratosthenes_sieve<small_int>(bound);
			if (prime.size() > static_cast<size_t>(primes))
				prime.erase(prime.begin() + primes, prime.end());
			for (auto p : prime)
				primes_.push_back(prime_t(p));
			primes_range_ = std::make_pair(primes_.begin(), primes_.begin() + std::min<int>(10, primes_.size()));
#ifdef HAVE_MULTIPLIER_STATS
			mstats_.open("mstats.txt");
#endif
		}
		~quadratic_sieve_cached_t(void)
		{
			stats_.resume();
		}
		small_int factor(const large_int &n, info_t &info)
		{
			small_int result = factor_imp(n, info, false);
			if (result == 1)
				result = factor_imp(n, info, true);
			if (result == 1)
				std::cout << n << ": found " << info.smooths.size() << ", required = " << info.base.size() << std::endl;
			return result;
		}
		small_int factor_imp(const large_int &n, info_t &info, bool with_multi = false)
		{
			small_int result = 1;
			init_info(n, info, with_multi);
			std::sort(info.poly.begin(), info.poly.end(), [](const poly_seed_t &lhs, const poly_seed_t &rhs) {
				return lhs.quality < rhs.quality;
			});
#ifdef HAVE_MULTIPLIER_STATS
			int poly_count = 0;
#endif
			for (const auto &poly_seed : info.poly)
			{
				stats_.polynomials();
				poly_t poly = create_poly(info, poly_seed);
				info.values.resize(poly.m * 2);
				std::fill(info.values.begin(), info.values.end(), 0);
				sieve(poly, info);
				result = collect_smooth(info, poly);
#ifdef HAVE_MULTIPLIER_STATS
				poly_count++;
#endif
				if (result != 1)
					break;
			}
#ifdef HAVE_MULTIPLIER_STATS
			mstats_ << multiplier_.quality(n, 1, primes_range_) << " " << poly_count << "\n" ;
#endif
			return result;
		}
	private:
		poly_t create_poly(info_t &info, const poly_seed_t &seed)
		{
			poly_t poly;

			if (seed.a1 > 0)
			{
				poly.a = seed.a1;
				poly.a2 = poly.a * poly.a;
				poly.b = quadratic_residue<small_int>(static_cast<small_int>(info.n % poly.a2), poly.a2, poly.a); // then take 'square root'...
				if (seed.a2 > 1)
				{
					small_int p1 = seed.a2;
					poly.a *= seed.a2;
					small_int g = seed.a2 * seed.a2;
					small_int r1 = quadratic_residue<small_int>(static_cast<small_int>(info.n % g), g, seed.a2); 

					small_int db = (r1 - poly.b);
					auto ext = extended_euclidean_algorithm(poly.a2, g);
					small_int h = (std::get<1>(ext) * db) % g;
					poly.b = poly.b + h * poly.a2;
					poly.a2 *= g;
					if (poly.b > poly.a2 / 2)
						poly.b = poly.a2 - poly.b;
				}
				poly.c = (poly.b * poly.b - info.n) / poly.a2;
				poly.m = 2 * (info.sqr2n / (2 * poly.a2)); // I want it even
				poly.m = std::max(std::min(poly.m, info.config.m2), 2 * (info.config.m2 / 4));
			}
			else
			{
				poly.a = 1;
				poly.a2 = 1 ;
				poly.b = 0;
				poly.c = info.n;
				poly.m = info.config.m2 * 4;
			}
			poly.m2 = poly.m * 2;
			return poly;
		}
		void init_info(large_int n, info_t &info, bool with_multi)
		{
			info.n = n ;
			if (with_multi)
			{
				auto multi = multiplier_.find_best(n, primes_range_);
				info.n = n * multi.first;
				info.multi = multi.first;
			}
			else
			{
				info.n = n;
				info.multi = 1;
			}
			info.base.clear();
			info.smooths.clear();
			info.poly.clear();
			info.factors.resize(info.config.base_size + 1);
			std::fill(info.factors.begin(), info.factors.end(), 0);
#ifdef HAVE_LARGE_PRIME
			info.large_primes.clear();
#endif
			info.sqrn  = static_cast<small_int>(std::sqrt(info.n));
			info.sqr2n = static_cast<small_int>(std::sqrt(info.n * 2));
			small_int m = info.config.m2;
			double a0 = sqrt(info.sqr2n / m);
			auto it = primes_.begin();
			auto end = primes_.end();
			double k = 0.85;
			for (; it != end; ++it)
				if (it->residue[info.n % it->value]) // keep it!
					if (info.base.size() < info.config.base_size)
					{
						info.base.push_back(&(*it));
						double t = k * it->value / a0;
						double q = t + 1 / t;
						if (q < info.config.pquality)
						{
							poly_seed_t poly{ q, it->value, 1};
							info.poly.push_back(poly);
						}
					}
					else
						break;

			if (info.config.more_polys)
			{
				size_t s = info.base.size() / 3;
				for (size_t i = 1; i < s; i++)
				{
					small_int p = info.base[i]->value;
					for (size_t j = 0; j < i; j++)
					{
						small_int q = info.base[j]->value;
						auto pq = q * p;
						double t = k * pq / a0;
						if (t > info.config.pquality)
							break;
						auto q1 = t + 1 / t + 1.0; // add a small penalty..
						if (q1 < info.config.pquality)
						{
							poly_seed_t poly{ q1, p, q };
							info.poly.push_back(poly);
						}
					}
				}
			}

			poly_seed_t poly{2 * info.config.pquality + 1, -1 }; // worse quality
			info.poly.push_back(poly);
		}
		void sieve(const poly_t &poly, info_t &info)
		{
			if (poly.a == 1)
			{
				for (auto p : info.base)
				{
					small_int t = p->residue[info.n % p->value];
					small_int r = 2 * p->value - ((info.sqrn - poly.m) % p->value) ;
					small_int i0 = (r + t) % p->value;
					if (p->value > 3)
				/*	{
						for (int i = i0; i < poly.m2 ; i += 2)
							info.values[i] += p->logp;
					}
					else */
					{
						small_int i1 = (r - t) % p->value;
						small_int index = std::min(i0, i1);
						small_int delta = std::max(i0, i1) - index;
						small_int m2 = poly.m2 - delta;
						small_int i = index;
						for (; i < m2; i += p->value)
						{
							info.values[i] += p->logp;
							info.values[i + delta] += p->logp;
						}
						if (i < poly.m2)
							info.values[i] += p->logp;
					}
				}
			}
			else
				for (auto p : info.base)
				{
					small_int a = poly.a % p->value;
					if (a == 0)
					{
						small_int b1 = (2 * poly.b) % p->value;
						if (b1 < 0)
							b1 += p->value;
						b1 = -(p->inverse[b1] * poly.c) % p->value ;
						small_int index = (p->value + poly.m + b1) % p->value;
						for (small_int i = index; i < poly.m2; i += p->value)
							info.values[i] += p->logp; 
					}
					else if (p->value > 2)
					{
						small_int a_1 = p->inverse[a];
						small_int a1 = (a_1 * a_1) % p->value;
						small_int t = info.n % p->value;
						small_int t0 = (a1 * (p->residue[t] - poly.b)) % p->value;
						small_int t1 = (a1 * (-p->residue[t] - poly.b)) % p->value;
						small_int index0 = (poly.m + t0 + p->value) % p->value;
						small_int index1 = (poly.m + t1 + p->value) % p->value;
						small_int index = std::min(index0, index1);
						small_int delta = std::max(index0, index1) - index;
						small_int m2 = poly.m2 - delta;
						small_int i = index;
						for (; i < m2; i += p->value)
						{
							info.values[i] += p->logp;
							info.values[i + delta] += p->logp;
						}
						if (i < poly.m2)
							info.values[i] += p->logp;
					}
					else
					{
						small_int m2 = poly.m2 - 1;
						for (int i = 1; i < m2; i += 2)
							info.values[i] += p->logp;
					}
				}
		}
		small_int collect_smooth(info_t &info, poly_t &poly)
		{
			real_t logp = static_cast<real_t>((info.config.sieve_offset + std::log(info.sqr2n) / 2 + std::log(poly.m)) * log_unit_e);
#ifdef HAVE_LARGE_PRIME
			small_int thrsp = info.base[info.config.base_size-1]->value;
			thrsp = thrsp * thrsp ;
			logp -= static_cast<real_t>(std::log(thrsp) * log_unit_e);
#endif
			for (small_int i = 0; i < poly.m2; i++)
				if (info.values[i] > logp)
				{
					large_int f = 0;
					large_int x = i - poly.m;
					if (poly.a == 1)
					{
						x += info.sqrn;
						f = x * x - info.n;
					}
					else
					{
						large_int f1 = poly.a2 * x + 2 * poly.b; // evaluate polynomial
						f = f1 * x + poly.c;
					}
					factors_t bit = 1;
					bit <<= info.factors_index;
					smooth_t smooth(poly.a2 * x + poly.b);
					if (f < 0)
					{
						f = -f;
						info.factors[0] |= bit;
					}
					for (size_t k = 0; k < info.config.base_size; k++)
					{
						int count = 0;
						small_int p = info.base[k]->value;
						while (f % p == 0)
						{
							f /= p;
							count++;
						}
						if (count & 1)
							info.factors[k + 1] |= bit;
						for (int j = 1; j < count; j += 2)
							smooth.sqr *= p;
					}
					stats_.processed();
					if (f == 1)
					{
						stats_.direct();
						info.smooths.push_back(smooth);
					}
#ifdef HAVE_LARGE_PRIME
					else if (f < thrsp)
					{
						stats_.large_primes();
						small_int f1 = static_cast<small_int>(f);
						auto it = info.large_primes.find(f1);
						if (it == info.large_primes.end())
							info.large_primes[f1] = smooth;
						else
						{
							info.smooths.push_back(smooth); // ok, some extra processing is required
							stats_.promoted();
						}
					}
#endif
					else
						stats_.discarded();
					if (info.smooths.size() > info.config.base_size + 2)
						return 2; ;
					//std::cout << i << ": value = " << static_cast<int>(info.values[i]) << " f = " << f << std::endl;
				}
			return 1;
		}
		typedef typename std::vector<prime_t>::iterator prime_it_t;
		typedef std::pair<prime_it_t, prime_it_t> primes_range_t;
		stat_t					stats_;
		std::vector<prime_t>	primes_;
		multiplier_t			multiplier_;
		primes_range_t			primes_range_;
#ifdef HAVE_MULTIPLIER_STATS
		mutable std::ofstream	mstats_;
#endif
	};

	template <class T>
	T mul_mod(const T &x, const T &y, const T &m)
	{
		return (x * y) % m;
	}

	uint64_t mul_mod1(const uint64_t &x, const uint64_t &y, const uint64_t &m)
	{
		uint32_t x1[2] = { (uint32_t)x, (uint32_t)(x >> 32)};
		uint64_t y1[2] = { (uint32_t)y, (uint32_t)(y >> 32)};
		uint64_t z0 = x1[0] * y1[0];
		uint64_t z1 = x1[0] * y1[1] + x1[1] * y1[0] ;
		uint64_t z00 = z0 + ((z1 & 0xFFFFFFFF) << 32); // overflow here?
		uint16_t *z01 = (uint16_t *)&z00;
		uint64_t z2 = x1[1] * y1[1] + (z1 >> 32) ;
		if (z00 < z0) // overflow on z00!
			z2++;
		z2 = z2 % m;
		z2 = ((z2 << 16) + z01[3]) % m;
		z2 = ((z2 << 16) + z01[2]) % m;
		z2 = ((z2 << 16) + z01[1]) % m;
		z2 = ((z2 << 16) + z01[0]) % m;
		return z2;
	}

	int64_t mul_mod(const int64_t &x, const int64_t &y, const int64_t &m)
	{
		return mul_mod1(x, y, m);
	}
	struct test_info_t
	{
		const char *file = nullptr;
		int algo = 0;
		int seeds = 0;
		int count = 0;
		int useed = 0;
		int bound = 10000;
		config_qs_t config_qs;
	};
	template <class large_int, class small_int>
	class prime_tester_t
	{
	public:
		struct input_item_t
		{
			large_int	n;
			bool		declared_prime;
			std::vector<small_int>	factors;
		};

		template <class stream_t>
		static input_item_t read_item(stream_t &ist)
		{
			char ch;
			input_item_t item;

			ist >> item.n >> ch;
			item.declared_prime = (ch == 'p');
			if (!item.declared_prime)
				do
				{
					small_int tmp;
					ist >> tmp;
					if (ist)
						item.factors.push_back(tmp);
					ist >> ch;
				} while (ist);

			return item;
		}
		typedef boost::multiprecision::cpp_int cpp_int_t;
		static bool fermat_test(const large_int &n)
		{
			large_int nm1 = n - 1;
#ifdef HAVE_CHECK
			large_int d = 8447990824LL;
			d = mul_mod(d, d, n);
			cpp_int_t n1 = n;
			cpp_int_t d1 = 8447990824LL;
			d1 = mul_mod(d1, d1, n1);
			cpp_int_t q1 = q;
			cpp_int_t x1 = 1;
#endif
			//
			// Begin with a single Fermat test - it excludes a lot of candidates:
			//
			large_int q(228); // We know n is greater than this, as we've excluded small factors
			large_int x = 1;
			for (; nm1; nm1 >>= 1)
			{
				if (nm1 & 1)
				{
					x = mul_mod(q, x, n);
#ifdef HAVE_CHECK
					x1 = mul_mod(q1, x1, n1);
#endif
				}
				q = mul_mod(q, q, n);
#ifdef HAVE_CHECK
				q1 = mul_mod(q1, q1, n1);
				std::cout << "q = " << q << ", q1 = " << q1 << ", x = " << x << ", x1 = " << x << std::endl;
#endif
			}
			return x == 1;
		}

		static bool custom_prime_test(const large_int &n, int *seeds, int count)
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
		//	int seeds[] = { 2, 5, 19, 41 };
			for (int i = 0; i < count; i++)
			{
				auto x = seeds[i];
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

		static std::vector<input_item_t> load(const char *file)
		{
			std::vector<input_item_t> result;
			std::ifstream ifile(file);
			while (ifile)
			{
				std::string line;
				std::getline(ifile, line);
				if (line.empty())
					break;
				std::istringstream ist(line);
				result.push_back(read_item(ist));
			}
			return result;
		}

		static void full_test(const std::vector<input_item_t> &items, large_int largest, large_int smallest, int count)
		{
			small_int bound = safe_cast<small_int>(sqrt(largest));
			auto primes = eratosthenes_sieve<small_int>(bound);
			small_int begin = safe_cast<small_int>(sqrt(smallest));
			auto it = std::lower_bound(primes.begin(), primes.end(), begin - 1);
			primes.erase(primes.begin(), it);
			int errors = 0;
			int composite = 0;
			int declared = 0;
			int processing = 0;
			for (const auto &item : items)
			{
				if (++processing % 1000 == 0)
					std::cout << "Processing " << (processing / (items.size() / 100)) << "%\r" << std::flush;
				if (item.declared_prime)
				{
					declared++;
					for (auto p : primes)
						if (item.n % p == 0)
						{
							errors++;
							break;
						}
				}
				else
					composite++;
				if (count && processing > count)
					break;
			}
			std::cout << "Examined: " << items.size() << "\n"
				<< "Declared prime = " << declared << "\n"
				<< "Errors = " << errors << "\n"
				<< "composite = " << composite << "\n";
		}
		static void brute_force(const std::vector<input_item_t> &items, large_int largest, large_int smallest, int count)
		{
			std::cout << "Brute force factorization\n";
			small_int bound = safe_cast<small_int>(sqrt(largest));
			auto primes = eratosthenes_sieve<small_int>(bound);
			small_int begin = safe_cast<small_int>(sqrt(smallest));
			auto it = std::lower_bound(primes.begin(), primes.end(), begin - 1);
			primes.erase(primes.begin(), it);
			int examined = 0;
			int factored = 0;
			for (const auto &item : items)
				if (!item.declared_prime)
				{
					examined++;
					for (auto p : primes)
						if (item.n % p == 0)
						{
							factored++;
							break;
						}
					if (examined % 1000 == 0)
						std::cout << "Examined " << examined << "%\r" << std::flush;
					if (count && examined >= count)
						break;
				}
			std::cout << "Examined: " << examined << ", factored = " << factored << "\n";
		}
		static int pollard_p1(const std::vector<input_item_t> &items, int *seeds,  const test_info_t &info)
		{
			pollard_p1_t<large_int, small_int> p1;
			p1.init(info.bound);
			int examined = 0;
			int factored = 0;
			std::cout << "Pollard p-1 factorization\n";
			if (info.useed > 0)
				seeds[0] = info.useed;
			for (const auto &item : items)
				if (!item.declared_prime)
				{
					examined++;
					if (examined % 1000 == 0)
						std::cout << "Processing " << (examined / (items.size() / 100)) << "%\r" << std::flush;
					for (int i = 0; i < info.seeds; i++)
					{
						auto f = p1.fact(item.n, seeds[i]);
						if (f != 1 && f != item.n)
						{
//							small_int f1 = safe_cast<small_int>(f);
//							small_int f2 = safe_cast<small_int>(item.n / f1);
							factored++;
							break;
						}
					}
					if (info.count && examined >= info.count)
						break;
				}
			std::cout << "Examined = " << examined << ", factored = " << factored << "\n";
			return examined;
		}
		static void pollard_rho(const std::vector<input_item_t> &items, int *seeds, const test_info_t &info)
		{
			int examined = 0;
			int factored = 0;
			std::cout << "Pollard rho factorization\n";
			if (info.useed > 0)
				seeds[0] = info.useed;
			for (const auto &item : items)
				if (!item.declared_prime)
				{
					examined++;
					if (examined % 1000 == 0)
						std::cout << "Processing " << (examined / (items.size() / 100)) << "%\r" << std::flush;
					for (int i = 0; i < info.seeds; i++)
					{
						large_int f = pollards_rho(item.n, info.bound, seeds[i]);
						if (f != 1 && f != item.n)
						{
							small_int f1 = safe_cast<small_int>(f);
							small_int f2 = safe_cast<small_int>(item.n / f1);
							factored++;
							break;
						}
					}
					if (info.count && examined >= info.count)
						break;
				}
			std::cout << "Examined = " << examined << ", factored = " << factored << "\n";
		}
		template <class stat>
		static int quadratic_sieve(const std::vector<input_item_t> &items, large_int largest, const test_info_t &info)
		{
			int examined = 0;
			int factored = 0;
			std::cout << "Quadratic sieve\n";
			quadratic_sieve_cached_t<stat, long long> qs(128);
			quadratic_sieve_cached_t<stat, long long>::info_t qsinfo;
			qsinfo.config = info.config_qs;
			for (const auto &item : items)
				if (!item.declared_prime)
				{
					examined++;
					if (examined % 1000 == 0)
						std::cout << "Processing " << (examined / (items.size() / 100)) << "%\r" << std::flush;
					if (qs.factor(static_cast<long long>(item.n), qsinfo) > 1)
						factored++;
					if (info.count && examined >= info.count)
						break;
				}
			std::cout << "\nExamined = " << examined << ", factored = " << factored << "\n";
			return examined;
		}
		static void elliptic_curve_1(const std::vector<input_item_t> &items, int *seeds, const test_info_t &info)
		{
			int examined = 0;
			int factored = 0;
			std::cout << "Elliptic curve factorization\n";
			if (info.useed > 0)
				seeds[0] = info.useed;
			elliptic_curve_t<large_int> ec(info.bound);
			for (const auto &item : items)
				if (!item.declared_prime)
				{
					examined++;
					if (examined % 1000 == 0)
						std::cout << "Processing " << (examined / (items.size() / 100)) << "%\r" << std::flush;
					elliptic_curve_t<large_int>::point_t pt{ seeds[1], seeds[2] };
					if (ec.run(item.n, seeds[0], pt))
						factored++;

					if (info.count && examined >= info.count)
						break;
				}
			std::cout << "Examined = " << examined << ", factored = " << factored << "\n";
		}
		static int prime_test(const std::vector<input_item_t> &items, int *seeds, int seed_count, int count)
		{
			int examined = 0;
			int errors = 0;
			std::cout << "prime test with " << seed_count << " seeds" << std::endl;
			for (const auto &item : items)
			{
				examined++;
				if (examined % 1000 == 0)
					std::cout << "Processing " << (examined / (items.size() / 100)) << "%\r" << std::flush;
				bool prime = custom_prime_test(item.n, seeds, seed_count);
				if (prime != item.declared_prime)
					errors++;
				if (count && examined >= count)
					break;
			}
			std::cout << "Examined = " << examined << ", errors = " << errors << "\n";
			return examined;
		}
		static int fermat_prime_test(const std::vector<input_item_t> &items, int count)
		{
			int examined = 0;
			int errors = 0;
			for (const auto &item : items)
			{
				examined++;
				if (examined % 1000 == 0)
					std::cout << "Processing " << (examined / (items.size() / 100)) << "%\r" << std::flush;
				bool prime = fermat_test(item.n);
				if (prime != item.declared_prime)
					errors++;
				if (count && examined >= count)
					break;
			}
			std::cout << "FExamined = " << examined << ", errors = " << errors << "\n";
			return examined;
		}

		static void process(const test_info_t &info)
		{
			std::cout << "Reading file.." << std::endl;
			auto items = load(info.file);
			large_int largest = 0;
			large_int smallest = std::numeric_limits<small_int>::max();
			for (const auto &item : items)
			{
				if (item.n > largest)
					largest = item.n;
				else if (item.n < smallest)
					smallest = item.n;
			}
			std::cout << "Start processing: range = [" << smallest << ", " << largest << "]" << std::endl;
			std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
			int seeds[] = { 2, 5, 19, 41, 67, 79 };
			double examined = static_cast<double>(info.count ? info.count : items.size());
			switch (info.algo)
			{
			case 0:
				if (info.seeds >= 0 && info.seeds <= 6)
					examined = prime_test(items, seeds, info.seeds, info.count);
				break;
			case 1:
				full_test(items, largest, smallest, info.count);
				break;
			case 2:
				brute_force(items, largest, smallest, info.count);
				break;
			case 3:
				examined = pollard_p1(items, seeds, info);
				break;
			case 4:
				pollard_rho(items, seeds, info);
				break;
			case 5:
				elliptic_curve_1(items, seeds, info);
				break;
			case 6:
				examined = quadratic_sieve<stats_t>(items, largest, info);
				break;
			case 61:
				examined = quadratic_sieve<stats_none_t>(items, largest, info);
				break;
			case 7:
				examined = fermat_prime_test(items, info.count);
				break;
			}
			auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
			std::cout << "Average iteration time = " << elapsed / examined << "us\n";
		}
	};
}

int main(int argc, char *argv[])
{
	bool use_cppint = true;
	zn::test_info_t info;
	for (int i = 0; i < argc; i++)
		if (strncmp(argv[i], "--file=", 7) == 0)
			info.file = argv[i] + 7;
		else if (strncmp(argv[i], "--seeds=", 8) == 0)
			info.seeds = atoi(argv[i] + 8);
		else if (strncmp(argv[i], "--useed=", 8) == 0)
			info.useed = atoi(argv[i] + 8);
		else if (strncmp(argv[i], "--bound=", 8) == 0)
			info.bound = info.config_qs.m2 = atoi(argv[i] + 8);
		else if (strncmp(argv[i], "--algo=", 7) == 0)
			info.algo = atoi(argv[i] + 7);
		else if (strncmp(argv[i], "--count=", 8) == 0)
			info.count = atoi(argv[i] + 8);
		else if (strncmp(argv[i], "--offset=", 9) == 0)
			info.config_qs.sieve_offset = atof(argv[i] + 9);
		else if (strncmp(argv[i], "--pquality=", 11) == 0)
			info.config_qs.pquality = atof(argv[i] + 11);
		else if (strncmp(argv[i], "--base-size=", 12) == 0)
			info.config_qs.base_size = atoi(argv[i] + 12);
		else if (strcmp(argv[i], "--more-polys") == 0)
			info.config_qs.more_polys = true;
		else if (strcmp(argv[i], "--use-long") == 0)
			use_cppint = false;
	if (info.file)
		if (use_cppint)
			zn::prime_tester_t<boost::multiprecision::cpp_int, long long>::process(info);
		else
			zn::prime_tester_t<long long, long>::process(info);
	else
		std::cout << "no file specified\n";
	return 0;
}