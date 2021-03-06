#ifndef znqssmall_H
#define znqssmall_H
/*
	A quadratic sieve implementation optimazed for small numbers (up to 60 bits)
*/
#include "znmultiplier.h"

namespace zn
{
	class stats_t
	{
	public:
		void processed(void) { processed_++; }
		void direct(void) { direct_++; }
		void discarded(void) { discarded_++; }
		void attempts(void) { attempts_++; }
		void attempts_log(void)
		{
			int delta = std::min<int>(attempts_ - attempts_ref_ - 1, sizeof(attempts_hist_) / sizeof(int) - 1);
			attempts_hist_[delta]++;
			attempts_ref_ = attempts_;
		}
		void promoted(void) { promoted_++; }
		void polynomials(void) { polynomials_++; }
		void large_primes(void) { large_primes_++; }
		template <class large_int, class small_int>
		void report(const large_int& n, const small_int& p)
		{
			//if (p == 1)
			//	std::cout << "\t" << n << std::endl;
		}
		void resume(void)
		{
			std::cout << "Processed = " << processed_ << "\n"
				<< "Discarded = " << discarded_ << "\n"
				<< "Attempts = " << attempts_ << "\n"
#ifdef HAVE_LARGE_PRIME
				<< "Direct = " << direct_ << "\n"
				<< "Large = " << large_primes_ << "\n"
				<< "Promoted = " << promoted_ << "\n"
#endif
				<< "Accepted = " << (processed_ - discarded_) << "\n"
				<< "Polynomials = " << polynomials_ << "\n"
				;
			for (auto attempts : attempts_hist_)
				std::cout << attempts << " ";
			std::cout << std::endl;
		}
	private:
		int processed_ = 0;
		int promoted_ = 0;
		int discarded_ = 0; // smooth detection failures
		int large_primes_ = 0;
		int direct_ = 0;
		int attempts_ = 0;
		int polynomials_ = 0;
		int attempts_ref_ = 0;
		int attempts_hist_[16] = { 0 };
	};
	struct stats_none_t
	{
	public:
		void processed(void) {}
		void direct(void) {}
		void discarded(void) {}
		void promoted(void) {}
		void attempts(void) { }
		void attempts_log(void) {}
		void polynomials(void) {}
		void large_primes(void) {}
		template <class large_int, class small_int>
		void report(const large_int&, const small_int&) {}
		void resume(void) {}
	};
	struct interpolator_t
	{
		float a;
		float b;
	};
	int interpolate(const interpolator_t& interp, int bits)
	{
		return static_cast<int>(interp.a + bits * interp.b);
	}


	template <class T, int N>
	struct bitmask_op_t
	{
		static void init(T* data, T value)
		{
			data[0] = value;
			bitmask_op_t<T, N - 1>::init(data + 1, value);
		}
		static bool is_zero(const T* rhs)
		{
			if (rhs[0])
				return false;
			return bitmask_op_t<T, N - 1>::is_zero(rhs + 1);
		}
		static unsigned int eval_msb(const T* data)
		{
			if (data[N - 1])
				return sizeof(T) * 8 + msb(data[N - 1]);
			else
				return bitmask_op_t<T, N - 1>::eval_msb(data);
		}
		static void eval_or(T* lhs, const T* rhs)
		{
			lhs[0] |= rhs[0];
			bitmask_op_t<T, N - 1>::eval_or(lhs + 1, rhs + 1);
		}
		static void eval_xor(T* lhs, const T* rhs)
		{
			lhs[0] ^= rhs[0];
			bitmask_op_t<T, N - 1>::eval_xor(lhs + 1, rhs + 1);
		}
		static void eval_and(T* lhs, const T* rhs)
		{
			lhs[0] &= rhs[0];
			bitmask_op_t<T, N - 1>::eval_and(lhs + 1, rhs + 1);
		}
		static void eval_not(T* data)
		{
			data[0] = ~data[0];
			bitmask_op_t<T, N - 1>::eval_not(data + 1);
		}
	};

	template <class T>
	struct bitmask_op_t<T, 0>
	{
		static void init(T*, T) {}
		static bool is_zero(const T*) { return true; }
		static unsigned int eval_msb(const T*) { return 0; }
		static void eval_or(T*, const T*) {}
		static void eval_xor(T*, const T*) {}
		static void eval_and(T*, const T*) {}
		static void eval_not(T*) {}
	};

	template <class T, int N>
	class bitmask_t
	{
		enum { bit_size = sizeof(T) * 8 };
	public:
		bitmask_t(T t = 0)
		{
			data_[0] = t;
			bitmask_op_t<T, N - 1>::init(data_ + 1, 0);
		}
		bool is_set(size_t bit) const
		{
			return (data_[bit / bit_size] & (1 << (bit % bit_size))) != 0;
		}
		bool empty(void) const
		{
			return bitmask_op_t<T, N>::is_zero(data_);
		}
		unsigned int msb(void) const
		{
			return bitmask_op_t<T, N>::eval_msb(data_);
		}
		void set(size_t bit)
		{
			data_[bit / bit_size] |= 1 << (bit % bit_size);
		}
		void unset(size_t bit)
		{
			data_[bit / bit_size] &= ~(1 << (bit % bit_size));
		}
		bitmask_t<T, N>& operator=(T value)
		{
			data_[0] = value;
			return *this;
		}
		bitmask_t<T, N>& operator|=(const bitmask_t<T, N>& rhs)
		{
			bitmask_op_t<T, N>::eval_or(data_, rhs.data_);
			return *this;
		}
		bitmask_t<T, N>& operator^=(const bitmask_t<T, N>& rhs)
		{
			bitmask_op_t<T, N>::eval_xor(data_, rhs.data_);
			return *this;
		}
		bitmask_t<T, N>& operator&=(const bitmask_t<T, N>& rhs)
		{
			bitmask_op_t<T, N>::eval_and(data_, rhs.data_);
			return *this;
		}
		bitmask_t<T, N>& operator~(void)
		{
			bitmask_op_t<T, N>::eval_not(data_);
			return *this;
		}
	private:
		T data_[N];
	};

	template <class T, int N>
	bitmask_t<T, N> operator&(const bitmask_t<T, N>& lhs, const bitmask_t<T, N>& rhs)
	{
		bitmask_t<T, N> result(lhs);
		return result &= rhs;

	}

	template <class T, int N>
	bool bit_test(const bitmask_t<T, N>& value, unsigned int bit)
	{
		return value.is_set(bit);
	}

	template <class T, int N>
	void bit_set(bitmask_t<T, N>& value, unsigned int bit)
	{
		value.set(bit);
	}

	template <class T, int N>
	void bit_unset(bitmask_t<T, N>& value, unsigned int bit)
	{
		value.unset(bit);
	}

	template <class T, int N>
	bool is_zero(const bitmask_t<T, N>& value)
	{
		return value.empty();
	}

	template <class T, int N>
	unsigned int msb(const bitmask_t<T, N>& value)
	{
		return value.msb();
	}

	struct config_qs_t
	{
		int		m2 = 1000;
		size_t	base_size = 32;
		double	sieve_offset = 0.6;
		double	pquality = 5.0; // to select primes for 'a' coefficient
		int		smooth_excess = 2; // number of smooths to find
		bool	more_polys = false;
		interpolator_t base_int = { 24.0f, 0.4f };// actual baseof n is a + b * msb(n), another good is {18, 0.18}
		interpolator_t m2_int = { 800.0f, 40.0f };
	};
	// A quadratic sieve class for factoring numbers up to (about) 19 (2^63) digits
	template <class stat_t = stats_none_t, class large_int = long long, class small_int = int>
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
					logk_[i * i] = 0;
			}
			template <class PrimeIt>
			std::pair<int, float> find_best(const large_int& n, const std::pair<PrimeIt, PrimeIt>& range) const
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
			float quality(const large_int& n, int k, const std::pair<PrimeIt, PrimeIt>& range) const
			{
				large_int kn = k * n;
				float fk = -logk_[k] / 2;
				for (auto iter = range.first; iter != range.second; ++iter)
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
		//typedef bitmask_t<unsigned int, 3> factors_t;
		struct smooth_t
		{
			large_int sqr; // product of prime with even exponents (/ 2)
			large_int axb;// a *x + b, the number squared
			factors_t factors;
			smooth_t(large_int ab = 1, large_int s = 1) : axb(ab), sqr(s), factors(0) {}
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
		struct factor_info_t
		{
			int			value;
			factors_t	mask;
			factor_info_t(int v) : value(v), mask(0) {}
		};
#define HAVE_MONTGOMERY
		struct info_t
		{
#ifdef HAVE_MONTGOMERY
			typedef montgomery_t<large_int> multiplier_t;
			multiplier_t				 mul;
#endif
			large_int					 n;
			int							 multi; // multiplier, most times 1
			std::vector<const prime_t*> base;
			std::vector<real_t>			 values;
			std::vector<smooth_t>		 smooths;
			std::vector<smooth_t>		 smooths_tmp;
			std::vector<factor_info_t>	 factors;
			std::vector<int>			 smooth_perm;
			std::vector<int>			 smooth_index;
			size_t						 smooth_required;
#ifdef HAVE_LARGE_PRIME
			std::map<small_int, smooth_t> large_primes;
#endif
			std::vector<poly_seed_t>	 poly;
			small_int					 sqr2n; // sqrt(n * 2)
			small_int					 sqrn; // sqrt(n)
			config_qs_t					 config;
			factors_t perm_bit(size_t j) const { return static_cast<factors_t>(1) << smooth_perm[j]; }
		};
		quadratic_sieve_cached_t(small_int primes = 256) : multiplier_(12)
		{
			small_int bound = static_cast<small_int>(primes * std::pow(std::log(primes), 1.3));
			std::vector<small_int> prime = eratosthenes_sieve<small_int>(bound);
			if (prime.size() > static_cast<size_t>(primes))
				prime.erase(prime.begin() + primes, prime.end());
			for (auto p : prime)
				primes_.push_back(prime_t(p));
			primes_range_ = std::make_pair(primes_.begin(), primes_.begin() + std::min<size_t>(10, primes_.size()));
#ifdef HAVE_MULTIPLIER_STATS
			mstats_.open("mstats.txt");
#endif
		}
		~quadratic_sieve_cached_t(void)
		{
			stats_.resume();
		}
		small_int factor(const large_int& n, info_t& info) const
		{
			small_int result = factor_imp(n, info, false);
			if (result == 1)
				result = factor_imp(n, info, true);
			stats_.report(n, result);
			return result;
		}
		small_int factor_imp(const large_int& n, info_t& info, bool with_multi = false) const
		{
			small_int result = 1;
			large_int n1 = static_cast<large_int>(sqrt(n + 0.5));
			if (n1 * n1 == n)
				return static_cast<small_int>(n1); // a perfect square: case handled here because below is a problem
			init_info(n, info, with_multi);
			std::sort(info.poly.begin(), info.poly.end(), [](const poly_seed_t& lhs, const poly_seed_t& rhs) {
				return lhs.quality < rhs.quality;
				});
#ifdef HAVE_MULTIPLIER_STATS
			int poly_count = 0;
#endif
			for (const auto& poly_seed : info.poly)
			{
				poly_t poly = create_poly(info, poly_seed);
				if (poly.b == 0 && poly.a > 1)
					continue;
				stats_.polynomials();
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
			mstats_ << multiplier_.quality(n, 1, primes_range_) << " " << poly_count << "\n";
#endif
			return result;
		}
	private:
		poly_t create_poly(info_t& info, const poly_seed_t& seed) const
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
				poly.a2 = 1;
				poly.b = 0;
				poly.c = info.n;
				poly.m = info.config.m2 * 4;
			}
			poly.m2 = poly.m * 2;
			return poly;
		}
		void init_info(large_int n, info_t& info, bool with_multi) const
		{
			info.n = n;
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
#ifdef HAVE_MONTGOMERY
			info.mul.init(info.n);
#endif
			int bits = msb(info.n);
			info.config.base_size = interpolate(info.config.base_int, bits);
			info.config.m2 = interpolate(info.config.m2_int, bits);
			info.base.clear();
			info.smooths.clear();
			info.poly.clear();
#ifdef HAVE_LARGE_PRIME
			info.large_primes.clear();
#endif
			info.sqrn = static_cast<small_int>(std::sqrt(info.n));
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
							poly_seed_t poly{ q, it->value, 1 };
							info.poly.push_back(poly);
						}
					}
					else
						break;
			info.factors.clear();
			info.factors.push_back(factor_info_t(-1));
			for (auto it : info.base)
				info.factors.push_back(it->value);
			info.smooth_required = info.config.base_size + info.config.smooth_excess;
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

			poly_seed_t poly{ 2 * info.config.pquality + 1, -1 }; // worse quality
			info.poly.push_back(poly);
		}
		void sieve(const poly_t& poly, info_t& info) const
		{
			if (poly.a == 1)
			{
				for (auto p : info.base)
				{
					small_int t = p->residue[info.n % p->value];
					small_int r = 2 * p->value - ((info.sqrn - poly.m) % p->value);
					small_int i0 = (r + t) % p->value;
					if (p->value > 2)
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
						b1 = -(p->inverse[b1] * poly.c) % p->value;
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
		small_int collect_smooth(info_t& info, poly_t& poly) const
		{
			real_t logp = static_cast<real_t>((info.config.sieve_offset + std::log(info.sqr2n) / 2 + std::log(poly.m)) * log_unit_e);
#ifdef HAVE_LARGE_PRIME
			small_int thrsp = info.base[info.config.base_size - 1]->value;
			thrsp = thrsp * thrsp;
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
					unsigned int bit = static_cast<unsigned>(info.smooths.size());
					smooth_t smooth(poly.a2 * x + poly.b, poly.a);
					if (f < 0)
					{
						f = -f;
						bit_set(info.factors[0].mask, bit);
						smooth.factors = 1;
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
						{
							bit_set(info.factors[k + 1].mask, bit);
							bit_set(smooth.factors, static_cast<unsigned>(k + 1));
						}
						else
							bit_unset(info.factors[k + 1].mask, bit);
						for (int j = 1; j < count; j += 2)
							smooth.sqr *= p;
					}
					stats_.processed();
					if (f == 1)
					{
						stats_.direct();
#ifdef _DEBUG
						smooth_invariant(info, smooth, false);
#endif
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
					if (info.smooths.size() >= static_cast<size_t>(info.smooth_required))
					{
						small_int result = solve(info);
						if (result != 1) // collect some more relations
							return result;
					}
					//std::cout << i << ": value = " << static_cast<int>(info.values[i]) << " f = " << f << std::endl;
				}
			return 1;
		}
		small_int build_result(info_t& info, int start_index) const
		{
			size_t sall = info.smooth_perm.size();
			size_t size = info.factors.size();
			for (size_t i = start_index; i < sall; i++)
			{
				int ip = info.smooth_perm[i];
				smooth_t s = info.smooths[ip]; // should pick up its list of factors
#ifdef _DEBUG
				smooth_invariant(info, s, false);
#endif
#ifdef HAVE_MONTGOMERY
				s.axb = info.mul.project(s.axb);
				s.sqr = info.mul.project(s.sqr);
#endif
				for (int j = 0; j < size; j++)
					if (bit_test(info.factors[j].mask, ip))
					{
						int sindex = info.smooth_index[j];
						const smooth_t& s1 = info.smooths[sindex];
#ifdef HAVE_MONTGOMERY
						s.axb = info.mul.prod(s.axb, info.mul.project(s1.axb));
						s.sqr = info.mul.prod(s.sqr, info.mul.project(s1.sqr));
						//s.axb = info.mul.mul(s.axb, s1.axb);
						//s.sqr = info.mul.mul(s.sqr, s1.sqr);
#else
						s.axb = mul_mod(s.axb, s1.axb, info.n);
						s.sqr = mul_mod(s.sqr, s1.sqr, info.n);
#endif
						factors_t common = s.factors;
						common &= s1.factors;
						int common_bits = msb(common);
						s.factors ^= s1.factors;
						for (int k = 0; k <= common_bits; k++)
							if (bit_test(common, k))
#ifdef HAVE_MONTGOMERY
								s.sqr = info.mul.prod(s.sqr, info.mul.project(info.factors[k].value));
								//s.sqr = info.mul.mul(s.sqr, info.factors[k].value);
#else
								s.sqr = mul_mod<large_int>(s.sqr, info.factors[k].value, info.n);
#endif

#if defined(_DEBUG)
						smooth_invariant(info, s, true);
#endif
					}
#ifdef HAVE_MONTGOMERY
				s.axb = info.mul.reduce(s.axb);
				s.sqr = info.mul.reduce(s.sqr);
#endif
				stats_.attempts();
				large_int n = info.n / info.multi;
				large_int q = euclidean_algorithm(n, s.sqr + s.axb);
				if (q != 1 && q != n)
				{
					stats_.attempts_log();
					auto p = n / q;
					return static_cast<small_int>(std::min(p, q));
				}
#if 0
				q = euclidean_algorithm(n, abs(s.sqr - s.axb));
				if (q != 1 && q != n)
					return static_cast<small_int>(q);
#endif
			}
			// some kind of backtracking: reuse first start_index smooth 
			// and set factors bitmask accordingly
			for (int i = 0; i < size; i++)
				info.factors[i].mask = 0;
			info.smooths_tmp.clear();
			for (int i = 0; i < start_index; i++)
			{
				auto& f = info.smooths[info.smooth_index[i]];
				info.smooths_tmp.push_back(f);
				factors_t bit(0);
				bit_set(bit, i);
				auto m = f.factors;
				int j0 = msb(m);
				for (int j = 0; j <= j0; j++)
					if (bit_test(m, j))
						info.factors[j].mask |= bit;
			}
			info.smooths = info.smooths_tmp;
			return 1;
		}
		void print(const info_t& info) const {}
		void print1(const info_t& info) const
		{
			std::cout << info.smooth_perm << "\n" << info.smooth_index << std::endl;
			int sall = info.smooths.size();
			int size = static_cast<int>(info.factors.size());
			for (int i = 0; i < size; i++)
			{
				factors_t f = info.factors[i].mask;
				for (int j = 0; j < sall; j++)
					std::cout << (((f & (1LL << j))) ? '1' : '0') << (j < sall - 1 ? " " : "\n");
			}
		}
		void smooth_invariant(const info_t& info, const smooth_t& smooth, bool red) const
		{
			large_int axb = smooth.axb;
			large_int sqr = smooth.sqr;
			if (red)
			{
#if defined(HAVE_MONTGOMERY)
				axb = info.mul.reduce(axb);
				sqr = info.mul.reduce(sqr);
#endif
			}
			large_int q1 = mul_mod(axb, axb, info.n);
			large_int q2 = mul_mod(sqr, sqr, info.n);
			auto mask = smooth.factors;
			int fmsb = msb(mask);
			for (int index = 0; index <= fmsb; index++)
				if (bit_test(mask, index))
					q2 = mul_mod<large_int>(q2, info.factors[index].value, info.n);
			large_int r = (q1 - q2) % info.n;
			if (r)
				std::cout << "Hmmm: " << info.n << ", " << info.multi << "\n";
		}
		void remove_unused_smooths(info_t& info, small_int& smooths_size) const
		{
			factors_t mask = 0;
			factors_t mask_end = static_cast<factors_t>(1) << (smooths_size - 1);
			for (auto& f : info.factors)
				if (bitcount(f.mask) == 1)
				{
					mask |= f.mask;
					f.mask = 0;
				}
			for (int i = 0; mask; i++, mask >>= 1)
				if (mask & 1)
				{
					while (mask & mask_end && i < smooths_size)
					{
						mask_end >>= 1;
						smooths_size--;
					}
					std::swap(info.smooth_perm[i], info.smooth_perm[--smooths_size]);
				}
			info.smooth_perm.erase(info.smooth_perm.begin() + smooths_size, info.smooth_perm.end());
		}
		small_int solve(info_t& info) const
		{
			small_int result = 1;
			int smooths_size = static_cast<int>(info.smooths.size());
			int base_size = static_cast<int>(info.factors.size());
			info.smooth_perm.resize(smooths_size);
			info.smooth_index.clear();
			for (int k = 0; k < smooths_size; k++)
				info.smooth_perm[k] = k;
			//remove_unused_smooths(info, smooths_size);
			int i = 0; // pointer to factors
			int j = 0; // pointer to smooths
			factors_t f;
			print(info);
			for (; i < base_size; i++)
				if (!is_zero(f = info.factors[i].mask))
				{
					unsigned int b;
					if (bit_test(f, b = info.smooth_perm[j]) == 0)
						for (int k = j + 1; k < smooths_size; k++)
							if (bit_test(f, b = info.smooth_perm[k]))
							{
								std::swap(info.smooth_perm[j], info.smooth_perm[k]);
								break;
							}
					// Gaussian elimination step
					for (int k = i + 1; k < base_size; k++)
						if (bit_test(info.factors[k].mask, b))
							info.factors[k].mask ^= f;
					info.smooth_index.push_back(info.smooth_perm[j]);
					j++;
					print(info);
				}
				else
					info.smooth_index.push_back(info.smooth_perm[j]);

			int j0 = j;
			for (--i, --j; i >= 0; i--)
				if (!is_zero(f = info.factors[i].mask))
				{
					int b = info.smooth_perm[j];
					for (int h = 0; h < i; h++)
						if (bit_test(info.factors[h].mask, b))
							info.factors[h].mask ^= f;
					j--;
					print(info);
				}
			return build_result(info, j0);
		}
		typedef typename std::vector<prime_t>::iterator prime_it_t;
		typedef std::pair<prime_it_t, prime_it_t> primes_range_t;
		mutable stat_t			stats_;
		std::vector<prime_t>	primes_;
		multiplier_t			multiplier_;
		primes_range_t			primes_range_;
#ifdef HAVE_MULTIPLIER_STATS
		mutable std::ofstream	mstats_;
#endif
	};

	template <class T>
	int bitcount(const T& t)
	{
		return 0;
	}
	const uint8_t bitcount_table[] = { 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 };
	int bitcount(const uint8_t& t)
	{
		return bitcount_table[t & 0x0F] + bitcount_table[t >> 4];
	}

	int bitcount(const uint16_t& t)
	{
		const uint8_t* t1 = reinterpret_cast<const uint8_t*>(&t);
		return bitcount(t1[0]) + bitcount(t1[1]);
	}

	int bitcount(const uint32_t& t)
	{
		const uint16_t* t1 = reinterpret_cast<const uint16_t*>(&t);
		return bitcount(t1[0]) + bitcount(t1[1]);
	}

	int bitcount(const long long& t)
	{
		const uint32_t* t1 = reinterpret_cast<const uint32_t*>(&t);
		return bitcount(t1[0]) + bitcount(t1[1]);
	}

	template <class T>
	struct mul_t
	{
		static T mod(const T& x, const T& y, const T& m)
		{
			return (x * y) % m;
		}
	};

	template <class T>
	T mul_mod(const T& x, const T& y, const T& m)
	{
		return mul_t<T>::mod(x, y, m);
	}

	template <>
	struct mul_t<uint64_t>
	{
#if 0
		static uint64_t mod(uint64_t a, uint64_t b, uint64_t c)
		{
			uint64_t d; /* to hold the result of a*b mod c */
						/* calculates a*b mod c, stores result in d */
			uint64_t dl = _umul128(a, b, &d);
			_asm {
				mov rax, qword ptr a;
				mul b;
				div c;
				mov[d], rdx;
			}
#if 0
			asm("mov %1, %%rax;"        /* put a into rax */
				"mul %2;"               /* mul a*b -> rdx:rax */
				"div %3;"               /* (a*b)/c -> quot in rax remainder in rdx */
				"mov %%rdx, %0;"        /* store result in d */
				:"=r"(d)                /* output */
				: "r"(a), "r"(b), "r"(c) /* input */
				: "%rax", "%rdx"         /* clobbered registers */
			);
#endif
			return d;
		}
#else
		static uint64_t mod(const uint64_t& x, const uint64_t& y, const uint64_t& m)
		{
			uint32_t x1[2] = { (uint32_t)x, (uint32_t)(x >> 32) };
			uint64_t y1[2] = { (uint32_t)y, (uint32_t)(y >> 32) };
			if (x1[1] == 0 && y1[1] == 0)
				return x1[0] * y1[0] % m;
			uint64_t z0 = x1[0] * y1[0];
			uint64_t z1 = x1[0] * y1[1] + x1[1] * y1[0];
			uint64_t z00 = z0 + ((z1 & 0xFFFFFFFF) << 32); // overflow here?
			uint64_t z2 = x1[1] * y1[1] + (z1 >> 32);
			if (z00 < z0) // overflow on z00!
				z2++;
			z2 = z2 % m;
			uint16_t* z21 = (uint16_t*)&m;
			if (z21[3] == 0)
			{
				uint16_t* z01 = (uint16_t*)&z00;
				z2 = ((z2 << 16) + z01[3]) % m;
				z2 = ((z2 << 16) + z01[2]) % m;
				z2 = ((z2 << 16) + z01[1]) % m;
				z2 = ((z2 << 16) + z01[0]) % m;
			}
			else
			{
				uint8_t* z01 = (uint8_t*)&z00;
				z2 = ((z2 << 8) + z01[7]) % m;
				z2 = ((z2 << 8) + z01[6]) % m;
				z2 = ((z2 << 8) + z01[5]) % m;
				z2 = ((z2 << 8) + z01[4]) % m;
				z2 = ((z2 << 8) + z01[3]) % m;
				z2 = ((z2 << 8) + z01[2]) % m;
				z2 = ((z2 << 8) + z01[1]) % m;
				z2 = ((z2 << 8) + z01[0]) % m;
			}
			return z2;
		}
#endif
	};

	template <>
	struct mul_t<int64_t>
	{
		static int64_t mod(const int64_t& x, const int64_t& y, const int64_t& m)
		{
			if (x >= 0)
				if (y >= 0)
					return mul_t<uint64_t>::mod(x, y, m);
				else
					return -static_cast<int64_t>(mul_t<uint64_t>::mod(x, -y, m));
			else
				if (y >= 0)
					return -static_cast<int64_t>(mul_t<uint64_t>::mod(-x, y, m));
				else
					return mul_t<uint64_t>::mod(-x, -y, m);
		}
	};

}
#endif
