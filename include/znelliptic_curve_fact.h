#ifndef znelliptic_curve_factH
#define znelliptic_curve_factH

#include "znbasic.h"

namespace zn
{
	template <class large_int, class small_int = long>
	class pollard_p1_t
	{
	public:
		pollard_p1_t(void) {}
		void init(small_int bound)
		{
			double logb = std::log(static_cast<double>(bound));
			int  b1 = safe_cast<int>(pow(bound, 0.5));
			std::vector<small_int> primes = eratosthenes_sieve<small_int, int>(b1);
			small_int exp1 = 1;
			for (auto prime : primes)
			{
				int exp = static_cast<int>(logb / std::log(static_cast<double>(prime)) + 0.5);
				for (int j = 0; j < exp; j++)
					if (exp1 < std::numeric_limits<small_int>::max() / prime)
						exp1 *= prime;
					else
					{
						exp_.push_back(exp1);
						exp1 = prime;
					}
			}
			if (exp1 != 1)
				exp_.push_back(exp1);
		}
		large_int fact(const large_int &n, small_int b)
		{
			large_int pw = b;
			size_t size = exp_.size();
			for (size_t i = 0 ; i < size ; i++)
			{
				pw = powm<large_int>(pw, exp_[i], n);
				if (i % 4 == 0)
				{
					large_int g = euclidean_algorithm<large_int>(pw - 1, n);
					if (g == n)
						return 1;
					else if (g > 1)
						return g;
				}
			}
			return euclidean_algorithm<large_int>(pw - 1, n);
		}
	private:
		std::vector<small_int> exp_;
	};

	template <class T, int N>
	struct rho_poly_t
	{
		static T apply(const T &t, const T &n)
		{
			return rho_poly_t<T, N - 1>::apply((t * t + 1) % n, n);
		}
	};
	template <class T>
	struct rho_poly_t<T, 1>
	{
		static T apply(const T &t, const T &n)
		{
			return (t * t + 1) % n;
		}
	};
	template <class large_int>
	large_int pollards_rho(const large_int &n, int count, int seed = 2)
	{
		large_int d;
		large_int x = seed;
		large_int y = x;
		const int iter = 5;
		for (int i = 0 ; i < count ; i++)
		{
			x = rho_poly_t<large_int, iter>::apply(x, n);
			y = rho_poly_t<large_int, iter * 2>::apply(y, n);
			d = gcd(abs(x - y), n);
			if (d != 1)
				break;
		}
		if (d < n)
			return d;
		else
			return 0;
	}

	template <class large_int>
	class elliptic_curve_projective_t
	{
	public:
		struct point_t
		{
			large_int x;
			large_int z;
		};
		elliptic_curve_projective_t(const large_int &n) : n_(n) 
		{}
		void init(const large_int &a, const point_t &pt)
		{
			a_ = a;
			auto ext = extended_euclidean_algorithm<large_int>(n_, 4);
			a24_ = ((a + 2) * std::get<2>(ext)) % n_ ;
#ifdef _DEBUG
			y_ = 1;
			auto x2 = (pt.x * pt.x) % n_;
			b_ = (x2 * (pt.x + a_) + pt.x) % n_;
#endif
		}
#ifdef _DEBUG
		void test(const point_t &pt)
		{
			point_t p2, p3, p4, p4a;
			//verify(pt, pt);
			add(pt, pt, p2);
			//verify(p2, pt);
			add(pt, p2, p3);
			//verify(p3, pt);
			add(pt, p3, p4);
			add(p2, p2, p4a);
		}
#endif
	private:
		void add(const point_t &p, const point_t &q, point_t &result)
		{
			if (p.x * q.z == q.x * p.z) // double
			{
				auto s = p.x + p.z;
				auto d = p.x - p.z;
				auto xz = (s * s - d * d) % n_;
				auto x1 = (s * d) % n_;
				auto z1 = (d * d + a24_ * xz) % n_;
				result.x = (x1 * x1) % n_;
				result.z = (xz * z1) % n_;
			}
			else  // add
			{
				auto dp = p.x - p.z;
				auto sp = p.x + p.z;
				auto dq = q.x - q.z;
				auto sq = q.x + q.z;
				auto dpsq = dp * sq;
				auto spdq = sp * dq;
				auto x1 = (dpsq + spdq) % n_;
				auto z1 = (dpsq - spdq) % n_;
				x1 = (x1 * x1) % n_;
				z1 = (z1 * z1) % n_;
				result.x = ((p.z - q.z) * x1) % n_;
				result.z = ((p.x - q.x) * z1) % n_;
			}
		}
#ifdef _DEBUG
		bool verify(const point_t &pt, const point_t &p1)
		{
			auto y1 = (pt.x * p1.x + 1) % n_;
			auto dx = pt.x - p1.x; 
			y1 = = y1 * (y1 + a_ - a) % n_ - 2 * a - 
			auto lhs = b_ * pt.z;
			auto xz = (pt.x * pt.z) % n_;
			auto x2 = (pt.x * pt.x) % n_;
			auto rhs = pt.x * x2 + xz * ((a_ * pt.x) % n_ + pt.z);
			return (lhs - rhs) % n_ == 0;
		}
		large_int b_;
		large_int y_;
#endif
		large_int a24_;
		large_int a_;
		large_int n_;
	};

	template <class large_int>
	class elliptic_curve_t
	{
	public:
		struct point_t
		{
			large_int x;
			large_int y;
		};
		elliptic_curve_t(const large_int &n) : n_(n){	}
		large_int factor(void) const { return f_; }
		bool run(const large_int &a, point_t &pt, const large_int &k)
		{
			point_t dummy;
			a_ = a;
			b_ = (pt.y * pt.y - pt.x * (pt.x * pt.x + a_)) % n_;
			bool result = kadd(pt, k, dummy);
			pt = dummy;
			return result;
		}
		// a = 5, n = 455839
		bool test(const large_int &a, const point_t &pt)
		{
			point_t p2, p4, p6;
			a_ = a;
			b_ = (pt.y * pt.y - pt.x * (pt.x * pt.x + a_)) % n_;
			add(pt, pt, p2);
			verify(p2);
			add(p2, p2, p4);
			verify(p4);
			add(p2, p4, p6);
			verify(p6);
			return p6.x != 1;
		}
	private:
		bool verify(const point_t &pt)
		{
			large_int y1 = pt.y * pt.y;
			large_int x1 = (pt.x * pt.x + a_) % n_;
			x1 = (x1 * pt.x + b_ - y1) % n_;
			return x1 == 0;
		}
		bool kadd(const point_t &pt, const large_int &k, point_t &result)
		{
			bool result_set = false;
			unsigned int bits = msb(k);
			point_t pwr = pt;
			for (unsigned int i = 0; i < bits; i++)
			{
				if (bit_test(k, i))
				{
					if (result_set)
					{
						if (add(pwr, result, result))
							return true;
						verify(result);
					}
					else
					{
						result_set = true;
						result = pwr;
					}
				}
				if (add(pwr, pwr, pwr))
					return true;
				verify(pwr);
			}
			return false;
		}
		bool add(const point_t &p, const point_t &q, point_t &result)
		{
			large_int s;
			if (p.x == q.x)
			{
				auto ext = extended_euclidean_algorithm<large_int>(n_, 2 * p.y);
				if ((f_ = abs(std::get<0>(ext))) != 1)
					return true;
				s = (3 * p.x * p.x + a_) % n_;
				s = (s * std::get<2>(ext)) % n_;
				if (std::get<0>(ext) < 0)
					s = -s;
			}
			else
			{
				auto ext = extended_euclidean_algorithm<large_int>(n_, p.x - q.x);
				if ((f_ = abs(std::get<0>(ext))) != 1)
					return true;
				s = ((p.y - q.y) * std::get<2>(ext)) % n_;
				if (std::get<0>(ext) < 0)
					s = -s;
			}
			auto resx = (s * s - p.x - q.x) % n_;
			result.y = (s * (p.x - resx) - p.y) % n_;
			result.x = resx;
			return false;
		}
		large_int f_;
		large_int a_;
		large_int b_;
		large_int n_; // number to factor
	};
}
#endif
