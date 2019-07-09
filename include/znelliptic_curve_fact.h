#ifndef znelliptic_curve_factH
#define znelliptic_curve_factH

#include "znbasic.h"

namespace zn
{
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
