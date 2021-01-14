#ifndef znmultiplier_H
#define znmultiplier_H

#ifdef _MSC_VER
#ifdef _M_X64 
#define HAVE_INTRINSIC
#endif
#endif
namespace zn
{
#ifdef HAVE_INTRINSIC
#include <intrin.h>
#endif

	template <class T>
	std::tuple<T, T> montgomery_inverse(const T& n)
	{
		T q;
		int i, j;
		T r0 = static_cast<T>(1) << (sizeof(T) * 8 - 1);
		T r[2] = { r0, n };
		T s[2] = { 1, 0 };
		T t[2] = { 0, 1 };

		for (i = 0, j = 1; r[j]; i = j, j = (j + 1) & 1)
		{
			q = r[i] / r[j];
			r[i] = r[i] - q * r[j];
			s[i] = s[i] - q * s[j];
			t[i] = t[i] - q * t[j];
		}
		if (s[i] & 1)
		{
			t[i] += r0;
			s[i] = (s[i] + n) / 2;
		}
		else
			s[i] /= 2;
		return std::tuple<T, T>(-r[i] * s[i], -r[i] * t[i]);
	}

#ifdef HAVE_INTRINSIC
	void multiplication_unsigned(const std::uint64_t& a, const std::uint64_t& b, std::uint64_t* result)
	{
		result[0] = _umul128(a, b, result + 1);
	}
#else
	void multiplication_unsigned(const std::uint64_t& a, const std::uint64_t& b, std::uint64_t* result)
	{
		const std::uint32_t* a1 = reinterpret_cast<const std::uint32_t*>(&a);
		const std::uint32_t* b1 = reinterpret_cast<const std::uint32_t*>(&b);
		const auto bits_h = sizeof(std::uint32_t) * 8;

		std::uint64_t c0 = static_cast<std::uint64_t>(a1[0])* b1[0];
		std::uint64_t c1 = static_cast<std::uint64_t>(a1[0])* b1[1] +
			static_cast<std::uint64_t>(a1[1])* b1[0];
		result[0] = c0 + ((c1 & 0xFFFFFFFF) << bits_h); // overflow here?
		result[1] = static_cast<std::uint64_t>(a1[1])* b1[1] + (c1 >> bits_h);
		result[1] += (result[0] < c0 ? 1 : 0); // overflow on result[0]!
	}
#endif
	void multiplication_unsigned(const std::uint32_t& a, const std::uint32_t& b, std::uint32_t* result)
	{
		std::uint64_t* res1 = reinterpret_cast< std::uint64_t*>(result);
		*res1 = static_cast<std::uint64_t>(a)* b;
	}

	void multiplication_signed(const std::int64_t& a, const std::int64_t& b, std::int64_t* result)
	{
		int sa = (a >= 0 ? 1 : -1);
		int sb = (b >= 0 ? 1 : -1);
		multiplication_unsigned(a * sa, b * sb, reinterpret_cast<std::uint64_t*>(result));
		if (sa * sb < 0)
		{
			result[0] = -result[0];
			if (result[0])
				result[1] = ~result[1];
			else
				result[1] = -result[1];
		}
	}

	std::uint64_t reminder_unsigned(const std::uint64_t* a, std::uint64_t b)
	{
		if (a[1] == 0)
			return a[0] % b;
		else
		{
			uint64_t c = a[1] % b;
			uint64_t a0 = a[0];
			std::uint64_t mask = std::uint64_t(1) << 63;
			for (int i = 0; i < 64; i++)
			{
				c = (c << 1) + (a0 & mask ? 1 : 0);
				if (c >= b)
					c -= b;
				a0 <<= 1;
			}
			return c;
		}
	}


	template <class T>
	class montgomery_t
	{
	public:
		typedef typename T signed_type;
		typedef typename std::make_unsigned<T>::type unsigned_type;
		montgomery_t(signed_type n = 65537)
		{
			init(n);
		}
		void init(signed_type n)
		{
			n_ = n;
			n1_ = montgomery_inverse(n);
			r2_ = compute_r2(n_);
		}
		unsigned_type project(signed_type a) const
		{
			return (a > 0 ? prod(r2_, a) : prod(r2_, a + n_));
		}
		unsigned_type mul(signed_type a, signed_type b) const 
		{
			if (a < 0)
				if (b < 0)
					return unred(prod(n_ - a, n_ - b));
				else
					return n_ - unred(prod(n_ - a, b));
			else
				if (b < 0)
					return n_ - unred(prod(a, n_ - b));
				else
					return unred(prod(a, b));
		}
		unsigned_type mul(unsigned_type a, unsigned_type b) const
		{
			return unred(prod(a, b));
		}
		unsigned_type mul_alt(unsigned_type a, unsigned_type b) const
		{
			return reduce(prod(project(a), project(b)));
		}
		unsigned_type power(unsigned_type base, unsigned_type exp) const
		{
			unsigned_type result1 = project(1);
			unsigned_type base1 = project(base);
			for (; exp; exp >>= 1)
			{
				if (exp & 1)
					result1 = prod(result1, base1);
				base1 = prod(base1, base1);
			}
			return prod(result1, 1);
		}
		unsigned_type reduce(const unsigned_type t) const
		{
			unsigned_type t1[2] = { t, 0 };
			return reduce(t1);
		}
		unsigned_type prod(const unsigned_type a1, const unsigned_type b1) const
		{
			unsigned_type t[2];
			multiplication_unsigned(a1, b1, t);
			return reduce(t);
		}
	private:
		unsigned_type unred(unsigned_type t) const // computes tR^-1 mod n
		{
			return prod(t, r2_);
		}
		unsigned_type reduce(const unsigned_type* t) const
		{
			unsigned_type mn[2];
			unsigned_type m = t[0] * n1_;
			multiplication_unsigned(m, n_, mn);
			if (mn[0] + t[0] < mn[0])
				mn[1]++;
			mn[1] += t[1];
			if (mn[1] > n_)
				mn[1] -= n_;
			return mn[1];
		}
		signed_type montgomery_inverse(const signed_type& n) const
		{
			signed_type q;
			int i, j;
			signed_type r0 = static_cast<signed_type>(1) << (sizeof(signed_type) * 8 - 1);
			signed_type r[2] = { r0, n };
			signed_type s[2] = { 1, 0 };
			signed_type t[2] = { 0, 1 };

			for (i = 0, j = 1; r[j]; i = j, j = (j + 1) & 1)
			{
				q = r[i] / r[j];
				r[i] = r[i] - q * r[j];
				s[i] = s[i] - q * s[j];
				t[i] = t[i] - q * t[j];
			}
			if (s[i] & 1)
			{
				t[i] += r0;
				s[i] = (s[i] + n) / 2;
			}
			else
				s[i] /= 2;
			return -r[i] * t[i];
		}
		unsigned_type compute_r2(const unsigned_type& n) const
		{
			unsigned_type c = static_cast<unsigned_type>(-1) % n;
			for (int i = 0; i < sizeof(unsigned_type) * 8; i++)
			{
				c = (c << 1) + 1;
				if (c >= n)
					c -= n;
			}
			return c + 1;
		}
		unsigned_type n_;
		unsigned_type n1_;
		unsigned_type r2_; //R^2 MOD n
	};



#ifdef HAVE_INTRINSIC
	class intrinsic_multiplier_t
	{
	public:
		intrinsic_multiplier_t(std::int64_t n = 65537) : n_(n) {}
		void init(std::int64_t n) { n_ = n; }
		std::uint64_t mul(std::uint64_t a, std::uint64_t b)
		{
			std::uint64_t h, l, r;
			l = _umul128(a, b, &h);
			_udiv128(h, l, n_, &r);
			return r;
		}
		std::uint64_t power(std::uint64_t base, std::uint64_t exp)
		{
			std::uint64_t result = 1;
			for (; exp; exp >>= 1)
			{
				if (exp & 1)
					result = mul(result, base);
				base = mul(base, base);
			}
			return result;
		}
	private:
		std::int64_t n_;
	};
#endif
}

#endif