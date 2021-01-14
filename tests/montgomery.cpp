#include "znbasic.h"
#include "zneratosthenes_sieve.h"

#include <iostream>
#include <chrono>

#include "znmultiplier.h"

using namespace zn;


template <class T, class N>
void speed_test(const char *tag, N n, N b)
{
	auto start = std::chrono::high_resolution_clock::now();
	T multiplier(n);
	for (int i = 0; i < 10000; i++, b++)
	{
		N h = multiplier.power(b, n - 1);
		if (h != 1)
			std::cout << "Error\n";
	}
	auto time = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() ;
	std::cout << tag << ", time = " << time << std::endl;
}

template <class T>
int simple_test(const std::vector<std::int64_t> &primes)
{
	int result = 0;
	for (auto n : primes)
		if (n > 10)
		{
			montgomery_t<T> mg(static_cast<T>(n));
			for (int i = 0; i < n; i++)
			{
				T p = std::rand() % n;
				T q = std::rand() % n;
				T b = mg.mul(p, q);
				T c = (p * q) % n;
				if (b != c)
					result++;
			}
		}
	if (result)
		std::cout << result << std::endl;
	return result;
}

namespace zn
{
	template <class T>
	class montgomery_old_t
	{
	public:
		typedef typename T signed_type;
		typedef typename std::make_unsigned<T>::type unsigned_type;
		montgomery_old_t(signed_type n = 65537)
		{
			init(n);
		}
		void init(signed_type n)
		{
			n_ = n;
			n1_ = montgomery_inverse(n);
			r2_ = compute_r2(n_);
		}
		unsigned_type mul(signed_type a, signed_type b)
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
		unsigned_type mul(unsigned_type a, unsigned_type b)
		{
			return unred(prod(a, b));
		}
		unsigned_type power(unsigned_type base, unsigned_type exp)
		{
			unsigned_type result = 1;
			for (; exp; exp >>= 1)
			{
				if (exp & 1)
					result = mul(result, base);
				base = mul(base, base);
			}
			return result;
		}
	private:
		unsigned_type unred(unsigned_type t) // computes tR^-1 mod n
		{
			return prod(t, r2_);
		}
		unsigned_type prod(const unsigned_type a1, const unsigned_type b1)
		{
			unsigned_type t[2];
			multiplication_unsigned(a1, b1, t);
			return reduce(t);
		}
		unsigned_type reduce(const unsigned_type* t)
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
		signed_type montgomery_inverse(const signed_type& n)
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
		unsigned_type compute_r2(const unsigned_type& n)
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

}


template <class T, bool alt>
int trivial_test(const char *tag, const std::vector<std::int64_t>& primes)
{
	auto start = std::chrono::high_resolution_clock::now();
	int result = 0;
	for (auto n : primes)
		if (n > 10)
		{
			montgomery_t<T> mg(static_cast<T>(n));
			for (int i = 0; i < n; i++)
			{
				T p = std::rand() % n;
				T q = std::rand() % n;
				T b = (alt ? mg.mul_alt(p, q) : mg.mul(p, q));
				T b1 = p * q % n;
				if (b != b1)
					result++;
			}
		}
	auto time = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count();
	std::cout << tag << time << std::endl;
	if (result)
		std::cout << result << std::endl;
	return result;
}

int main(int argc, char* argv[])
{
	std::vector<std::int64_t> primes = zn::eratosthenes_sieve<std::int64_t>(4000);
	
	trivial_test<std::int64_t, false>("mul = ", primes);
	trivial_test<std::int64_t, true>("alt = ", primes);
	simple_test<std::int64_t>(primes);
	simple_test<std::int32_t>(primes);

	std::uint64_t p = 1269093085800313;
	std::uint64_t q = 3971584980549;
	montgomery_t<std::int64_t> mgp(p);
	for (int i = 0; i < 100; i++)
	{
		std::uint64_t t = i + q;
		auto p1 = mgp.power(t, p - 1);
		if (p1 != 1)
			std::cout << "Hmmm\n";
	}
	speed_test<montgomery_t<std::int64_t>, std::int64_t>("Montgomery ", p, q);
	speed_test<montgomery_old_t<std::int64_t>, std::int64_t>("Montgomery old ", p, q);
#ifdef HAVE_INTRINSIC
	speed_test<intrinsic_multiplier_t, std::int64_t>("Intrinsic", p, q);
#endif
	return 0;
}