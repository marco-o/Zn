#include "znbasic.h"
#include "zneratosthenes_sieve.h"

#include <iostream>

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
        q = r[i] / r[j] ;
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

void multiplication_unsigned(const std::uint64_t& a, const std::uint64_t& b, std::uint64_t*result)
{
	const std::uint32_t* a1 = reinterpret_cast<const std::uint32_t*>(&a);
	const std::uint32_t* b1 = reinterpret_cast<const std::uint32_t*>(&b);
	const auto bits_h = sizeof(std::uint32_t) * 8;

	if ((a1[1] == 0) && (b1[1] == 0))
	{
		result[0] = static_cast<std::uint64_t>(a1[0])* b1[0];
		result[1] = 0;
	}
	else
	{
		std::uint64_t c0 = static_cast<std::uint64_t>(a1[0]) * b1[0];
		std::uint64_t c1 = static_cast<std::uint64_t>(a1[0]) * b1[1] + 
			               static_cast<std::uint64_t>(a1[1]) * b1[0];
		result[0] = c0 + ((c1 & 0xFFFFFFFF) << bits_h); // overflow here?
		result[1] = static_cast<std::uint64_t>(a1[1])* b1[1] + (c1 >> bits_h);
		if (result[0] < c0) // overflow on result[0]!
			result[1]++;
	}
}

void multiplication_signed(const std::int64_t& a, const std::int64_t& b, std::int64_t* result)
{
	int sa = (a >= 0 ? 1 : -1);
	int sb = (b >= 0 ? 1 : -1);
	multiplication_unsigned(a * sa, b * sb, reinterpret_cast<std::uint64_t *>(result));
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



// assumes r = 1 << 64
class montgomery_t
{
public:
	montgomery_t(std::int64_t n) : n_(n)
	{
		n1_ = montgomery_inverse(n);
		r2_ = compute_r2(n_);
	}
	std::uint64_t unred(std::uint64_t t) // computes tR^-1 mod n
	{
		return prod(t, r2_);
	}
	std::uint64_t prod(const std::uint64_t a1, const std::uint64_t b1)
	{
		std::uint64_t t[2];
		multiplication_unsigned(a1, b1, t);
		return reduce(t);
	}
	std::uint64_t reduce(const std::uint64_t* t)
	{
		std::uint64_t mn[2];
		std::uint64_t m = t[0] * n1_;
		multiplication_unsigned(m, n_, mn);
		if (mn[0] + t[0] < mn[0])
			mn[1]++;
		mn[1] += t[1];
		if (mn[1] > n_)
			mn[1] -= n_;
		return mn[1];
	}
private:
	template <class T>
	T montgomery_inverse(const T& n)
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
		return -r[i] * t[i];
	}
	template <class T>
	T compute_r2(const T &n)
	{
		uint64_t c = static_cast<T>(-1) % n;
		for (int i = 0; i < sizeof(T) * 8; i++)
		{
			c = (c << 1) + 1;
			if (c >= n)
				c -= n;
		}
		return c + 1;
	}
	std::uint64_t n_;
	std::uint64_t n1_ ;
	std::uint64_t r2_; //R^2 MOD n
};

int main(int argc, char* argv[])
{
	std::vector<std::int64_t> primes = zn::eratosthenes_sieve<std::int64_t>(2000);
	
	for (auto n : primes)
		if (n > 10)
		{
			montgomery_t mg(n);
			for (int i = 0; i < n; i++)
			{
				std::uint64_t p = std::rand() % n;
				std::uint64_t q = std::rand() % n;
				std::uint64_t a = mg.prod(p, q);
				std::uint64_t b = mg.unred(a);
				std::uint64_t c = (p * q) % n;
				if (b != c)
					std::cout << b << std::endl;
			}
		}
	return 0;
}