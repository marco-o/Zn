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
	return std::tuple<T, T>(r[i] * s[i], r[i] * t[i]);
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


int main(int argc, char* argv[])
{
	std::vector<std::int64_t> primes = zn::eratosthenes_sieve<std::int64_t>(200);
	
	for (auto n : primes)
		if (n > 2)
		{
			auto m = montgomery_inverse(n);
			std::int64_t n1 = std::get<1>(m);
			std::int64_t r1 = static_cast<std::uint64_t>(std::get<0>(m));
			std::uint64_t res[2];
			multiplication_unsigned(n, n1, res);
			if (res[0] != 1)
				std::cout << "Hmmm\n";
		}
	return 0;
}