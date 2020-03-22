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

int main(int argc, char* argv[])
{
	std::vector<std::int64_t> primes = zn::eratosthenes_sieve<std::int64_t>(2000);
	
	for (auto n : primes)
		if (n > 10)
		{
			montgomery_t<std::int64_t> mg(n);
			for (int i = 0; i < n; i++)
			{
				std::uint64_t p = std::rand() % n;
				std::uint64_t q = std::rand() % n;
				std::uint64_t b = mg.mul(p, q);
				std::uint64_t c = (p * q) % n;
				if (b != c)
					std::cout << b << std::endl;
			}
		}
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
#ifdef HAVE_INTRINSIC
	speed_test<intrinsic_multiplier_t, std::int64_t>("Intrinsic", p, q);
#endif
	return 0;
}