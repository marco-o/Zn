//---------------------------------------------------------------------------------
//
//  Zn 
//  Copyright Marco Oman 2019
//
// Distributed under the Boost Software License, Version 1.0. 
// (See accompanying file LICENSE_1_0.txt or copy at 
// http://www.boost.org/LICENSE_1_0.txt)
//
#ifndef znquadratic_sieve_H
#define znquadratic_sieve_H

namespace zn
{
	//
	// this function does an approximate reverse of estimation
	// of prime numbers pi(n) = n / log(n)
	//
	template <class small_int>
	small_int primes_range(small_int base_size)
	{
		double result = base_size;
		for (int i = 0; i < 10; i++)
			result = base_size * std::log(result);
		return static_cast<int>(result);
	}

	template <class large_int, class small_int = int>
	large_int quadratic_sieve(const large_int &n, small_int base_size)
	{
		// double the range; half of them won't be a quadratic residue
		small_int range = primes_range(base_size * 2);
		auto primes = eratosthenes_sieve<small_int>(range);
		std::vector<std::pair<small_int, large_int>> base;
		large_int q;
		for (auto p : primes)
			if ((q = quadratic_residue<large_int>(n, p)))
				base.push_back(std::make_pair(p, q));
		std::cout << "base size = " << base.size() << std::endl;
		return n;
	}

}

#endif
