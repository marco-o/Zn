#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/miller_rabin.hpp>
#include <iostream>
#include <iomanip>
#include "znbasic.h"

namespace zn
{
	typedef boost::multiprecision::cpp_int large_int;

	bool custom_prime_test(const large_int &n)
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
		int seeds[] = { 2, 5, 19, 41 };
		for (auto x : seeds)
		{
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
				y = powm(y, 2, n);
			}
		}
		return true;  // Yeheh! probably prime.
	}

	bool rb_test(const boost::multiprecision::cpp_int &n)
	{
		bool result = !boost::multiprecision::miller_rabin_test(n, 4);
		if (result != custom_prime_test(n))
			custom_prime_test(n);
		
		return result;
	}
}
