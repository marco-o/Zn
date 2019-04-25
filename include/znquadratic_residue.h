//---------------------------------------------------------------------------------
//
//  Zn 
//  Copyright Marco Oman 2019
//
// Distributed under the Boost Software License, Version 1.0. 
// (See accompanying file LICENSE_1_0.txt or copy at 
// http://www.boost.org/LICENSE_1_0.txt)
//
#ifndef znquadratic_residue_H
#define znquadratic_residue_H

namespace zn
{
	//
	// Tonelli Shanks algorithm as described here
	// https://en.wikipedia.org/wiki/Tonelli-Shanks_algorithm
	// Parameter ps is required when p is actually a power of prime p^k 
	// And in that case its value must be p^(k-1)
	//
	template <class Int>
	Int quadratic_residue_odd(Int n, Int p, Int ps = 1)
	{
		int s = 0;
		Int q = p - ps;
		Int p2 = q / 2;
		Int two = 2;

		if (p == 2)
			return n % p;
		for (Int r = 0 ; bit_test(q, 0) == 0 ; s++)
			divide_qr(q, two, q, r);

		// find quadratic non residue z using Euler's criterion
		Int z = 2;
		while (powm(z, p2, p) == 1) // if satisfied z *is* a quadratic residue
			z++;

		int M = s;
		Int c = powm(z, q, p);
		Int t = powm(n, q, p);
		Int R = powm(n, (q + 1) / 2, p);
		while (t > 1 && M > 0)
		{
			Int t1 = t;
			for (int i = 1; i < M; i++)
			{
				t1 = (t1 * t1) % p;
				if (t1 == 1)
				{
					Int c1 = power(two, M - i - 1);
					Int b = powm(c, c1, p);
					M = i;
					c = b * b % p; 
					t = (t * c) % p;
					R = (R * b) % p;
				}
			}
			if (t1 != 1)
				return 0;
		} 
		if (t == 1)
			return R;
		else 
			return 0;
	}

	template <class Int>
	Int quadratic_residue(Int n, Int p, Int ps = 1)
	{
		if (bit_test(p, 0))
			return quadratic_residue_odd(n, p, ps);
		else // job done only for odd n
		{
			if ((n & 7) != 1)
				return 0; 
			Int n1 = n / 2;
			for (Int r = 1; r < n1; r += 2)
				if ((r * r % p) == n)
					return r;
			return 0;
		}
	}
}

#endif
