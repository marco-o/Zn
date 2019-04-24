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
	//
	template <class Int>
	Int quadratic_residue(Int n, Int p)
	{
		int s = 0;
		Int q = p - 1;
		Int p2 = q / 2;
		Int two = 2;

		if (p == 2)
			return n % p;
		for (Int r = 0 ; bit_test(q, 0) == 0 ; s++)
			divide_qr(q, two, q, r);

		// find quadratic non residue z using eulere criterion
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
}

#endif
