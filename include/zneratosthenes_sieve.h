//---------------------------------------------------------------------------------
//
//  Zn 
//  Copyright Marco Oman 2019
//
// Distributed under the Boost Software License, Version 1.0. 
// (See accompanying file LICENSE_1_0.txt or copy at 
// http://www.boost.org/LICENSE_1_0.txt)
//
#ifndef eratosthenes_sieve_H
#define eratosthenes_sieve_H

#include <vector>

namespace zn
{
    template <class T, class N = int, class Slot = unsigned int>
    std::vector<T> eratosthenes_sieve(N bound)
    {
        const N bits_per_slot = sizeof(Slot) * 8 ;
        N slots = (bound + bits_per_slot - 1) / bits_per_slot ;
        Slot s0 = 0xAA ;
        for (N i = 0 ;i < sizeof(Slot) ; i++)
            s0 = s0 | (s0 << i * 8) ;
        std::vector<Slot> sieve(slots, s0) ;
        sieve[0] ^= 6 ; // remove 1, add 2
        for (N p = 3 ; p * p <= bound ; p += 2)
        {
            N slot = p / bits_per_slot ;
            N pos = p % bits_per_slot ;
            if (sieve[slot] & (1 << pos))
                for (N j = p * p ; j < bound ; j += p)
                {
                    slot = j / bits_per_slot ;
                    pos = j % bits_per_slot ;
                    sieve[slot] &= ~(1 << pos) ;
                }
        }
        // removing stuff past bound
        N slot = (bound - 1) / bits_per_slot ;
        N pos = (bound - 1) % bits_per_slot ;
        sieve[slot] &= static_cast<Slot>(-1) >> (bits_per_slot - pos - 1) ;
        // collect results
        std::vector<T> result ;
        for (N i = 0 ; i < slots ; i++)
            if (sieve[i])
                for (Slot j = 0, k = sieve[i] ; k ; j++, k >>= 1)
                    if (k & 1)
                        result.push_back(i * bits_per_slot + j) ;
        return result ;
    }
}

#endif