#
# Testing script
# Here some primes used for tests
# 79 223 557 12377 2309 3109 3461 5107 5857 6571 7297
#
using Printf
include("eratostenes.jl")

function test_quadratic_residue(limit)
    primes = eratosthenes_sieve(100)
    errors = 0
    for p in primes
        for n =1:p-1
            r = quadratic_residue_odd(n, p)
            # do the same 'by hand'
            r1 = 0
            for j=1:p÷2
                if j*j%p == n
                    r1 = j
                    break
                end
            end
            if (r1 != r) && (r1 != p-r)
                @printf("Error on %d: %d instead of %d\n", p, r, r1)
                errors += 1
            end
        end
    end
    return errors ;
end


#err = test_quadratic_residue(100)
#@printf("tested quadratic residue, found %d errors\n", err)

#inverse_prime_count(1000)

quadratic_sieve(2309 * 3109, 50)