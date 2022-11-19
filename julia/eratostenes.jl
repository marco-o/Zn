#
#
#

function power_module(base, exp, mod)
    exponent = exp
	power = base 
	result = convert(typeof(base), 1)
	while exponent != 0
	    if (exponent & 1) != 0
			result = result * power % mod
		end
		power = power * power % mod 
		exponent >>= 1
	end
	return result 
end


#
# Tonelli Shanks algorithm as described here
# https://en.wikipedia.org/wiki/Tonelli-Shanks_algorithm
# Parameter ps is required when p is actually a power of prime p^k 
# and in that case its value must be p^(k-1)
#
function quadratic_residue_odd(n, p, ps = 1)
	if p == 2
		return n % p
	end
	
	r = zero(typeof(p))
	q = p - ps 
	p2 = q ÷ 2 
	s = 0 
	while ((q & 1) == 0)
	    q >>= 1
		s += 1
	end
	
	# find non-quadratic residue z using Euler's criterion
    z = convert(typeof(p), 2) 
	while power_module(z, p2, p) == 1
	   z += 1 
	end
	   
	M = s 
	c = power_module(z, q, p) 
	t = power_module(n, q, p)
	R = power_module(n, (q+1) >> 1, p)
	while (t > 1) && (M > 0)
		t1 = t
		i = 1
		while i < M
			t1 = (t1 * t1) % p
			if t1 == 1
				c1 = 2 ^ (M - i - 1)
				b = power_module(c, c1, p)
				M = i 
				c = (b * b) % p
				t = (t * c) % p
				R = (R * b) % p
			end
			i += 1
		end
		if t1 != 1
			return 0
		end
	end
	if t == 1
	    return R
	else
		return 0
	end
end

function eratosthenes_sieve(limit)
	sieve=[true for n=0:limit]
	p = convert(typeof(limit), 2)
	result = Vector{typeof(limit)}()
	while (p * p <= limit)
		for j=p * p:p:limit
		   sieve[j] = false 
		end
		append!(result, p)
		for j=p+1:limit
			if (sieve[j])
			    p = j 
				break
			end
		end
	end
	for j=p:limit
		if (sieve[j])
			append!(result, j)
		end
	end
	return result 
end

function find_prime_base(N, limit)
	p = eratosthenes_sieve(limit) 
	result = Vector{typeof((limit, limit))}()
	for p1 ∈ p
		r = quadratic_residue_odd(N, p1) 
		if r != 0
			push!(result, (p1, r))
		end
	end
	return result
end

function prime_count(limit)
	return limit / log(limit) 
end

#
# Find out a number N such that N/logN = count
#
function inverse_prime_count(count)
	n = count * log(count)
	for i=1:10
		c0 = n / log(n)
		logn = log(n)
		dn = (count - c0) * logn ^ 2 / (logn - 1)
		n = n + dn
		if dn < 1
			break
		end
	end
	return trunc(Int, n)
end

function quadratic_sieve(N, base_size)
	sieve_limit = inverse_prime_count(base_size * 2)
	primes = find_prime_base(N, sieve_limit)
	a = 0
	return a
end

