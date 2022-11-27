#
#
#
using Printf
using OffsetArrays

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

#
# Informations about a prime packed together
#
struct PrimeInfo
	prime::Int64		# the actual prime
	residue::Int64	# quadratic residue of N mod p
end

function find_prime_base(N, limit)
	p = eratosthenes_sieve(limit) 
	result = Vector{PrimeInfo}()
	for p1 ∈ p
		r = quadratic_residue_odd(N, p1) 
		if r != 0
			push!(result, PrimeInfo(p1, r))
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
#
# All stuff required to sieve
#
struct PrimeInfoEx
	base::PrimeInfo
	idx1::Int64    # this is the smallest value that is a quadratic residue within sieving range
	idx2::Int64    # this is the largest value that is a quadratic residue 
	delta::Int64	 # The one corresponding to the other root															
	logp::UInt8		# log(p), used for sievin
end

struct SieveInfo
	primes::Vector{PrimeInfoEx}
	start::BigInt 
	size::Int64
	threshold::UInt8
	large_prime::Int64
end

function mylog(p) 
	return convert(UInt8, round(5.0 * log(p)))
end

function build_sieve_info(primes::Vector{PrimeInfo}, N, size)
	result = Vector{PrimeInfoEx}()
	start = convert(typeof(N), round(sqrt(N))) - (size ÷ 2) 
	for p ∈ primes
		st = p.prime * 2 - convert(Int64, start % p.prime)
		idx1 = (st + p.residue) % p.prime
		idx2 = (st - p.residue) % p.prime
		#i1 = (start + idx1) % p.prime
		#i2 = (start + idx2) % p.prime
		idx_min = min(idx1, idx2) 
		idx_max = max(idx1, idx2)
		logp = mylog(p.prime)
		info_ex = PrimeInfoEx(p, idx_min, idx_max, idx_max - idx_min, logp)
		push!(result, info_ex)
	end	
	last_prime = last(primes).prime 
	threshold = convert(UInt8, round(mylog(N) - 2 * mylog(last_prime)))
	return SieveInfo(result, start, size, threshold, last_prime*last_prime)
end

#
# Attempts the factorization of the polynomial y = x^2-N.
# Value i gets passed as is a valid hint in factorization
#
struct FactorInfo
	x
	y
	factors
end

function attempt_factorization(sieve_info::SieveInfo, N, i)
	x = sieve_info.start + i
	y = x * x - N
	y1 = y 
	factors = Vector{typeof((0, 0))}()
	@printf("Factorization of %d [%d]\n", y1, i)
	if y1 < 0
		push!(factors, (0, -1))
		y1 = -y1
	end
	for index ∈ eachindex(sieve_info.primes)
		p = sieve_info.primes[index]
		i1 = i % p.base.prime 
		if ((i1 == p.idx1) || (i1 == p.idx2))
			exponent = 0 
			while y1 % p.base.prime == 0
				exponent += 1
				y1 ÷= p.base.prime
			end
			if exponent % 2 == 1
				push!(factors, (index, p.base.prime))
			else
				if exponent == 0
					x1 = convert(Int64, x % p.base.prime)
					@printf("Error in division by %d [res = %d, x1 = %d]\n", p.base.prime, p.base.residue, x1)
				end
			end
		else
			if y1 % p.base.prime == 0
				@printf("Undetected division by %d\n", p.base.prime)
			end
		end
	end
	if (y1 < sieve_info.large_prime)
		return (FactorInfo(x, y, factors), convert(Int64, y1))
	end
	return (nothing, 0)
end

function sieve(sieve_info::SieveInfo, N)
	#
	# Phase 1: do actual sieving
	#
	base=OffsetVector([zero(UInt8) for n=1:sieve_info.size], 0:sieve_info.size-1)
	for p ∈ sieve_info.primes
		i = p.idx1 
		limit = sieve_info.size - p.delta 
		while i < limit
			base[i]         += p.logp
			base[i+p.delta] += p.logp 
			i += p.base.prime
		end
		if i < sieve_info.size
			base[i] += p.logp
		end
	end
	#
	# Phase 2: look for candidates for full factorization
	#
	result = Vector{FactorInfo}()
	reminders = Dict{Int64, FactorInfo}()
	for i = 0:sieve_info.size-1
		if base[i] > sieve_info.threshold
			(factorization, large_prime) = attempt_factorization(sieve_info, N, i)
			if large_prime != 0
				if large_prime == 1
					push!(result, factorization)
				else
					prev_factor = get(reminders, large_prime, nothing)
					if (isnothing(prev_factor))
						get!(reminders, large_prime, factorization)
					else
						push(result, merge_factors(factorization, prev_factor))
					end
				end
			end
		end
	end
	return result 
end

function quadratic_sieve(N, base_size, sieve_range)
	primes = find_prime_base(N, inverse_prime_count(base_size * 2))
	sieve_info = build_sieve_info(primes, N, sieve_range)
	factors = sieve(sieve_info, N)
	return factors
end

#
# TODO 
# Add -1 to base in factorization
# Handle large primes with a dictionary
#
quadratic_sieve(2309 * 3109, 20, 1000)
