#ifndef znpolynomialH
#define znpolynomialH

namespace zn
{

	struct polynomial_seed_t
	{
		std::vector<int> index; // base used for the polynomial
		float target_log;
		polynomial_seed_t(void) : target_log(-1) {}
		bool is_null(void) const { return target_log < -0.5; }
	};

	template <class small_int>
	struct prime_info_t
	{
		small_int prime0;
		small_int prime;
		small_int residue;
		prime_info_t(small_int p, small_int p2, small_int r) : prime0(p), prime(p2), residue(r) {}
	};

	template <class large_int, class small_int>
	struct polynomial_t
	{
		typedef small_int small_int_t;
		std::vector<int> index; // base used for the polynomial
		large_int a0;
		large_int a;  // a = a0 * a0
		large_int b;
		large_int c;
		small_int x1;
		small_int x2;
		bool valid;

		polynomial_t(const std::vector<int> &idx,
						const std::vector<prime_info_t<small_int>> &base,
						const large_int &n) : index(idx), valid(true)
		{
			large_int r = 0;
			const prime_info_t<small_int> &bp0 = base[idx[0]];
			a0 = bp0.prime;
			a = a0 * a0;
			b = quadratic_residue<large_int>(n, a, a0); // actually a power of prime
			valid = (b != 0);
			size_t index_count = index.size();
			for (size_t i = 1; valid && (i < index_count); i++)
			{
				const prime_info_t<small_int> &bp = base[idx[i]];
				large_int p1 = bp.prime0;
				a0 *= p1;
				large_int g = bp.prime;
				large_int r1 = bp.residue;
				if (r1 == 0)
					valid = false;
				large_int db = (r1 - b);
				auto ext = extended_euclidean_algorithm(a, g);
#if DBG_SIEVE >= DBG_SIEVE_TRACE
				if (std::get<0>(ext) != 1)
					throw not_relatively_prime_t<large_int>(db, a);
#endif
				large_int h = (std::get<1>(ext) * db) % g;
				b = b + h * a;
				a *= g;
				if (b > a / 2)
					b = a - b;
#if DBG_SIEVE >= DBG_SIEVE_TRACE
				if (((b * b) % a) != (n % a))
					throw std::runtime_error("Quadratic residue composition error");
				if (((n - b * b) % a) != 0)
					throw std::runtime_error("Quadratic residue internal error");
#endif
			}
			c = (b * b - n) / a;
			compute_zeros(n);
		}

		large_int eval(const large_int &x) const
		{
			large_int x1 = a * x + 2 * b;
			return x1 * x + c;
		}
		template <class real>
		real eval_log(small_int x) const
		{
#if 1
			large_int y = abs(eval(x));
			return real_op_t<real>::log1(y);
#else
			return loga + log(abs((x - z1) * (x - z2)) + 1);
#endif
		}
		std::pair<small_int, small_int> zeros(void) const
		{
			return std::make_pair(x1, x2);
		}
		// polynomial is negative in the range of the roots, including bounds
		void compute_zeros(const large_int &n)
		{
			large_int d = safe_cast<large_int>(sqrt(n));
			x1 = safe_cast<small_int>((-b - d) / a);
			x2 = safe_cast<small_int>((-b + d) / a);
#if DBG_SIEVE >= DBG_SIEVE_TRACE
			large_int y1 = eval(x1 - 1);
			large_int y0 = eval(x1);
			if (y1 < 0 || y0 > 0)
				throw std::runtime_error("Error in finding x1");
			y1 = eval(x2 + 1);
			y0 = eval(x2);
			if (y1 < 0 || y0 > 0)
				throw std::runtime_error("Error in finding x2");
#endif
		}
	};

	template <class large_int, class small_int>
	struct polynomial_siqs_t
	{
		typedef small_int small_int_t;
		std::vector<int> index; // base used for the polynomial
		large_int a0;
		large_int a;  // a = a0 * a0
		std::vector<large_int> b1;
		std::vector<large_int> c1;
		large_int b;
		large_int c;
		small_int x1;
		small_int x2;

		polynomial_siqs_t(const std::vector<int> &idx,
						  const std::vector<prime_info_t<small_int>> &primes, // product of g gives a
						  const large_int &n)
		{
			a = 1;
			a0 = 1;
			for (auto i : idx)
			{
				const prime_info_t<small_int> &p = primes[i];
				a  *= p.prime;
				a0 *= p.prime0;
			}
			std::vector<large_int> base;
			large_int bsum = 0;
			for (auto i : idx)
			{
				const prime_info_t<small_int> &p = primes[i];
				large_int ak = a / p.prime;
				small_int ak1 = safe_cast<small_int>(ak % p.prime);
				small_int ak2 = std::get<2>(extended_euclidean_algorithm(p.prime, ak1));
				if (ak2 < 0)
					ak2 += p.prime;
				ak = p.residue * ak * ak2;
				ak = ak % a;
				base.push_back(ak);
				bsum += ak;
			}
			size_t order = idx.size() - 1;
			size_t count = static_cast<size_t>(1) << order;
			size_t code = 0;
			std::vector<large_int> result;
			for (size_t i = 0; ; )
			{
				large_int bx = bsum % a;
				if (bit_test(bx, 0) == 0)
					if (bx > 0)
						bx -= a;
					else
						bx += a;
				b1.push_back(bx);
				c1.push_back((bx * bx - n) / a);
				if (++i == count)
					break;
				size_t next_code = i ^ (i >> 1);
				size_t changed_bit = next_code ^ code;
				size_t index = lsb(changed_bit);
				if (changed_bit & next_code)
					bsum -= 2 * base[index];
				else
					bsum += 2 * base[index];
				code = next_code;
			}
			compute_zeros(n);
		}
		size_t count(void) const { return b1.size(); }
		void select(size_t index)
		{
			b = b1[index];
			c = c1[index];
		}

		large_int eval(const large_int &x) const
		{
			large_int x1 = a * x + 2 * b;
			return x1 * x + c;
		}
		template <class real>
		real eval_log(small_int x) const
		{
#if 1
			large_int y = abs(eval(x));
			return real_op_t<real>::log1(y);
#else
			return loga + log(abs((x - z1) * (x - z2)) + 1);
#endif
		}
		std::pair<small_int, small_int> zeros(void) const
		{
			return std::make_pair(x1, x2);
		}
		// polynomial is negative in the range of the roots, including bounds
		void compute_zeros(const large_int &n)
		{
			large_int d = safe_cast<large_int>(sqrt(n));
			//x1 = safe_cast<small_int>((-b - d) / a);
			//x2 = safe_cast<small_int>((-b + d) / a);
			x1 = safe_cast<small_int>(- d / a);
			x2 = safe_cast<small_int>(  d / a);
		}
	};

	template <class large_int, class small_int>
	class polynomial_generator_t
	{
		enum { first_base_e = 2 };
	public:
		struct base_t : prime_info_t<small_int>
		{
			float	logp;
			base_t(const prime_info_t<small_int> &p) : prime_info_t<small_int>(p), 
				                     logp(static_cast<float>(std::log(p.prime0))) {}
		};
		polynomial_generator_t(const large_int &n, small_int m, 
			                   const std::vector<prime_info_t<small_int>> &base)
		{
			large_int n1 = 2 * n;
			large_int a2 = safe_cast<large_int>(sqrt(n1)) / m;
			polynomial_seed_t result;
			target_ = real_op_t<float>::log1(a2) / 2;
			/*
		    the point is keep enough values with log around target / min_order
			order_init could take care of that
			*/
			small_int largest = base.rbegin()->prime;
			auto largest_log = real_op_t<float>::log1(largest);
			int min_order = static_cast<int>(std::ceil(target_ / largest_log));
			size_t size = base.size();
			for (size_t i = 0 ; i < size ; i++)
				base_.push_back(base_t(base[i]));
			order_init(min_order);
		}
		polynomial_seed_t operator()(void)
		{
			increment();
			while (!within_target(index_, true))
				increment();
			return index_;
		}
		void order_init(int order)
		{
			index_.index.clear();
			for (int i = 1; i < order; i++)
				index_.index.push_back(i - 1);
			index_.index.push_back(static_cast<int>(base_.size() - 1));
			coarse_init();
		}
	private:
		bool within_target(const polynomial_seed_t &seed, bool update_err)
		{
			float err = abs(seed.target_log - target_);
			if (update_err)
			{
				average_count_++;
				const float alpha = 0.1f + 1.0f / (average_count_ + 0.2f);
				average_target_error_ = average_target_error_ * (1.0f - alpha) + alpha * err;
			}
			return err < average_target_error_ ; // this means keep approximately half of them
		}
		void increment(void)
		{
			size_t order = index_.index.size();
			if (index_.target_log < target_) // short of target, increase first element
			{
				size_t idx = order - 2;
				index_.target_log -= base_[index_.index[idx]].logp;
				index_.index[idx]++;
				index_.target_log += base_[index_.index[idx]].logp;
			}
			else
			{
				size_t idx = order - 1;
				index_.target_log -= base_[index_.index[idx]].logp;
				index_.index[idx]--;
				index_.target_log += base_[index_.index[idx]].logp;
			}
			auto rit = index_.index.rbegin();
			if (rit[0] != rit[1])
				return;

			if (order > 2)
			{
				size_t i = order - 2;
				for (; i > 0; i--)
					if (index_.index[i] + 1 < index_.index[i + 1])
						break;
				index_.index[i]++;
				for (i++; i < order - 1; i++)
					index_.index[i] = index_.index[i - 1] + 1;
				index_.index[order - 1] = static_cast<int>(base_.size() - 1);
			}
			// promote to a larger order
			if (index_.index[order - 2] >= index_.index[order - 1])
			{
				for (size_t i = 0; i < order; i++)
					index_.index[i] = static_cast<int>(i + first_base_e);
				index_.index.push_back(static_cast<int>(base_.size() - 1));
			}
			coarse_init();
		}
		void coarse_init(void)
		{
			init_seed(index_);
			int limit_mask = 3; // exit when both limits have been touched
			auto it = base_.begin();
			for (size_t idx = index_.index.size() - 1; idx > 0 && limit_mask && !within_target(index_, false) ; idx--)
				if (index_.target_log < target_) // too light; increase first
				{
					limit_mask &= 2; // remove lower bit
					index_.target_log -= base_[index_.index[idx - 1]].logp;
					it = std::lower_bound(base_.begin() + index_.index[idx - 1] + 1,
						base_.begin() + index_.index[idx] - 1,
						target_ - index_.target_log,
						[](const base_t &base, float p) {
						return base.logp < p;
					});
					index_.index[idx - 1] = static_cast<int>(it - base_.begin());
					index_.target_log += base_[index_.index[idx - 1]].logp;
				}
				else // too heavy, lower the high
				{
					limit_mask &= 1; // remove higher bit
					index_.target_log -= base_[index_.index[idx]].logp;
					it = std::upper_bound(base_.begin() + index_.index[idx - 1] + 1,
						base_.begin() + index_.index[idx] - 1,
						target_ - index_.target_log,
						[](float p, const base_t &base) {
						return p < base.logp;
					});
					index_.index[idx] = static_cast<int>(it - base_.begin());
					index_.target_log += base_[index_.index[idx]].logp;
				}
		}
		void init_seed(polynomial_seed_t &seed)
		{
			seed.target_log = 0;
			size_t size = seed.index.size();
			for (size_t i = 0; i < size; i++)
				seed.target_log += base_[seed.index[i]].logp;
		}
		polynomial_seed_t	index_;
		float				target_;
		float				average_target_error_ = 0.1f ;
		int					average_count_ = 0;
		std::vector<base_t> base_;
	};

	template <class poly_t, class real_t>
	class sieve_range_t
	{
		typedef typename poly_t::small_int_t small_int;
		const poly_t   &poly_;
		small_int		m_;
	public:
		typedef long run_int_t;
		struct sieve_run_t
		{
			run_int_t p; // may be a power of prime
			run_int_t a; // poly.a % prime
			run_int_t r;
			run_int_t a1; // inverse of a
			real_t    lg; // logarithm to subtract
			// the following values change for each sub-polynomial
			run_int_t b;  // poly.b % prime
			run_int_t x;
			int bix;
		};
		sieve_range_t(const poly_t &poly, small_int m) : poly_(poly), m_(m) {}
		std::vector<real_t> fill(real_t offset = 4) const 
		{
			std::vector<real_t> values(2 * static_cast<size_t>(m_));
			auto zeros = poly_.zeros();
			if (-m_ < zeros.first)
			{
				real_t u = poly_.eval_log<real_t>(zeros.first);
				real_t v = poly_.eval_log<real_t>(zeros.second);
				fill_range(values, -m_, poly_.eval_log<real_t>(m_), zeros.first, u);
				fill_range(values, zeros.first, u, zeros.second, v);
				fill_range(values, zeros.second, v, m_, poly_.eval_log<real_t>(m_));
			}
			else // enters here only during tests
				fill_range(values, -m_, poly_.eval_log<real_t>(-m_), m_, poly_.eval_log<real_t>(m_));
			for (auto &v : values)
				v += offset;
			return values;
		};
		template <class base_t>
		static void build_run(const poly_t &poly,
								const std::vector<base_t> &bases,
								std::vector<sieve_run_t> &runs)
		{
			sieve_run_t run;
			size_t bases_size = bases.size();
			for (size_t k = 0 ; k < bases_size ; k++)
			{
				auto &base = bases[k];
				size_t powers = base.powers();
				run.lg = base.logp();
				run.bix = k;
				for (size_t i = 0; i < powers; i++)
				{
					run.p = static_cast<run_int_t>(base.prime(i));
					run.a = safe_cast<run_int_t>(poly.a % run.p);
					if (run.a == 0)
						break; // we are hitting a divisor of a
					run.a1 = std::get<1>(extended_euclidean_algorithm<run_int_t>(run.a, run.p));
					run.b = safe_cast<run_int_t>(poly.b % run.p);
					run.r = static_cast<run_int_t>(base.residue(i));
					run.x = ((run.p - run.b + run.r) * run.a1) % run.p;
					runs.push_back(run);
					if (run.p > 2)
					{
						run.r = run.p - run.r;
						run.x = ((run.p - run.b + run.r) * run.a1) % run.p;
						runs.push_back(run);
					}
					run.bix = -1;
				}
			}
		}
		static void update_run(const poly_t &poly,
							   std::vector<sieve_run_t> &runs)
		{
			for (auto &run : runs)
			{
				run.b = safe_cast<run_int_t>(poly.b % run.p);
				run.x = ((run.p - run.b + run.r) * run.a1) % run.p;
			}
		}
	private:
		void fill_range(std::vector<real_t> &values,
					    small_int begin,	 real_t t_1,
						small_int end,	     real_t t1) const 
		{
			small_int mid = (begin + end) / 2;
			real_t t0 = poly_.eval_log<real_t>(mid);
			auto q1 = t_1 / static_cast<double>(t0);
			auto q0 = static_cast<double>(t0) / t1;
			auto t = q1 / q0 + q0 / q1 - 2;
			if (t < 0.002)
			{
				fill_linear(values, begin, t_1, mid, t0);
				fill_linear(values, mid, t0, end, t1);
			}
			else
				if (end - begin < 16)
					fill_exact(values, begin, end);
				else
				{
					fill_range(values, begin, t_1, mid, t0);
					fill_range(values, mid, t0, end, t1);
				}
		}
		void fill_linear(std::vector<real_t> &values,
						 const small_int &x1, real_t t1,
						 const small_int &x2, real_t t2) const
		{
			size_t size = static_cast<size_t>(x2 - x1);
			size_t offset = static_cast<size_t>(x1 + m_);
			double mx = static_cast<double>(t2 - t1) / size;
			for (size_t i = 0; i < size; i++)
				values[i + offset] = static_cast<real_t>(t1 + mx * i);
		}
		void fill_exact(std::vector<real_t> &values,
			            small_int begin, small_int end) const
		{
			size_t size = safe_cast<size_t>(end - begin);
			size_t offset = safe_cast<size_t>(begin + m_);
			for (size_t i = 0; i < size; i++)
				values[offset + i] = poly_.eval_log<real_t>(begin + i);
		}
	};


}
#endif
