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

	template <class large_int, class small_int>
	struct polynomial_t
	{
		std::vector<int> index; // base used for the polynomial
		large_int a0;
		large_int a;  // a = a0 * a0
		large_int b;
		large_int c;
		small_int x1;
		small_int x2;
		bool valid;
		template <class base_ref_t>
		large_int residue(const base_ref_t &base, const large_int &a, const large_int &a2, const large_int &n)
		{
			if (base.powers() > 1)
				return base.residue(1);
			else
				return quadratic_residue<large_int>(n, a2, a); // actually a power of prime
		}
		template <class base_ref_t>
		polynomial_t(const std::vector<int> &idx,
						const std::vector<base_ref_t> &base,
						const large_int &n) : index(idx), valid(true)
		{
			large_int r = 0;
			const base_ref_t &bp0 = base[idx[0]];
			a0 = bp0.prime(0);
			a = a0 * a0;
			b = residue(bp0, a0, a, n);
			valid = (b != 0);
			size_t index_count = index.size();
			for (size_t i = 1; valid && (i < index_count); i++)
			{
				const base_ref_t &bp = base[idx[i]];
				large_int p1 = bp.prime(0);
				a0 *= p1;
				large_int g = p1 * p1;
				large_int r1 = residue(bp, p1, g, n);
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

	template <class large_int, class small_int, class base_ref_t>
	class polynomial_generator_t
	{
		enum { first_base_e = 2 };
	public:
		struct base_t
		{
			small_int prime;
			float	logp;
			base_t(small_int p) : prime(p), logp(static_cast<float>(std::log(p))) {}
		};
		polynomial_generator_t(const large_int &n, small_int m, const std::vector<base_ref_t> &base)
		{
			large_int n1 = 2 * n;
			large_int a2 = safe_cast<large_int>(sqrt(n1)) / m;
			polynomial_seed_t result;
			target_ = real_op_t<float>::log1(a2) / 2;
			small_int largest = *base.rbegin();
			auto largest_log = real_op_t<float>::log1(largest);
			int min_order = static_cast<int>(std::ceil(target_ / largest_log));
			size_t size = base.size();
			for (size_t i = size / (min_order * 2 + 1) ; i < size ; i++)
				base_.push_back(base_t(base[i]));
			order_init(min_order);
		}
		polynomial_seed_t operator()(void)
		{
			increment();
			while (!within_target(index_))
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
		bool within_target(const polynomial_seed_t &seed)
		{
			return abs(seed.target_log - target_) < 0.1;
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
			auto it = base_.begin();
			for (size_t idx = index_.index.size() - 1; idx > 0 && !within_target(index_); idx--)
				if (index_.target_log < target_) // too light; increase first
				{
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
		std::vector<base_t> base_;
	};

}
#endif
