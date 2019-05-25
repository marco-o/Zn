#ifndef quadratic_sieve_base_H
#define quadratic_sieve_base_H
namespace zn
{
	template <class large_int, class small_int, class real>
	class quadratic_sieve_base_t
	{
	public:
		enum smooth_status_e { smooth_idle_e, smooth_valid_e, smooth_candidate_e };
		struct base_t
		{
			std::vector<small_int>	prime_;
			std::vector<small_int>	residue_; // quadratic residue
			real					logp_;
			base_t(small_int p) : logp_(static_cast<real>(-std::log(p)))
			{
			}
			bool valid_for_polynomial(void) const { return prime_.size() > 1; }
			size_t powers(void) const { return prime_.size(); }
			small_int prime(size_t idx) const { return prime_[idx]; }
			small_int residue(size_t idx) const { return residue_[idx]; }
			real	  logp(void) const { return logp_; }
			bool eval_residue(const large_int &n)
			{
				large_int n1 = n % prime;
				residue_ = static_cast<small_int>(quadratic_residue<large_int>(n1, prime, prime / prime0)); // actually a power of prime
				return residue_ != 0;
			}
			static base_t build(small_int prime, const large_int &n)
			{
				base_t result(prime);
				small_int prime_pwr = prime;
				small_int prime1 = 1;
				small_int prime_power_end = static_cast<small_int>(std::numeric_limits<int>::max() / prime);
				for (int i = 0; (i < 10) && (prime_pwr < prime_power_end); i++)
				{
					small_int n1 = safe_cast<small_int>(n % prime_pwr);
					small_int residue = quadratic_residue<small_int>(n1, prime_pwr, prime1); // actually a power of prime
					if (residue == 0)
						break;
					result.prime_.push_back(prime_pwr);
					if (residue > prime_pwr / 2)
						residue = prime_pwr - residue;
					result.residue_.push_back(residue);
					// what about employ some kind of lifting?
					prime1 = prime_pwr;
					prime_pwr *= prime;
				}
				return result;
			}
		};
		struct base_ref_t : public base_t
		{
			base_ref_t(const base_t &base) : base_t(base) {}

			std::vector<int> smooths; // indexes
		};
		template <class smooth_t>
		void erase_base(base_ref_t &base, std::vector<smooth_t> &smooths)
		{
			size_t count = base.smooths.size();
			smooth_t &ref = smooths[base.smooths[0]];
			for (size_t i = 1; i < count; i++)
			{
				smooth_t &value = smooths[base.smooths[i]];
				value.compose(ref);
			}
			ref.invalidate();
		}
		template <class smooth_t>
		void erase_base(std::vector<base_ref_t> &base, std::vector<smooth_t> &smooths)
		{
			std::vector<int> base_remapping;
			int base_size = static_cast<int>(base.size());
			int base_map = 0;
#if DBG_SIEVE >= DBG_SIEVE_INFO
			std::cout << "Removing unused bases: start from " << base_size << " and " << smooths.size() << std::endl;
#endif
			for (int i = 0; i < base_size; i++)
			{
				auto &base_item = base[i];
				if (base_item.smooths.size() < 2)
				{
					base_remapping.push_back(-1);
					if (!base_item.smooths.empty())
						smooths[base_item.smooths[0]].invalidate();
					continue; // remove it
				}
				base_remapping.push_back(base_map);
				base_item.smooths.clear();
				if (base_map != i)
					base[base_map] = base_item;
				base_map++;
			}
			base.erase(base.begin() + base_map, base.end());
			// compact smooths
			int smooths_size = static_cast<int>(smooths.size());
			int smooths_map = 0;
			for (int i = 0; i < smooths_size; i++)
			{
				auto &smooth = smooths[i];
				if (!smooth.type() == smooth_valid_e)
					continue;
				smooth.remap_factors(base_remapping);
				if (i != smooths_map)
					smooths[smooths_map] = smooth;
				smooths_map++;
			}
			smooths.erase(smooths.begin() + smooths_map, smooths.end());
#if DBG_SIEVE >= DBG_SIEVE_INFO
			std::cout << "Base reduced to " << base.size() << ", smooth to " << smooths.size() << std::endl;
#endif
		}
		template <class smooth_t>
		void process_candidates_chunk(std::vector<base_ref_t> &base,
			std::map<large_int, smooth_t> &candidates,
			std::vector<smooth_t> &chunk,
			std::vector<smooth_t> &smooths,
			const large_int &n)
		{
			for (auto &s : chunk)
			{
				if (s.type() == smooth_candidate_e)
				{
					auto it = candidates.find(s.reminder());
					if (it != candidates.end())
					{
						s.compose(it->second, n, base);
#if 0 // DBG_SIEVE >= DBG_SIEVE_ERROR
						if (!s.invariant(n_, base_))
							std::cout << "Hmmm";
#endif // DBG_SIEVE	
					}
					else // just put it aside
					{
						candidates[s.reminder()] = s;
						continue;
					}
				}
				int si = static_cast<int>(smooths.size());
				for (auto f : s.factors())
					base[f].smooths.push_back(si);
				smooths.push_back(s);
			}
		}
		//
		// this function does an approximate reverse of estimation
		// of prime numbers pi(n) = n / log(n)
		//
		static small_int primes_range(small_int base_size)
		{
			double result = static_cast<double>(base_size);
			for (int i = 0; i < 10; i++)
				result = base_size * std::log(result);
			return static_cast<small_int>(result);
		}

	};
}
#endif
