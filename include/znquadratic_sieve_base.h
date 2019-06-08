#ifndef quadratic_sieve_base_H
#define quadratic_sieve_base_H
namespace zn
{
	template <class real>
	struct real_op_t
	{
		static real unit(void) { return 1; }
		template <class T>
		static real log1(T value)
		{
			return static_cast<real>(log(abs(safe_cast<double>(value))));
		}
	};

	template <>
	struct real_op_t<short>
	{
		static short unit(void) { return 120; }
		template <class T>
		static short log1(T value)
		{
			return static_cast<short>(unit() * safe_cast<double>(log(abs(safe_cast<double>(value)))));
		}
	};

	template <class large_int, class small_int>
	small_int premultiplier(const large_int &n, const std::vector<small_int> &primes)
	{
		int next_square_root = 2;
		double max_value = 0;
		int k_max = 3;
		double clog2 = std::log(2);
		for (int k = 3; k < 100; k++)
		{
			if (k == next_square_root * next_square_root)
			{
				next_square_root++;
				continue;
			}
			large_int kn = k * n;
			double fk = -std::log(k) / 2;
			for (auto p : primes)
				if (p == 2)
				{
					switch (safe_cast<int>(kn % 8))
					{
					case 1:
						fk += 2 * clog2;
						break;
					case 2:
					case 3:
					case 6:
					case 7:
						fk += clog2 / 2;
						break;
					case 5:
						fk += clog2;
						break;
					default:
						break;
					}
				}
				else
				{
					small_int knp = safe_cast<small_int>(kn % p);
					if (quadratic_residue(knp, p) > 0)
						if (k % p == 0)
							fk += std::log(p) / p;
						else
							fk += 2 * std::log(p) / (p - 1);
					if (p > 1000)
						break;
				}
			if (fk > max_value)
			{
				max_value = fk;
				k_max = k;
			}
		}
		return k_max;
	}

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
			base_t(small_int p) : logp_(-real_op_t<real>::log1(p))
			{
			}
			bool valid_for_polynomial(void) const { return prime_.size() > 1; }
			size_t powers(void) const { return prime_.size(); }
			small_int prime(size_t idx) const { return prime_[idx]; }
			small_int residue(size_t idx) const { return residue_[idx]; }
			real	  logp(void) const { return logp_; }
			static base_t build(small_int prime, const large_int &n)
			{
				base_t result(prime);
				small_int prime_pwr = prime;
				small_int prime1 = 1;
				small_int prime_power_end = static_cast<small_int>(std::numeric_limits<int>::max());
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
		void erase_smooth_ref(std::vector<base_ref_t> &base, const std::vector<int> &base_indexes, int smooth_index)
		{
			for (auto bindex : base_indexes)
			{
				auto &erased_base = base[bindex];
				auto jt = std::lower_bound(erased_base.smooths.begin(), erased_base.smooths.end(), smooth_index);
				if (jt != erased_base.smooths.end())
					erased_base.smooths.erase(jt);
				else
					LOG_INFO << "Cannot erase..." << log_base_t::newline_t();
			}
		}
		void insert_smooth_ref(std::vector<base_ref_t> &base, const std::vector<int> &base_indexes, int smooth_index)
		{
			for (auto bindex : base_indexes)
			{
				auto &added_base = base[bindex];
				auto jt = std::lower_bound(added_base.smooths.begin(), added_base.smooths.end(), smooth_index);
				added_base.smooths.insert(jt, smooth_index);
			}
		}
		template <class smooth_t>
		void erase_base_items(std::vector<base_ref_t> &base, std::vector<smooth_t> &smooths, const large_int &n)
		{
			// erase bases with less references (smooths)
			size_t size = base.size();
			std::vector<std::pair<size_t, size_t>> base_index(size);
			for (size_t i = 0; i < size; i++)
				base_index[i] = std::make_pair(base[i].smooths.size(), i);
			std::sort(base_index.begin(), base_index.end());
			int reduce = static_cast<int>(base.size() / 5);

			for (int k = 0 ; k < reduce; k++)
			{
				int index = static_cast<int>(size - k - 1); // base_index[k].second);
				base_ref_t &item = base[index];
				if (item.smooths.size() < 2)
				{
					if (item.smooths.size() == 1)
					{
						smooth_t &smooth = smooths[item.smooths[0]];
						erase_smooth_ref(base, smooth.factors(), item.smooths[0]);
						smooth.invalidate();
					}
					continue;
				}
				smooth_t &smooth = smooths[item.smooths[0]];
				size_t count = item.smooths.size();
				for (size_t i = count - 1 ; i > 0 ; i--)
				{
					int smooth_index = item.smooths[i];
					smooth_t &value = smooths[smooth_index];
					if (smooth.type() != smooth_valid_e)
						continue;
					std::vector<int> erased;
					std::vector<int> added;
					value.compose(smooth, n, base, &erased, &added);
					erase_smooth_ref(base, erased, smooth_index);
					insert_smooth_ref(base, added, smooth_index);
				}
				// remove all rferences to this smooth
				erase_smooth_ref(base, smooth.factors(), item.smooths[0]);
				smooth.invalidate();

				//item.smooths.clear();
			}
			LOG_INFO << "Removed all bases up to " << (base.size() - reduce) << log_base_t::newline_t() ;
		}
		static size_t actual_base_size(std::vector<base_ref_t> &base)
		{
			size_t result = 0;
			for (const auto &item : base)
				if (!item.smooths.empty())
					result++;
			return result;
		}
		template <class smooth_t>
		void erase_base(std::vector<base_ref_t> &base, std::vector<smooth_t> &smooths, const large_int &n)
		{
			std::vector<int> base_remapping;
			int base_size = static_cast<int>(base.size());
			int base_map = 0;
			LOG_INFO << "Removing unused bases: start from " << base_size 
				     << " and " << smooths.size() << log_base_t::newline_t();
			erase_base_items(base, smooths, n);
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
				if (smooth.type() != smooth_valid_e)
					continue;
				smooth.remap_factors(base_remapping);
				if (i != smooths_map)
					smooths[smooths_map] = smooth;
				smooths_map++;
			}
			smooths.erase(smooths.begin() + smooths_map, smooths.end());
			LOG_INFO << "Base reduced to " << base.size() << ", smooth to " << smooths.size() << log_base_t::newline_t();
		}
		template <class smooth_t>
		int  process_candidates_chunk(std::vector<base_ref_t> &base,
									  std::map<large_int, smooth_t> &candidates,
									  std::vector<smooth_t> &chunk,
									  std::vector<smooth_t> &smooths,
									  const large_int &n)
		{
			int promoted = 0;
			for (auto &s : chunk)
			{
				if (s.type() == smooth_candidate_e)
				{
					auto it = candidates.find(s.reminder());
					if (it != candidates.end())
					{
#ifdef HAVE_CANDIDATE_ANALYSYS
						size_t value = static_cast<size_t>(std::round(log(safe_cast<double>(s.reminder()))));
						if (value >= composed_hist_.size())
							composed_hist_.resize(value + 1);
						composed_hist_[value]++;
#endif
						promoted++;
						s.compose(it->second, n, base);
#if DBG_SIEVE >= DBG_SIEVE_ERROR
						if (!s.invariant(n, base))
							std::cout << "Hmmm";
#endif // DBG_SIEVE	
					}
					else // just put it aside
					{
#ifdef HAVE_CANDIDATE_ANALYSYS
						size_t value = static_cast<size_t>(std::round(log(safe_cast<double>(s.reminder()))));
						if (value >= candidate_hist_.size())
							candidate_hist_.resize(value + 1);
						candidate_hist_[value]++;
#endif
						candidates[s.reminder()] = s;
						continue;
					}
				}
				int si = static_cast<int>(smooths.size());
				for (auto f : s.factors())
					base[f].smooths.push_back(si);
				smooths.push_back(s);
			}
			return promoted;
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
#ifdef HAVE_CANDIDATE_ANALYSYS
		void print_analysis(small_int largest_prime)
		{
			size_t size = candidate_hist_.size();
			if (composed_hist_.size() < size)
				composed_hist_.resize(size);
			size_t begin = static_cast<size_t>(round(log(largest_prime))) - 1;
			int candidates = std::accumulate(candidate_hist_.begin(), candidate_hist_.end(), 0);
			int composed = std::accumulate(composed_hist_.begin(), composed_hist_.end(), 0);
			int candidates1 = candidates;
			int composed1 = composed;
			std::cout << "Composed = " << composed << std::endl;
			for (size_t i = begin; i < size; i++)
			{
				std::cout << candidate_hist_[i] << "  \t" 
					      << (candidates1 * 100) / candidates << " \t"					       
					      << composed_hist_[i] << "  \t" 
						  << (composed1 * 100) / composed << std::endl ;
				composed1 -= composed_hist_[i];
				candidates1 -= candidate_hist_[i];
			}
		}
		std::vector<int> candidate_hist_;
		std::vector<int> composed_hist_;
#endif
	};

}
#endif
