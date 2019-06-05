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
			int reduce = static_cast<int>(base.size() / 5);
			for (int k = 0 ; k < reduce; k++)
			{
				int index = static_cast<int>(base.size() - k - 1);
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
#ifdef HAVE_CANDIDATE_ANALYSYS
						size_t value = static_cast<size_t>(std::round(log(safe_cast<double>(s.reminder()))));
						if (value >= composed_hist_.size())
							composed_hist_.resize(value + 1);
						composed_hist_[value]++;
#endif
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

	class linear_solver_t
	{
	public:
		typedef unsigned int slot_t;
		enum { bits_per_slot = 8 * sizeof(slot_t) };
		// linear system: rows is #base_, cols is number of smooths
		template <class smooth_t>
		std::vector<std::vector<int>> solve(std::vector<smooth_t> &smooths, size_t base_size)
		{
			size_t smooth_size = smooths.size();
			size_t slots = (smooth_size + bits_per_slot - 1) / bits_per_slot;
			std::vector<std::vector<slot_t>> matrix(base_size, std::vector<slot_t>(slots, 0));
			for (size_t i = 0; i < smooth_size; i++)
			{
				const smooth_t &smooth = smooths[i];
				int slot = static_cast<int>(i / bits_per_slot);
				slot_t mask = 1 << static_cast<slot_t>(i % bits_per_slot);
				if (smooth.sign_neg())
					matrix[0][slot] |= mask;
				for (const auto idx : smooth.factors())
					matrix[idx + 1][slot] |= mask;
			}
#if DBG_SIEVE >= DBG_SIEVE_DEBUG
			std::cout << "Initial status\n";
			trace_matrix(matrix, smooth_size);
#endif
			std::vector<int> smooth_perm(smooth_size);
			for (size_t i = 0; i < smooth_size; i++)
				smooth_perm[i] = static_cast<int>(i);
			int i = 0; // actual row used as pivot
			std::vector<int> rows;
			std::vector<int> rows_idx;
			for (size_t k = 0; k < base_size; k++)
			{
				rows_idx.push_back(static_cast<int>(i));
#if DBG_SIEVE >= DBG_SIEVE_DEBUG
				std::cout << "\nPass " << i << " out of " << base_size << "\n";
				trace_matrix(matrix, smooth_size);
#endif
				int	   slot = i / bits_per_slot;
				slot_t mask = 1 << static_cast<slot_t>(i % bits_per_slot);
				auto &v = matrix[k];
				int j = find_nonzero(v, static_cast<int>(smooth_size), i);
				if (j < 0)
					continue;
				if (j != i)
				{
					swap_bit(matrix, i, j);
					std::swap(smooth_perm[i], smooth_perm[j]);
#if DBG_SIEVE >= DBG_SIEVE_DEBUG
					std::cout << "Swap columns " << i << ", " << j << std::endl;
#endif
				}
				for (j = static_cast<int>(k + 1); j < static_cast<int>(base_size); j++)
					if (matrix[j][slot] & mask)
						for (size_t h = slot; h < slots; h++)
							matrix[j][h] ^= v[h];
				rows.push_back(static_cast<int>(k));
				i++;
			}
			// backsubstitution
			size_t rows_size = rows.size();
			for (size_t c = rows_size - 1; c > 0; c--) // indexing on column, which is != from row, i.e. 'i'
			{
				int i = rows[c];
#if DBG_SIEVE >= DBG_SIEVE_DEBUG
				std::cout << "\nBack " << i << " out of " << base_size << "\n";
				trace_matrix(matrix, smooth_size);
#endif
				auto &v = matrix[i];
				int	   slot = static_cast<int>(c / bits_per_slot);
				slot_t mask = 1 << static_cast<slot_t>(c % bits_per_slot);
				if (v[slot] & mask) // may be zero...
					for (int j = i - 1; j >= 0; j--)
						if (matrix[j][slot] & mask)
						{
							auto &w = matrix[j];
							for (size_t k = slot; k < slots; k++)
								w[k] ^= v[k];
						}
			}

			std::vector<std::vector<int>> result;
#if DBG_SIEVE >= DBG_SIEVE_DEBUG
			std::cout << "Result:\n";
			trace_matrix(matrix, smooth_size);
#endif
			for (size_t i = rows_size; i < smooth_size; i++)
			{
				std::vector<int> idx;
				slot_t slot = static_cast<slot_t>(i / bits_per_slot);
				slot_t mask = 1 << static_cast<slot_t>(i % bits_per_slot);

				for (size_t j = 0; j < base_size; j++)
					if (matrix[j][slot] & mask)
						idx.push_back(smooth_perm[rows_idx[j]]);  // make use of smooths_perm
				idx.push_back(smooth_perm[i]);
				result.push_back(idx);
			}
			return result;
		}
	private:
		int find_nonzero(const std::vector<slot_t> &v, int size, int i)
		{
			int slot = static_cast<int>(i / bits_per_slot);
			int j = static_cast<int>(i % bits_per_slot);
			slot_t mask = 1 << j;
			for (i = slot * bits_per_slot; i < size; i += bits_per_slot, slot++)
			{
				slot_t s = v[slot];
				if (s)
					for (; j < bits_per_slot; j++, mask <<= 1)
						if (s & mask)
							return i + j;
				j = 0;
				mask = 1;
			}
			return -1;
		}
		void swap_bit(std::vector<std::vector<slot_t>> &matrix, int i, int j) // assume j > i
		{
			const int sloti = i / bits_per_slot;
			const int slotj = j / bits_per_slot;
			const size_t size = matrix.size();
			const slot_t di = static_cast<slot_t>(i % bits_per_slot);
			const slot_t dj = static_cast<slot_t>(j % bits_per_slot);
			const slot_t maski = 1 << di;
			const slot_t maskj = 1 << dj;
			if (sloti == slotj)
			{
				const slot_t delta = dj - di;
				const slot_t mask = maski | maskj;
				const slot_t maskn = ~mask;
				for (size_t k = 0; k < size; k++)
				{
					slot_t x = matrix[k][sloti];
					slot_t mi = (x & maski) << delta;
					slot_t mj = (x & maskj) >> delta;
					matrix[k][sloti] = (x &maskn) | mi | mj;
				}
			}
			else
			{
				int delta = dj - di;
				const slot_t maskni = ~maski;
				const slot_t masknj = ~maskj;
				if (delta > 0)
					for (size_t k = 0; k < size; k++)
					{
						slot_t &xi = matrix[k][sloti];
						slot_t &xj = matrix[k][slotj];
						slot_t mi = (xi & maski) << delta;
						slot_t mj = (xj & maskj) >> delta;
						xi = (xi &maskni) | mj;
						xj = (xj &masknj) | mi;
					}
				else
				{
					delta = -delta;
					for (size_t k = 0; k < size; k++)
					{
						slot_t &xi = matrix[k][sloti];
						slot_t &xj = matrix[k][slotj];
						slot_t mi = (xi & maski) >> delta;
						slot_t mj = (xj & maskj) << delta;
						xi = (xi &maskni) | mj;
						xj = (xj &masknj) | mi;
					}
				}
			}
		}
		void trace_matrix(const std::vector<std::vector<slot_t>> &m, size_t hsize)
		{
			int lcount = 0;
			if (hsize > 80)
				return;
			for (auto &v : m)
			{
				size_t scount = 0;
				for (auto s : v)
					for (size_t i = 0; i < bits_per_slot && scount < hsize; i++, scount++)
						std::cout << (s & (1 << i) ? '1' : '0') << (i % 8 == 7 ? " " : "");
				std::cout << "\n" << (++lcount % 8 == 0 ? "\n" : "");
			}
			std::cout << std::flush;
		}
	};
}
#endif
