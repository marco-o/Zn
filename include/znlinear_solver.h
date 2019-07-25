#ifndef znlinear_solver_H
#define znlinear_solver_H
namespace zn
{
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
#if 0 // DBG_SIEVE >= DBG_SIEVE_DEBUG
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
#if 0 // DBG_SIEVE >= DBG_SIEVE_DEBUG
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
#if 0 // DBG_SIEVE >= DBG_SIEVE_DEBUG
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
#if 0 // DBG_SIEVE >= DBG_SIEVE_DEBUG
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
