/*
	Test program to develop small number factoriation
	(for partial-partial relations)
*/

#include "znbasic.h"
#include "zneratosthenes_sieve.h"
#include "znelliptic_curve_fact.h"
#include "znquadratic_residue.h"
#include "boost/multiprecision/cpp_int.hpp"
#include <vector>
#include <fstream>
#include <string>
#include <numeric>
#include <map>

namespace zn
{
	// A quadratic sieve class for factoring numbers up to (anout) 19 (2^63) digits
	template <class large_int = long long, class small_int = int>
	class quadratic_sieve_cached_t
	{
	public:
		typedef unsigned char real_t;
		enum { log_unit_e = 8 };
		struct prime_t
		{
			small_int value;
			real_t	  logp;
			std::vector<small_int> inverse;
			std::vector<small_int> residue; // quadratic residues
			prime_t(small_int v = 2) : value(v), residue(v, 0), inverse(1, 0)
			{
				logp = static_cast<real_t>(std::log(v) * log_unit_e);
				std::map<small_int, small_int> squares;
				for (small_int i = 1; i < v; i++)
				{
					inverse.push_back(std::get<2>(extended_euclidean_algorithm<small_int>(v, i)));
					small_int sq = i * i % v;
					squares.insert(std::make_pair(sq, i));
				}
				for (auto it : squares)
					residue[it.first] = it.second;
			}
		};
		struct info_t
		{
			std::vector<const prime_t *> base;
			std::vector<real_t> values;
			small_int m2;
			size_t base_size;
			small_int a;
			small_int a2;
			small_int b;
			large_int c;
			info_t(size_t m = 10000, size_t bs = 40) : values(m * 2), m2(m * 2), base_size(bs) {}
		};
		quadratic_sieve_cached_t(small_int primes = 256)
		{
			small_int bound = static_cast<small_int>(primes * std::pow(std::log(primes), 1.3));
			std::vector<small_int> prime = eratosthenes_sieve<small_int>(bound);
			if (prime.size() > static_cast<size_t>(primes))
				prime.erase(prime.begin() + primes, prime.end());
			for (auto p : prime)
				primes_.push_back(prime_t(p));
		}
		small_int factor(large_int n, info_t &info)
		{
			info.base.clear() ;
			auto it = primes_.begin();
			auto end = primes_.end();
			info.a = 1;
			small_int sqr2n = static_cast<small_int>(std::sqrt(n * 2));
			small_int m = 40000;
			small_int a0 = static_cast<small_int>(sqrt(sqr2n / m));
			for ( ; it != end ; ++it)
				if (it->residue[n % it->value]) // keep it!
					if (info.a == 1 && it->value > a0)
					{
						info.a = it->value;
						info.a2 = info.a * info.a;
						info.b = quadratic_residue<small_int>(static_cast<small_int>(n % info.a2), info.a2, info.a); // then take 'square root'...
						info.c = (info.b * info.b - n) / info.a2;
					}
					else if (info.base.size() < info.base_size)
						info.base.push_back(&(*it));
					else
						break;
			info.m2 = m * 2;
			info.values.resize(info.m2);
			std::fill(info.values.begin(), info.values.end(), 0);
			for (auto p : info.base)
			{
				small_int a = info.a % p->value;
				small_int a_1 = p->inverse[a];
				small_int a1 = (a_1 * a_1) % p->value;
				small_int t = n % p->value;
				small_int t0 = (a1 * (p->residue[t] - info.b)) % p->value;
				small_int t1 = (a1 * (-p->residue[t] - info.b)) % p->value;
				small_int index0 = (m + t0 + p->value) % p->value;
				small_int index1 = (m + t1 + p->value) % p->value;
				small_int index = std::min(index0, index1);
				small_int delta = std::max(index0, index1) - index ;
				small_int m2 = info.m2 - delta;
				small_int i = index;
				for (; i < m2 ; i += p->value)
				{
					info.values[i] += p->logp;
					info.values[i + delta] += p->logp;
				}
				if (i < info.m2)
					info.values[i] += p->logp;
			}
			real_t logp = static_cast<real_t>((std::log(sqr2n / 8) + std::log(m)) * log_unit_e);
			int count = 0;
			for (small_int i = 0 ; i < info.m2 ; i++)
				if (info.values[i] > logp)
				{
					small_int x = i - m;
					large_int f1 = info.a2 * x + 2 * info.b; // evaluate polynomial
					large_int f = f1 * x + info.c;
					if (f < 0)
						f = -1;
					for (auto p : info.base)
						while (f % p->value == 0)
						{
							f /= p->value;
						}
					if (f == 1)
						count++;
					//std::cout << i << ": value = " << static_cast<int>(info.values[i]) << " f = " << f << std::endl;
				}
//			std::cout << "found " << count << ", required = " << info.base.size() << std::endl;
			return 1;
		}
	private:
		std::vector<prime_t> primes_;
	};

	struct test_info_t
	{
		const char *file = nullptr;
		int algo = 0;
		int seeds = 0;
		int count = 0;
		int useed = 0;
		int bound = 10000;
	};
	template <class large_int, class small_int>
	class prime_tester_t
	{
	public:
		struct input_item_t
		{
			large_int	n;
			bool		declared_prime;
			std::vector<small_int>	factors;
		};

		template <class stream_t>
		static input_item_t read_item(stream_t &ist)
		{
			char ch;
			input_item_t item;

			ist >> item.n >> ch;
			item.declared_prime = (ch == 'p');
			if (!item.declared_prime)
				do
				{
					small_int tmp;
					ist >> tmp;
					if (ist)
						item.factors.push_back(tmp);
					ist >> ch;
				} while (ist);

			return item;
		}

		static bool custom_prime_test(const large_int &n, int *seeds, int count)
		{
			large_int nm1 = n - 1;
			//
			// Begin with a single Fermat test - it excludes a lot of candidates:
			//
			large_int q(228), x, y; // We know n is greater than this, as we've excluded small factors
			x = powm(q, nm1, n);
			if (x != 1u)
				return false;

			q = n - 1;
			unsigned k = lsb(q);
			q >>= k;

			//
			// Execute the trials:
			//
		//	int seeds[] = { 2, 5, 19, 41 };
			for (int i = 0; i < count; i++)
			{
				auto x = seeds[i];
				y = powm<large_int>(x, q, n);
				unsigned j = 0;
				while (true)
				{
					if (y == nm1)
						break;
					if (y == 1)
					{
						if (j == 0)
							break;
						return false; // test failed
					}
					if (++j == k)
						return false; // failed
					y = (y * y) % n; // powm<large_int>(y, 2, n);
				}
			}
			return true;  // Yeheh! probably prime.
		}

		static std::vector<input_item_t> load(const char *file)
		{
			std::vector<input_item_t> result;
			std::ifstream ifile(file);
			while (ifile)
			{
				std::string line;
				std::getline(ifile, line);
				if (line.empty())
					break;
				std::istringstream ist(line);
				result.push_back(read_item(ist));
			}
			return result;
		}

		static void full_test(const std::vector<input_item_t> &items, large_int largest, large_int smallest, int count)
		{
			small_int bound = safe_cast<small_int>(sqrt(largest));
			auto primes = eratosthenes_sieve<small_int>(bound);
			small_int begin = safe_cast<small_int>(sqrt(smallest));
			auto it = std::lower_bound(primes.begin(), primes.end(), begin - 1);
			primes.erase(primes.begin(), it);
			int errors = 0;
			int composite = 0;
			int declared = 0;
			int processing = 0;
			for (const auto &item : items)
			{
				if (++processing % 1000 == 0)
					std::cout << "Processing " << (processing / (items.size() / 100)) << "%\r" << std::flush;
				if (item.declared_prime)
				{
					declared++;
					for (auto p : primes)
						if (item.n % p == 0)
						{
							errors++;
							break;
						}
				}
				else
					composite++;
				if (count && processing > count)
					break;
			}
			std::cout << "Examined: " << items.size() << "\n"
				<< "Declared prime = " << declared << "\n"
				<< "Errors = " << errors << "\n"
				<< "composite = " << composite << "\n";
		}
		static void brute_force(const std::vector<input_item_t> &items, large_int largest, large_int smallest, int count)
		{
			std::cout << "Brute force factorization\n";
			small_int bound = safe_cast<small_int>(sqrt(largest));
			auto primes = eratosthenes_sieve<small_int>(bound);
			small_int begin = safe_cast<small_int>(sqrt(smallest));
			auto it = std::lower_bound(primes.begin(), primes.end(), begin - 1);
			primes.erase(primes.begin(), it);
			int examined = 0;
			int factored = 0;
			for (const auto &item : items)
				if (!item.declared_prime)
				{
					examined++;
					for (auto p : primes)
						if (item.n % p == 0)
						{
							factored++;
							break;
						}
					if (examined % 1000 == 0)
						std::cout << "Examined " << examined << "%\r" << std::flush;
					if (count && examined >= count)
						break;
				}
			std::cout << "Examined: " << examined << ", factored = " << factored << "\n";
		}
		static int pollard_p1(const std::vector<input_item_t> &items, int *seeds,  const test_info_t &info)
		{
			pollard_p1_t<large_int, small_int> p1;
			p1.init(info.bound);
			int examined = 0;
			int factored = 0;
			std::cout << "Pollard p-1 factorization\n";
			if (info.useed > 0)
				seeds[0] = info.useed;
			for (const auto &item : items)
				if (!item.declared_prime)
				{
					examined++;
					if (examined % 1000 == 0)
						std::cout << "Processing " << (examined / (items.size() / 100)) << "%\r" << std::flush;
					for (int i = 0; i < info.seeds; i++)
					{
						auto f = p1.fact(item.n, seeds[i]);
						if (f != 1 && f != item.n)
						{
//							small_int f1 = safe_cast<small_int>(f);
//							small_int f2 = safe_cast<small_int>(item.n / f1);
							factored++;
							break;
						}
					}
					if (info.count && examined >= info.count)
						break;
				}
			std::cout << "Examined = " << examined << ", factored = " << factored << "\n";
			return examined;
		}
		static void pollard_rho(const std::vector<input_item_t> &items, int *seeds, const test_info_t &info)
		{
			int examined = 0;
			int factored = 0;
			std::cout << "Pollard rho factorization\n";
			if (info.useed > 0)
				seeds[0] = info.useed;
			for (const auto &item : items)
				if (!item.declared_prime)
				{
					examined++;
					if (examined % 1000 == 0)
						std::cout << "Processing " << (examined / (items.size() / 100)) << "%\r" << std::flush;
					for (int i = 0; i < info.seeds; i++)
					{
						large_int f = pollards_rho(item.n, info.bound, seeds[i]);
						if (f != 1 && f != item.n)
						{
							small_int f1 = safe_cast<small_int>(f);
							small_int f2 = safe_cast<small_int>(item.n / f1);
							factored++;
							break;
						}
					}
					if (info.count && examined >= info.count)
						break;
				}
			std::cout << "Examined = " << examined << ", factored = " << factored << "\n";
		}
		static int quadratic_sieve(const std::vector<input_item_t> &items, large_int largest, const test_info_t &info)
		{
			int examined = 0;
			int factored = 0;
			std::cout << "Quadratic sieve\n";
			quadratic_sieve_cached_t<long long> qs(128);
			quadratic_sieve_cached_t<long long>::info_t qsinfo;
			for (const auto &item : items)
				if (!item.declared_prime)
				{
					examined++;
					if (examined % 1000 == 0)
						std::cout << "Processing " << (examined / (items.size() / 100)) << "%\r" << std::flush;
					qs.factor(static_cast<long long>(item.n), qsinfo);
					if (info.count && examined >= info.count)
						break;
				}
			std::cout << "Examined = " << examined << ", factored = " << factored << "\n";
			return examined;
		}
		static void elliptic_curve_1(const std::vector<input_item_t> &items, int *seeds, const test_info_t &info)
		{
			int examined = 0;
			int factored = 0;
			std::cout << "Elliptic curve factorization\n";
			if (info.useed > 0)
				seeds[0] = info.useed;
			elliptic_curve_t<large_int> ec(info.bound);
			for (const auto &item : items)
				if (!item.declared_prime)
				{
					examined++;
					if (examined % 1000 == 0)
						std::cout << "Processing " << (examined / (items.size() / 100)) << "%\r" << std::flush;
					elliptic_curve_t<large_int>::point_t pt{ seeds[1], seeds[2] };
					if (ec.run(item.n, seeds[0], pt))
						factored++;

					if (info.count && examined >= info.count)
						break;
				}
			std::cout << "Examined = " << examined << ", factored = " << factored << "\n";
		}
		static int prime_test(const std::vector<input_item_t> &items, int *seeds, int seed_count, int count)
		{
			int examined = 0;
			int errors = 0;
			std::cout << "prime test with " << seed_count << " seeds" << std::endl;
			for (const auto &item : items)
			{
				examined++;
				if (examined % 1000 == 0)
					std::cout << "Processing " << (examined / (items.size() / 100)) << "%\r" << std::flush;
				bool prime = custom_prime_test(item.n, seeds, seed_count);
				if (prime != item.declared_prime)
					errors++;
				if (count && examined >= count)
					break;
			}
			std::cout << "Examined = " << examined << ", errors = " << errors << "\n";
			return examined;
		}

		static void process(const test_info_t &info)
		{
			std::cout << "Reading file.." << std::endl;
			auto items = load(info.file);
			large_int largest = 0;
			large_int smallest = std::numeric_limits<small_int>::max();
			for (const auto &item : items)
			{
				if (item.n > largest)
					largest = item.n;
				else if (item.n < smallest)
					smallest = item.n;
			}
			std::cout << "Start processing: range = [" << smallest << ", " << largest << "]" << std::endl;
			std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
			int seeds[] = { 2, 5, 19, 41, 67, 79 };
			double examined = (info.count ? info.count : items.size());
			switch (info.algo)
			{
			case 0:
				if (info.seeds >= 0 && info.seeds <= 6)
					examined = prime_test(items, seeds, info.seeds, info.count);
				break;
			case 1:
				full_test(items, largest, smallest, info.count);
				break;
			case 2:
				brute_force(items, largest, smallest, info.count);
				break;
			case 3:
				examined = pollard_p1(items, seeds, info);
				break;
			case 4:
				pollard_rho(items, seeds, info);
				break;
			case 5:
				elliptic_curve_1(items, seeds, info);
				break;
			case 6:
				examined = quadratic_sieve(items, largest, info);
				break;
			}
			auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
			std::cout << "Average iteration time = " << elapsed / examined << "us\n";
		}
	};
}

int main(int argc, char *argv[])
{
	bool use_cppint = true;
	zn::test_info_t info;
	for (int i = 0; i < argc; i++)
		if (strncmp(argv[i], "--file=", 7) == 0)
			info.file = argv[i] + 7;
		else if (strncmp(argv[i], "--seeds=", 8) == 0)
			info.seeds = atoi(argv[i] + 8);
		else if (strncmp(argv[i], "--useed=", 8) == 0)
			info.useed = atoi(argv[i] + 8);
		else if (strncmp(argv[i], "--bound=", 8) == 0)
			info.bound = atoi(argv[i] + 8);
		else if (strncmp(argv[i], "--algo=", 7) == 0)
			info.algo = atoi(argv[i] + 7);
		else if (strncmp(argv[i], "--count=", 8) == 0)
			info.count = atoi(argv[i] + 8);
		else if (strcmp(argv[i], "--use-long") == 0)
			use_cppint = false;
	if (info.file)
		if (use_cppint)
			zn::prime_tester_t<boost::multiprecision::cpp_int, long long>::process(info);
		else
			zn::prime_tester_t<long long, long>::process(info);
	else
		std::cout << "no file specified\n";
	return 0;
}