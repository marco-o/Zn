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
#include "znqssmall.h"

#include <intrin.h>
//#define HAVE_LARGE_PRIME
//#define HAVE_MULTIPLIER_STATS
template <class stream_t>
stream_t &operator<<(stream_t &stream, const std::vector<int> &v)
{
	for (auto x : v)
		stream << x << " ";
	return stream;
}

namespace zn
{
	template <>
	class montgomery_t<boost::multiprecision::cpp_int>
	{
	public:
		montgomery_t(const boost::multiprecision::cpp_int&) {}
		int power(const boost::multiprecision::cpp_int&, 
			      const boost::multiprecision::cpp_int&)
		{
			return 1;
		}
	};


	struct test_info_t
	{
		const char *file = nullptr;
		int algo = 0;
		int seeds = 0;
		int count = 0;
		int useed = 0;
		int bound = 10000;
		config_qs_t config_qs;
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
		typedef boost::multiprecision::cpp_int cpp_int_t;
		static bool fermat_test(const large_int &n)
		{
			large_int nm1 = n - 1;
#ifdef HAVE_CHECK
			large_int d = 8447990824LL;
			d = mul_mod(d, d, n);
			cpp_int_t n1 = n;
			cpp_int_t d1 = 8447990824LL;
			d1 = mul_mod(d1, d1, n1);
			cpp_int_t q1 = q;
			cpp_int_t x1 = 1;
#endif
			//
			// Begin with a single Fermat test - it excludes a lot of candidates:
			//
			large_int q(228); // We know n is greater than this, as we've excluded small factors
			large_int x = 1;
			for (; nm1; nm1 >>= 1)
			{
				if (nm1 & 1)
				{
					x = mul_mod(q, x, n);
#ifdef HAVE_CHECK
					x1 = mul_mod(q1, x1, n1);
#endif
				}
				q = mul_mod(q, q, n);
#ifdef HAVE_CHECK
				q1 = mul_mod(q1, q1, n1);
				std::cout << "q = " << q << ", q1 = " << q1 << ", x = " << x << ", x1 = " << x << std::endl;
#endif
			}
			return x == 1;
		}
		static bool fermat_test_montgomery(const large_int& n)
		{
			montgomery_t<large_int> mg(n);
			return mg.power(228, n - 1) == 1;
		}
		static bool fermat_test_montg(const large_int &n)
		{
			large_int nm1 = n - 1;
			//
			// Begin with a single Fermat test - it excludes a lot of candidates:
			//
			large_int q1(228); // We know n is greater than this, as we've excluded small factors
			large_int q2 = mul_mod(q1, q1, n);
			large_int b1 = static_cast<large_int>(1) << (msb(nm1) - 1);
			for (; b1; b1 >>= 1)
			{
				if (nm1 & b1)
				{
					q1 = mul_mod(q1, q2, n);
					q2 = mul_mod(q2, q2, n);
				}
				else
				{
					q2 = mul_mod(q1, q2, n);
					q1 = mul_mod(q1, q1, n);
				}
			}
			return q1 == 1;
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
		static void brute_force(const std::vector<input_item_t> &items, large_int largest, large_int /* smallest*/, int count)
		{
			std::cout << "Brute force factorization\n";
			small_int bound = safe_cast<small_int>(sqrt(largest));
			auto primes = eratosthenes_sieve<small_int>(bound);
			int examined = 0;
			int factored = 0;
			for (const auto &item : items)
				if (!item.declared_prime)
				{
					examined++;
					for (auto p : primes)
						if (item.n % p == 0)
						{
							std::cout << item.n << ": " << p << " x " << item.n / p << std::endl;
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
		template <class stat>
		static int quadratic_sieve(const std::vector<input_item_t> &items, large_int largest, const test_info_t &info)
		{
			int examined = 0;
			int factored = 0;
			std::cout << "Quadratic sieve\n";
			quadratic_sieve_cached_t<stat, long long> qs(160);
			quadratic_sieve_cached_t<stat, long long>::info_t qsinfo;
			qsinfo.config = info.config_qs;
			for (const auto &item : items)
				if (!item.declared_prime)
				{
					examined++;
					if (examined % 1000 == 0)
						std::cout << "Processing " << (examined / (items.size() / 100)) << "%\r" << std::flush;
					if (qs.factor(static_cast<long long>(item.n), qsinfo) > 1)
						factored++;
					if (info.count && examined >= info.count)
						break;
				}
			std::cout << "\nExamined = " << examined << ", factored = " << factored << "\n";
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
		static int fermat_prime_test(const std::vector<input_item_t> &items, int count)
		{
			int examined = 0;
			int errors = 0;
			for (const auto &item : items)
			{
				examined++;
				if (examined % 1000 == 0)
					std::cout << "Processing " << (examined / (items.size() / 100)) << "%\r" << std::flush;
				bool prime = fermat_test(item.n);
				if (prime != item.declared_prime)
					errors++;
				if (count && examined >= count)
					break;
			}
			std::cout << "FExamined = " << examined << ", errors = " << errors << "\n";
			return examined;
		}
		static int fermat_prime_test_montag(const std::vector<input_item_t>& items, int count)
		{
			int examined = 0;
			int errors = 0;
			for (const auto& item : items)
			{
				examined++;
				if (examined % 1000 == 0)
					std::cout << "Processing " << (examined / (items.size() / 100)) << "%\r" << std::flush;
				bool prime = fermat_test_montg(item.n);
				if (prime != item.declared_prime)
					errors++;
				if (count && examined >= count)
					break;
			}
			std::cout << "FExamined = " << examined << ", errors = " << errors << "\n";
			return examined;
		}
		static int fermat_prime_test_montgomery(const std::vector<input_item_t>& items, int count)
		{
			int examined = 0;
			int errors = 0;
			for (const auto& item : items)
			{
				examined++;
				if (examined % 1000 == 0)
					std::cout << "Processing " << (examined / (items.size() / 100)) << "%\r" << std::flush;
				bool prime = fermat_test_montgomery(item.n);
				if (prime != item.declared_prime)
					errors++;
				//std::cout << item.n << " " << (prime ? "p\n" : "c\n");
				if (count && examined >= count)
					break;
			}
			std::cout << "FExamined = " << examined << ", errors = " << errors << "\n";
			return examined;
		}
		static void process(const test_info_t &info)
		{
			std::cout << "Reading file.." << std::endl;
			auto items = load(info.file);
			large_int largest = 0;
			large_int smallest = std::numeric_limits<long long>::max();
			for (const auto &item : items)
			{
				if (item.n > largest)
					largest = item.n;
				else if (item.n < smallest)
					smallest = item.n;
			}
			std::cout << "Start processing: range = [" << smallest << ", " << largest << "]"
				      << " (" << msb(smallest) << ", " << msb(largest) << ")" << std::endl;
			std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
			int seeds[] = { 2, 5, 19, 41, 67, 79 };
			double examined = static_cast<double>(info.count ? info.count : items.size());
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
				examined = quadratic_sieve<stats_t>(items, largest, info);
				break;
			case 61:
				examined = quadratic_sieve<stats_none_t>(items, largest, info);
				break;
			case 7:
				examined = fermat_prime_test(items, info.count);
				break;
			case 71:
				examined = fermat_prime_test_montag(items, info.count);
				break;
			case 72:
				examined = fermat_prime_test_montgomery(items, info.count);
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
			info.bound = info.config_qs.m2 = atoi(argv[i] + 8);
		else if (strncmp(argv[i], "--algo=", 7) == 0)
			info.algo = atoi(argv[i] + 7);
		else if (strncmp(argv[i], "--count=", 8) == 0)
			info.count = atoi(argv[i] + 8);
		else if (strncmp(argv[i], "--base.a=", 9) == 0)
			info.config_qs.base_int.a = static_cast<float>(atof(argv[i] + 9));
		else if (strncmp(argv[i], "--base.b=", 9) == 0)
			info.config_qs.base_int.b = static_cast<float>(atof(argv[i] + 9));
		else if (strncmp(argv[i], "--m2.a=", 7) == 0)
			info.config_qs.m2_int.a = static_cast<float>(atof(argv[i] + 7));
		else if (strncmp(argv[i], "--m2.b=", 7) == 0)
			info.config_qs.m2_int.b = static_cast<float>(atof(argv[i] + 7));
		else if (strncmp(argv[i], "--offset=", 9) == 0)
			info.config_qs.sieve_offset = atof(argv[i] + 9);
		else if (strncmp(argv[i], "--excess=", 9) == 0)
			info.config_qs.smooth_excess = atoi(argv[i] + 9);
		else if (strncmp(argv[i], "--pquality=", 11) == 0)
			info.config_qs.pquality = atof(argv[i] + 11);
		else if (strncmp(argv[i], "--base-size=", 12) == 0)
			info.config_qs.base_size = atoi(argv[i] + 12);
		else if (strcmp(argv[i], "--more-polys") == 0)
			info.config_qs.more_polys = true;
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