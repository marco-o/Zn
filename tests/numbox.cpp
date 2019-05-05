//---------------------------------------------------------------------------------
//
//  Zn 
//  Copyright Marco Oman 2019
//
// Distributed under the Boost Software License, Version 1.0. 
// (See accompanying file LICENSE_1_0.txt or copy at 
// http://www.boost.org/LICENSE_1_0.txt)
//
#include <iostream>
#include <string.h>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/miller_rabin.hpp>
#include "zngroup.h"
#include "zneratosthenes_sieve.h"
#include "znquadratic_residue.h"

#include <map>
#include <string>

template <class K, class T>
class register_entry_t
{
public:
	register_entry_t(const K &k, T t)
	{
		entries()[k] = t;
	}
	static bool get(const K &k, T &t)
	{
		auto &m = entries();
		auto it = m.find(k);
		if (it == m.end())
			return false;
		t = it->second;
		return true;
	}
private:
	static std::map<K, T> &entries(void)
	{
		static std::map<K, T> the_instance;
		return the_instance;
	}
};


typedef std::map<std::string, std::string> param_map_t;

struct entry_t
{
	//std::string help;
	int(*funct)(const param_map_t &);
};
typedef int(*funct_t)(const param_map_t &);

typedef register_entry_t<std::string, funct_t> register_cmd_t;

typedef boost::multiprecision::cpp_int large_int_t;


template <class Int>
Int find_next_prime(Int n)
{
	if (n == 2)
		return n;
	if (!bit_test(n, 0))
		++n;
	for (; ; n += 2)
		if (miller_rabin_test(n, 25))
			return n;
	return 0;
}

template <class Int>
Int find_next_prime(const param_map_t &params)
{
	auto it = params.find("n");
	if (it == params.end())
		throw std::runtime_error("Missing parameter n");
	Int n(it->second);
	return find_next_prime(n);
}

template <class Int>
int print_next_prime(const param_map_t &params)
{
	auto result = find_next_prime<Int>(params);
	std::cout << result << std::endl;
	return 0;
}

template <class Int>
Int random_number(int dgt)
{
	std::string text;
	for (int i = 0; i < dgt; i++)
		text.push_back('0' + std::rand() % 10);
	return Int(text);
}


template <class Int>
int find_large_composite(const param_map_t &params)
{
	auto it = params.find("dgt");
	if (it == params.end())
		throw std::runtime_error("Missing parameter dgt (digits)");
	auto seed = params.find("seed");
	if (seed != params.end())
		std::srand(atoi(seed->second.c_str()));
	int digits = atoi(it->second.c_str());
	Int p1 = find_next_prime<Int>(random_number<Int>(digits));
	Int p2 = find_next_prime<Int>(random_number<Int>(digits));
	std::cout << (p1 * p2) << " = " << p1 << " * " << p2 << std::endl;

	return 0;
}


register_cmd_t next_prime(std::string("next-prime"), print_next_prime<large_int_t>);
register_cmd_t large_composite(std::string("large-composite"), find_large_composite<large_int_t>);

int main(int argc, char *argv[])
{
	funct_t funct;
	param_map_t params;
	for (int i = 0; i < argc; i++)
	{
		if (register_cmd_t::get(argv[i], funct))
			try
			{
				funct(params);
			}
			catch (std::exception &exc)
			{
				std::cerr << exc.what() << std::endl;
			}
		else if (strncmp(argv[i], "--", 2) == 0)
		{
			std::string param(argv[i] + 2);
			auto pos = param.find_first_of("=");
			if (pos != std::string::npos)
				params[param.substr(0, pos)] = param.substr(pos + 1);
		}
	}
	return 0;
}