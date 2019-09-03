#include "znbasic.h"
#include <iostream>
#include <iomanip>
#include <string.h>



template <class stream_t>
class log_stream_t : public log_base_t
{
public:
	log_stream_t(stream_t &stream) : stream_(stream) {}
	log_stream_t &operator << (int value) override
	{
		stream_ << value; return *this;
	}
	log_stream_t &operator << (long long value) override
	{
		stream_ << value; return *this;
	}
	log_stream_t &operator << (size_t value) override
	{
		stream_ << value; return *this;
	}
	log_base_t &operator << (double value) override
	{
		stream_ << std::setprecision(4) << value; return *this;
	}
	log_stream_t &operator << (const std::string &text) override
	{
		stream_ << text; return *this;
	}
	log_base_t &operator << (const flush_t &) override 
	{
		stream_ << std::flush;
		return *this;
	}
	log_base_t &operator << (const newline_t &) override
	{
		stream_ << std::endl;
		return *this;
	}
private:
	stream_t &stream_;
};

class log_null_t : public log_base_t
{
public:
	log_null_t &operator << (int) override
	{
		return *this;
	}
	log_null_t &operator << (size_t) override
	{
		return *this;
	}
	log_base_t &operator << (long long) override
	{
		return *this;
	}
	log_base_t &operator << (double) override
	{
		return *this;
	}
	log_null_t &operator << (const std::string &) override
	{
		return *this;
	}
	log_base_t &operator << (const flush_t &) override
	{
		return *this;
	}
	log_base_t &operator << (const newline_t &) override
	{
		return *this;
	}
};

void log_base_t::init(int argc, char *argv[])
{
	for (int i = 0; i < argc; i++)
		if (strncmp(argv[i], "--log-level=", 12) == 0)
			if (strcmp(argv[i], "debug") == 0)
				level() = debug_e;
			else if (strcmp(argv[i], "info") == 0)
				level() = info_e;
			else if (strcmp(argv[i], "warning") == 0)
				level() = warning_e;
			else if (strcmp(argv[i], "trace") == 0)
				level() = trace_e;
}


log_base_t &log_base_t::instance(level_e l)
{
	static log_stream_t<std::ostream> stdlog(std::cout);
	static log_null_t nullog;
	if (static_cast<int>(l) <= static_cast<int>(level()))
		return stdlog;
	else
		return nullog;
}
