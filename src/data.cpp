#include "data.hpp"

#include <fstream>
#include <string>

void save_search_data(const std::string& method, size_t n, size_t m, size_t l, size_t num_nanoseconds, size_t num_repetitions, size_t min_par_compare, bool is_genetic)
{
	std::string PREFIX = get_data_directory(is_genetic) + HOSTNAME + "-" + SEARCH_DATA_FILENAME;

	std::string filename = PREFIX + method + CSV_EXTENSION;

	std::ofstream file(filename, std::ios::app);
	file << n << "," << m << "," << l << "," << num_nanoseconds << "," << num_repetitions << "," << min_par_compare << std::endl;
}

void save_construction_data(const std::string& method, size_t n, size_t m, size_t l, size_t num_nanoseconds, size_t num_repetitions, size_t min_par_compare, bool is_genetic)
{
	std::string PREFIX = get_data_directory(is_genetic) + HOSTNAME + "-" + CONSTRUCTION_DATA_FILENAME;

	std::string filename = PREFIX + method + CSV_EXTENSION;

	std::ofstream file(filename, std::ios::app);
	file << n << "," << m << "," << l << "," << num_nanoseconds << "," << num_repetitions << "," << min_par_compare << std::endl;
}

void save_removal_data(const std::string& method, size_t n, size_t m, size_t l, size_t num_nanoseconds, size_t num_repetitions, size_t min_par_compare, bool is_genetic)
{
	std::string PREFIX = get_data_directory(is_genetic) + HOSTNAME + "-" + REMOVAL_DATA_FILENAME;

	std::string filename = PREFIX + method + CSV_EXTENSION;

	std::ofstream file(filename, std::ios::app);
	file << n << "," << m << "," << l << "," << num_nanoseconds << "," << num_repetitions << "," << min_par_compare << std::endl;
}

std::string get_hostname()
{
	std::string hostname;
	std::ifstream file("/proc/sys/kernel/hostname");
	std::getline(file, hostname);
	return hostname;
}

std::string pretty_print(size_t nanoseconds)
{
	size_t milliseconds = nanoseconds / 1e6;
	size_t seconds = nanoseconds / 1e9;
	size_t minutes = seconds / 60;
	size_t hours = minutes / 60;

	if (hours > 0)
	{
		return std::to_string(hours) + "h " + std::to_string(minutes % 60) + "m " + std::to_string(seconds % 60) + "s " + std::to_string(milliseconds % 1000) + "ms";
	}
	else if (minutes > 0)
	{
		return std::to_string(minutes) + "m " + std::to_string(seconds % 60) + "s " + std::to_string(milliseconds % 1000) + "ms";
	}
	else if (seconds > 0)
	{
		return std::to_string(seconds) + "s " + std::to_string(milliseconds % 1000) + "ms";
	}
	else
	{
		return std::to_string(milliseconds) + "ms";
	}
}