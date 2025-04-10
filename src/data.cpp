/**
 * @file data.cpp
 * @brief Implements utilities for saving performance data and formatting time.
 */
#include "data.hpp"

#include <fstream>
#include <string>
#include <iostream> // Required for std::endl used indirectly by save_* functions

void save_search_data(const std::string& method, size_t n, size_t m, size_t l, size_t num_nanoseconds, size_t num_repetitions, size_t min_par_compare, bool is_genetic)
{
	// Construct the full path and filename for the CSV data file.
	std::string PREFIX = get_data_directory(is_genetic) + HOSTNAME + "-" + SEARCH_DATA_FILENAME;
	std::string filename = PREFIX + method + CSV_EXTENSION;

	// Open the file in append mode. Creates the file if it doesn't exist.
	std::ofstream file(filename, std::ios::app);
	// Write the data as a single line, comma-separated.
	file << n << "," << m << "," << l << "," << num_nanoseconds << "," << num_repetitions << "," << min_par_compare << std::endl;
	// File is automatically closed when 'file' goes out of scope.
}

void save_construction_data(const std::string& method, size_t n, size_t m, size_t l, size_t num_nanoseconds, size_t num_repetitions, size_t min_par_compare, bool is_genetic)
{
	// Construct the full path and filename for the CSV data file.
	std::string PREFIX = get_data_directory(is_genetic) + HOSTNAME + "-" + CONSTRUCTION_DATA_FILENAME;
	std::string filename = PREFIX + method + CSV_EXTENSION;

	// Open the file in append mode. Creates the file if it doesn't exist.
	std::ofstream file(filename, std::ios::app);
	// Write the data as a single line, comma-separated.
	file << n << "," << m << "," << l << "," << num_nanoseconds << "," << num_repetitions << "," << min_par_compare << std::endl;
	// File is automatically closed when 'file' goes out of scope.
}

void save_removal_data(const std::string& method, size_t n, size_t m, size_t l, size_t num_nanoseconds, size_t num_repetitions, size_t min_par_compare, bool is_genetic)
{
	// Construct the full path and filename for the CSV data file.
	std::string PREFIX = get_data_directory(is_genetic) + HOSTNAME + "-" + REMOVAL_DATA_FILENAME;
	std::string filename = PREFIX + method + CSV_EXTENSION;

	// Open the file in append mode. Creates the file if it doesn't exist.
	std::ofstream file(filename, std::ios::app);
	// Write the data as a single line, comma-separated.
	file << n << "," << m << "," << l << "," << num_nanoseconds << "," << num_repetitions << "," << min_par_compare << std::endl;
	// File is automatically closed when 'file' goes out of scope.
}

std::string get_hostname()
{
	std::string hostname;
	// Attempt to open the hostname file provided by the kernel (Linux-specific).
	std::ifstream file("/proc/sys/kernel/hostname");
	if (file.is_open()) // Check if the file was successfully opened
	{
		// Read the entire line (hostname) from the file.
		std::getline(file, hostname);
		// File is automatically closed when 'file' goes out of scope.
	}
	// Return the retrieved hostname (or an empty string if reading failed).
	return hostname;
}

std::string pretty_print(size_t nanoseconds)
{
	size_t milliseconds = nanoseconds / 1e6;
	size_t seconds = nanoseconds / 1e9;
	size_t minutes = seconds / 60;
	size_t hours = minutes / 60;

	// Build the output string based on the largest significant unit.
	if (hours > 0)
	{
		// Format as h, m, s, ms
		return std::to_string(hours) + "h " + std::to_string(minutes % 60) + "m " + std::to_string(seconds % 60) + "s " + std::to_string(milliseconds % 1000) + "ms";
	}
	else if (minutes > 0)
	{
		// Format as m, s, ms
		return std::to_string(minutes) + "m " + std::to_string(seconds % 60) + "s " + std::to_string(milliseconds % 1000) + "ms";
	}
	else if (seconds > 0)
	{
		// Format as s, ms
		return std::to_string(seconds) + "s " + std::to_string(milliseconds % 1000) + "ms";
	}
	else // Only milliseconds or less
	{
		// Format as ms
		return std::to_string(milliseconds) + "ms";
	}
}
