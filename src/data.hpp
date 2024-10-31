#pragma once

#include <chrono>
#include <ctime>
#include <string>

static const std::string DATA_DIRECTORY = "data-v5/";
static const std::string SEARCH_DATA_FILENAME = "search-data-";
static const std::string CONSTRUCTION_DATA_FILENAME = "construction-data-";
static const std::string REMOVAL_DATA_FILENAME = "removal-data-";
static const std::string CSV_EXTENSION = ".csv";

// Save search data to a file
// n: number of elements in skip list
// m: length of search string
// l: lcp of search string
void save_search_data(const std::string& method, size_t n, size_t m, size_t l, size_t num_nanoseconds, size_t num_repetitions = 1, size_t min_par_compare = 0);

// Save construction data to a file
// n: number of elements to be inserted
// m: string length
// l: average lcp length
void save_construction_data(const std::string& method, size_t n, size_t m, size_t l, size_t num_nanoseconds, size_t num_repetitions = 1, size_t min_par_compare = 0);

// Save removal data to a file
// n: number of elements to be removed
// m: string length
void save_removal_data(const std::string& method, size_t n, size_t m, size_t num_nanoseconds, size_t num_repetitions = 1, size_t min_par_compare = 0);

std::string get_hostname();
static const std::string HOSTNAME = get_hostname();

std::string pretty_print(size_t nanoseconds);

struct WallTimer
{
	std::chrono::high_resolution_clock::time_point start_time;

	WallTimer()
	{
	}

	size_t elapsed_nanoseconds() const noexcept
	{
		auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start_time);
	
		return duration.count();
	}

	void start(const std::string& prompt = "", bool reset = true) noexcept
	{
		if (reset)
		{
			start_time = std::chrono::high_resolution_clock::now();
		}

		if (!prompt.empty())
		{
			printf("%s... \t", prompt.c_str());
		}
	}

	size_t print() noexcept
	{
		auto nanoseconds = elapsed_nanoseconds();
		printf("done (%s)\n", pretty_print(nanoseconds).c_str());
		return nanoseconds;
	}
};

struct CPUTimer
{
	std::clock_t start_time;

	CPUTimer()
	{
	}

	size_t elapsed_nanoseconds() const noexcept
	{
		auto end = std::clock();
		return (end - start_time) * 1e9 / CLOCKS_PER_SEC;
	}

	void start(const std::string& prompt = "", bool reset = true) noexcept
	{
		if (reset)
		{
			start_time = std::clock();
		}

		if (!prompt.empty())
		{
			printf("%s... \t", prompt.c_str());
		}
	}

	size_t print() noexcept
	{
		auto nanoseconds = elapsed_nanoseconds();
		printf("done (%s)\n", pretty_print(nanoseconds).c_str());
		fflush(stdout);
		return nanoseconds;
	}
};
