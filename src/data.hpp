/**
 * @file data.hpp
 * @brief Defines utilities for managing performance data logging and timing.
 *
 * This file provides constants for data directories and filenames, functions
 * for saving performance metrics (search, construction, removal times) to CSV files,
 * and timer classes (WallTimer, CPUTimer) for measuring execution time.
 * It helps in systematically collecting and organizing benchmark results.
 */
#pragma once

#include <chrono> // For WallTimer
#include <ctime>  // For CPUTimer
#include <string>
#include <fstream> // Included for get_hostname implementation detail, though not strictly needed for declarations
#include <iostream> // Included for WallTimer/CPUTimer print methods

/** @brief Directory path for storing data related to genetic datasets. */
static const std::string GENETIC_DATA_DIRECTORY = "data-genetic/";
/** @brief Directory path for storing data related to synthetic datasets. */
static const std::string SYNTHETIC_DATA_DIRECTORY = "data-synthetic/";
/** @brief Base filename prefix for search performance data files. */
static const std::string SEARCH_DATA_FILENAME = "search-data-";
/** @brief Base filename prefix for construction performance data files. */
static const std::string CONSTRUCTION_DATA_FILENAME = "construction-data-";
/** @brief Base filename prefix for removal performance data files. */
static const std::string REMOVAL_DATA_FILENAME = "removal-data-";
/** @brief Standard file extension for Comma Separated Value files. */
static const std::string CSV_EXTENSION = ".csv";

/**
 * @brief Gets the appropriate data directory path based on the data type.
 * @param is_genetic If true, returns the genetic data directory; otherwise, returns the synthetic data directory.
 * @return const std::string& A constant reference to the selected directory path string.
 */
inline const std::string& get_data_directory(bool is_genetic = false)
{
	return is_genetic ? GENETIC_DATA_DIRECTORY : SYNTHETIC_DATA_DIRECTORY;
}

/**
 * @brief Saves search performance data to a CSV file.
 *
 * Appends a record to a CSV file named based on the hostname, method, and data type.
 * The file stores parameters and timing results for search operations.
 *
 * @param method A string identifying the search method or data structure used (e.g., "SkipTrie", "ParallelZipTrie").
 * @param n The number of elements in the data structure when the search was performed.
 * @param m The length of the search key (e.g., number of characters or bits).
 * @param l The Longest Common Prefix (LCP) length between the search key and its path/result in the structure.
 * @param num_nanoseconds The total time taken for the search operation(s), in nanoseconds.
 * @param num_repetitions The number of times the search operation was repeated to get the total time. Defaults to 1.
 * @param min_par_compare A parameter potentially related to parallel comparison thresholds (specific meaning depends on context). Defaults to 0.
 * @param is_genetic Flag indicating whether the data pertains to genetic (true) or synthetic (false) datasets. Defaults to false.
 */
void save_search_data(const std::string& method, size_t n, size_t m, size_t l, size_t num_nanoseconds, size_t num_repetitions = 1, size_t min_par_compare = 0, bool is_genetic = false);

/**
 * @brief Saves construction performance data to a CSV file.
 *
 * Appends a record to a CSV file named based on the hostname, method, and data type.
 * The file stores parameters and timing results for data structure construction operations.
 *
 * @param method A string identifying the construction method or data structure used (e.g., "SkipTrie", "ParallelZipTrie").
 * @param n The number of elements inserted during the construction.
 * @param m The characteristic length of the elements being inserted (e.g., average or fixed string length).
 * @param l A parameter potentially representing the average LCP length of the inserted data, or another characteristic.
 * @param num_nanoseconds The total time taken for the construction operation(s), in nanoseconds.
 * @param num_repetitions The number of times the construction operation was repeated or batched. Defaults to 1.
 * @param min_par_compare A parameter potentially related to parallel comparison thresholds used during construction. Defaults to 0.
 * @param is_genetic Flag indicating whether the data pertains to genetic (true) or synthetic (false) datasets. Defaults to false.
 */
void save_construction_data(const std::string& method, size_t n, size_t m, size_t l, size_t num_nanoseconds, size_t num_repetitions = 1, size_t min_par_compare = 0, bool is_genetic = false);

/**
 * @brief Saves removal performance data to a CSV file.
 *
 * Appends a record to a CSV file named based on the hostname, method, and data type.
 * The file stores parameters and timing results for data structure removal operations.
 *
 * @param method A string identifying the removal method or data structure used (e.g., "SkipTrie", "ParallelZipTrie").
 * @param n The number of elements present before removal or the number of elements removed.
 * @param m The characteristic length of the elements being removed (e.g., average or fixed string length).
 * @param l A parameter potentially representing the LCP or other characteristic related to the removed elements.
 * @param num_nanoseconds The total time taken for the removal operation(s), in nanoseconds.
 * @param num_repetitions The number of times the removal operation was repeated or batched. Defaults to 1.
 * @param min_par_compare A parameter potentially related to parallel comparison thresholds used during removal. Defaults to 0.
 * @param is_genetic Flag indicating whether the data pertains to genetic (true) or synthetic (false) datasets. Defaults to false.
 */
void save_removal_data(const std::string& method, size_t n, size_t m, size_t l, size_t num_nanoseconds, size_t num_repetitions = 1, size_t min_par_compare = 0, bool is_genetic = false);

/**
 * @brief Retrieves the hostname of the machine executing the code.
 * @details Implemented by reading from `/proc/sys/kernel/hostname` on Linux systems.
 * @return std::string The hostname.
 */
std::string get_hostname();
/** @brief Constant string storing the hostname, retrieved once at program startup. */
static const std::string HOSTNAME = get_hostname();

/**
 * @brief Converts a duration in nanoseconds to a human-readable string format.
 * @details Formats the time into hours (h), minutes (m), seconds (s), and milliseconds (ms) as appropriate.
 * @param nanoseconds The duration in nanoseconds.
 * @return std::string A formatted string representing the duration (e.g., "1m 30s 500ms", "120ms").
 */
std::string pretty_print(size_t nanoseconds);

/**
 * @struct WallTimer
 * @brief A timer class measuring wall-clock time using `std::chrono::high_resolution_clock`.
 * @details Provides methods to start the timer, get the elapsed time, and print the elapsed time.
 * Suitable for measuring the real time elapsed, including potential waits or system overhead.
 */
struct WallTimer
{
	std::chrono::high_resolution_clock::time_point start_time; /**< Stores the time point when the timer was last started/reset. */

	/** @brief Default constructor. */
	WallTimer() = default; // Use default constructor

	/**
	 * @brief Calculates the elapsed time since the timer was started.
	 * @return size_t The elapsed time in nanoseconds.
	 */
	size_t elapsed_nanoseconds() const noexcept
	{
		auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start_time);
		return duration.count();
	}

	/**
	 * @brief Starts or resets the timer and optionally prints a starting message.
	 * @param prompt A message to print to stdout indicating what is being timed. If empty, no message is printed.
	 * @param reset If true (default), resets the start time. If false, keeps the previous start time (useful for cumulative timing).
	 */
	void start(const std::string& prompt = "", bool reset = true) noexcept
	{
		if (reset)
		{
			start_time = std::chrono::high_resolution_clock::now();
		}

		if (!prompt.empty())
		{
			// Use printf for potentially better performance in tight loops compared to std::cout
			printf("%s... \t", prompt.c_str());
			fflush(stdout); // Ensure the prompt is displayed immediately
		}
	}

	/**
	 * @brief Prints the elapsed time since the last start/reset to stdout.
	 * @details Calculates the elapsed time, formats it using `pretty_print`, and prints "done (formatted_time)\n".
	 * @return size_t The elapsed time in nanoseconds that was printed.
	 */
	size_t print() noexcept
	{
		auto nanoseconds = elapsed_nanoseconds();
		printf("done (%s)\n", pretty_print(nanoseconds).c_str());
		fflush(stdout); // Ensure the output is displayed immediately
		return nanoseconds;
	}
};

/**
 * @struct CPUTimer
 * @brief A timer class measuring CPU time used by the current process using `std::clock()`.
 * @details Provides methods to start the timer, get the elapsed CPU time, and print the elapsed CPU time.
 * Suitable for measuring the amount of processor time consumed by the code, excluding time spent waiting.
 * @note The resolution and accuracy of `std::clock()` can vary between systems.
 */
struct CPUTimer
{
	std::clock_t start_time; /**< Stores the clock tick count when the timer was last started/reset. */

	/** @brief Default constructor. */
	CPUTimer() = default; // Use default constructor

	/**
	 * @brief Calculates the elapsed CPU time since the timer was started.
	 * @return size_t The elapsed CPU time in nanoseconds.
	 */
	size_t elapsed_nanoseconds() const noexcept
	{
		auto end = std::clock();
		// Calculate duration in nanoseconds based on clock ticks and CLOCKS_PER_SEC
		return static_cast<size_t>(static_cast<double>(end - start_time) * 1e9 / CLOCKS_PER_SEC);
	}

	/**
	 * @brief Starts or resets the CPU timer and optionally prints a starting message.
	 * @param prompt A message to print to stdout indicating what is being timed. If empty, no message is printed.
	 * @param reset If true (default), resets the start time. If false, keeps the previous start time.
	 */
	void start(const std::string& prompt = "", bool reset = true) noexcept
	{
		if (reset)
		{
			start_time = std::clock();
		}

		if (!prompt.empty())
		{
			printf("%s... \t", prompt.c_str());
			fflush(stdout); // Ensure the prompt is displayed immediately
		}
	}

	/**
	 * @brief Prints the elapsed CPU time since the last start/reset to stdout.
	 * @details Calculates the elapsed CPU time, formats it using `pretty_print`, and prints "done (formatted_time)\n".
	 * @return size_t The elapsed CPU time in nanoseconds that was printed.
	 */
	size_t print() noexcept
	{
		auto nanoseconds = elapsed_nanoseconds();
		printf("done (%s)\n", pretty_print(nanoseconds).c_str());
		fflush(stdout); // Ensure the output is displayed immediately
		return nanoseconds;
	}
};
