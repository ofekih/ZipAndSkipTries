/**
 * @file synthetic_benchmark.cu
 * @brief Benchmarking program for trie data structures using synthetic data.
 *
 * @details This program benchmarks various trie implementations (SkipTrie, ZipTrie,
 * ParallelSkipTrie, ParallelZipTrie, and C-Trie++) on synthetic data with controlled
 * properties. It measures and compares construction, search, and removal performance
 * across different implementations, with configurable parameters for word length,
 * number of words, and average LCP length. Results are saved to CSV files for further
 * analysis and visualization.
 *
 * @see genetic_benchmark.cu
 * @see src/synthetic.hpp
 * @see src/BitString.cuh
 * @see src/SkipTrie.hpp
 * @see src/ZipTrie.hpp
 * @see src/ParallelSkipTrie.cuh
 * @see src/ParallelZipTrie.cuh
 */

#include "src/BitString.cuh"
#include "src/SkipTrie.hpp"
#include "src/ParallelSkipTrie.cuh"
#include "src/ZipTrie.hpp"
#include "src/ParallelZipTrie.cuh"
#include "src/synthetic.hpp"
#include "src/data.hpp"

#include "ctriepp/ctriepp/CTriePP.hpp"
#include "ctriepp/ctriepp/LongString.hpp"

#include <chrono>
#include <algorithm>
#include <random>
#include <string>
#include <vector>
#include <locale.h>

#include <iostream>

using namespace ctriepp;

/**
 * @brief Creates a BitString from a standard string.
 *
 * @details Utility function that converts a standard string to a BitString<char>.
 * This is used for creating BitString objects from the synthetic data strings.
 *
 * @param str The standard string to convert.
 * @return BitString<char> The resulting BitString object.
 */
BitString<char> from_string(const std::string& str)
{
	return BitString<char>(str);
}

/**
 * @struct DataPair
 * @brief A structure holding both BitString and LongString representations of synthetic data.
 *
 * @details This structure is used to store the same string in two different formats:
 * as a BitString (for our trie implementations) and as a LongString (for C-Trie++),
 * allowing for fair comparison between different data structures.
 */
struct DataPair
{
	BitString<char> bit_string; ///< String represented as BitString (for SkipTrie/ZipTrie variants)
	LongString long_string; ///< String represented as LongString (for C-Trie++)
};

/**
 * @brief Generates synthetic data with controlled properties.
 *
 * @details Creates a vector of DataPair objects containing randomly generated strings
 * with specified properties. The strings are generated using the get_random_words function
 * from synthetic.hpp, which creates strings with a controlled average LCP length.
 * Each string is stored both as a BitString and as a LongString.
 *
 * @param num_words The number of words to generate.
 * @param word_length The length of each word.
 * @param mean_lcp_length The average Longest Common Prefix length between words.
 * @return std::vector<DataPair> A vector of DataPair objects containing the generated data.
 *
 * @see get_random_words
 */
inline std::vector<DataPair> generate_data(size_t num_words, size_t word_length, double mean_lcp_length)
{
	static std::vector<std::string> words;

	words = get_random_words(word_length, num_words, mean_lcp_length);

	std::vector<DataPair> data;
	data.reserve(num_words);

	std::transform(words.begin(), words.end(), std::back_inserter(data), [](const std::string& str) {
		return DataPair{ BitString<char>(str), LongString(&str) };
	});

	return data;
}

/**
 * @brief Runs construction benchmarks on various trie implementations using synthetic data.
 *
 * @details Measures the time taken to construct different trie data structures (SkipTrie,
 * ParallelSkipTrie, memory-intensive ZipTrie, memory-efficient ZipTrie, memory-intensive
 * ParallelZipTrie, memory-efficient ParallelZipTrie, and C-Trie++) with the same synthetic
 * data. For each implementation, the benchmark is repeated multiple times and the results
 * are saved to CSV files.
 *
 * @param num_words The number of words to use in the benchmark.
 * @param word_length The length of each word.
 * @param mean_lcp The average Longest Common Prefix length between words.
 * @param num_repetitions The number of repetitions for each benchmark (default: 1000).
 *
 * @see save_construction_data
 */
void run_construction_benchmark(size_t num_words, size_t word_length, double mean_lcp, unsigned num_repetitions = 1000)
{
	CPUTimer timer;

	auto data = generate_data(num_words, word_length, mean_lcp);

	timer.start();

	for (unsigned i = 0; i < num_repetitions; ++i)
	{
		SkipTrie<char> trie;

		for (const auto& word : data)
		{
			trie.insert(&word.bit_string);
		}
	}

	save_construction_data("skip-trie", num_words, word_length, mean_lcp, timer.elapsed_nanoseconds(), num_repetitions, 0, false);

	timer.start();

	for (unsigned i = 0; i < num_repetitions; ++i)
	{
		ParallelSkipTrie<char> trie(word_length);

		for (const auto& word : data)
		{
			trie.insert(&word.bit_string);
		}
	}

	save_construction_data("parallel-skip-trie", num_words, word_length, mean_lcp, timer.elapsed_nanoseconds(), num_repetitions, MIN_PAR_COMPARE_WORD_SIZE, false);

	timer.start();

	for (unsigned i = 0; i < num_repetitions; ++i)
	{
		ZipTrie<char, false> trie(num_words, word_length);

		for (const auto& word : data)
		{
			trie.insert(&word.bit_string);
		}
	}

	save_construction_data("memory-intensive-zip-trie", num_words, word_length, mean_lcp, timer.elapsed_nanoseconds(), num_repetitions, 0, false);

	timer.start();

	for (unsigned i = 0; i < num_repetitions; ++i)
	{
		ParallelZipTrie<char, false> trie(num_words, word_length);

		for (const auto& word : data)
		{
			trie.insert(&word.bit_string);
		}
	}

	save_construction_data("parallel-memory-intensive-zip-trie", num_words, word_length, mean_lcp, timer.elapsed_nanoseconds(), num_repetitions, MIN_PAR_COMPARE_WORD_SIZE, false);

	timer.start();

	for (unsigned i = 0; i < num_repetitions; ++i)
	{
		ZipTrie<char, true> trie(num_words, word_length);

		for (const auto& word : data)
		{
			trie.insert(&word.bit_string);
		}
	}

	save_construction_data("zip-trie", num_words, word_length, mean_lcp, timer.elapsed_nanoseconds(), num_repetitions, 0, false);

	timer.start();

	for (unsigned i = 0; i < num_repetitions; ++i)
	{
		ParallelZipTrie<char, true> trie(num_words, word_length);

		for (const auto& word : data)
		{
			trie.insert(&word.bit_string);
		}
	}

	save_construction_data("parallel-zip-trie", num_words, word_length, mean_lcp, timer.elapsed_nanoseconds(), num_repetitions, MIN_PAR_COMPARE_WORD_SIZE, false);

	timer.start();

	for (unsigned i = 0; i < num_repetitions; ++i)
	{
		CTriePP<char, false> ctriepp;

		for (const auto& word : data)
		{
			ctriepp.insert(word.long_string, true);
		}
	}

	save_construction_data("ctrie++", num_words, word_length, mean_lcp, timer.elapsed_nanoseconds(), num_repetitions, 0, false);
}

/**
 * @brief Runs both positive and negative search benchmarks.
 *
 * @details This is a wrapper function that calls both run_contains_true_benchmark and
 * run_contains_false_benchmark to measure search performance for both existing and
 * non-existing elements in the tries.
 *
 * @param num_words The number of words to use in the benchmark.
 * @param word_length The length of each word.
 * @param mean_lcp The average Longest Common Prefix length between words.
 * @param num_repetitions The number of repetitions for each benchmark (default: 1000).
 *
 * @see run_contains_true_benchmark
 * @see run_contains_false_benchmark
 */
void run_search_benchmark(size_t num_words, size_t word_length, size_t mean_lcp, unsigned num_repetitions = 1000)
{
	auto data = generate_data(num_words + num_repetitions, word_length, mean_lcp);

	SkipTrie<char> skip_trie;
	ParallelSkipTrie<char> parallel_skip_trie(word_length);
	ZipTrie<char, false> memory_intensive_zip_trie(num_words, word_length);
	ParallelZipTrie<char, false> parallel_memory_intensive_zip_trie(num_words, word_length);
	ZipTrie<char, true> zip_trie(num_words, word_length);
	ParallelZipTrie<char, true> parallel_zip_trie(num_words, word_length);
	CTriePP<char, false> ctriepp;

	for (size_t i = 0; i < num_words; ++i)
	{
		skip_trie.insert(&data[i].bit_string);
		parallel_skip_trie.insert(&data[i].bit_string);
		memory_intensive_zip_trie.insert(&data[i].bit_string);
		parallel_memory_intensive_zip_trie.insert(&data[i].bit_string);
		zip_trie.insert(&data[i].bit_string);
		parallel_zip_trie.insert(&data[i].bit_string);
		ctriepp.insert(data[i].long_string, true);
	}

	CPUTimer timer;

	bool all_found = true;

	timer.start();

	for (unsigned i = 0; i < num_repetitions; ++i)
	{
		all_found &= skip_trie.contains(&data[num_words + i].bit_string);
	}

	save_search_data("skip-trie", num_words, word_length, mean_lcp, timer.elapsed_nanoseconds(), num_repetitions, 0, false);

	timer.start();

	for (unsigned i = 0; i < num_repetitions; ++i)
	{
		all_found &= parallel_skip_trie.contains(&data[num_words + i].bit_string);
	}

	save_search_data("parallel-skip-trie", num_words, word_length, mean_lcp, timer.elapsed_nanoseconds(), num_repetitions, MIN_PAR_COMPARE_WORD_SIZE, false);

	timer.start();

	for (unsigned i = 0; i < num_repetitions; ++i)
	{
		all_found &= memory_intensive_zip_trie.contains(&data[num_words + i].bit_string);
	}

	save_search_data("memory-intensive-zip-trie", num_words, word_length, mean_lcp, timer.elapsed_nanoseconds(), num_repetitions, 0, false);

	timer.start();

	for (unsigned i = 0; i < num_repetitions; ++i)
	{
		all_found &= parallel_memory_intensive_zip_trie.contains(&data[num_words + i].bit_string);
	}

	save_search_data("parallel-memory-intensive-zip-trie", num_words, word_length, mean_lcp, timer.elapsed_nanoseconds(), num_repetitions, MIN_PAR_COMPARE_WORD_SIZE, false);

	timer.start();

	for (unsigned i = 0; i < num_repetitions; ++i)
	{
		all_found &= zip_trie.contains(&data[num_words + i].bit_string);
	}

	save_search_data("zip-trie", num_words, word_length, mean_lcp, timer.elapsed_nanoseconds(), num_repetitions, 0, false);

	timer.start();

	for (unsigned i = 0; i < num_repetitions; ++i)
	{
		all_found &= parallel_zip_trie.contains(&data[num_words + i].bit_string);
	}

	save_search_data("parallel-zip-trie", num_words, word_length, mean_lcp, timer.elapsed_nanoseconds(), num_repetitions, MIN_PAR_COMPARE_WORD_SIZE, false);

	timer.start();

	for (unsigned i = 0; i < num_repetitions; ++i)
	{
		all_found &= ctriepp.contains(data[num_words + i].long_string);
	}

	save_search_data("ctrie++", num_words, word_length, mean_lcp, timer.elapsed_nanoseconds(), num_repetitions, 0, false);

	if (all_found)
	{
		printf("Error: all elements were found\n");
		printf("These print statements are here to prevent the compiler from optimizing out the search\n");
	}
}

/**
 * @brief Runs benchmarks with varying LCP lengths.
 *
 * @details Runs construction and search benchmarks for different average LCP lengths,
 * starting from 4 and doubling until reaching the word length. This helps analyze how
 * the performance of different trie implementations varies with the LCP length.
 *
 * @param num_words The number of words to use in the benchmark.
 * @param word_length The length of each word.
 * @param num_repetitions The number of repetitions for each benchmark (default: 1000).
 */
void run_variable_lcp_benchmarks(size_t num_words, size_t word_length, unsigned num_repetitions = 1000)
{
	WallTimer timer;

	for (double mean_lcp = 4; mean_lcp <= word_length; mean_lcp *= 2)
	{
		timer.start("Running benchmarks for mean LCP length " + std::to_string(mean_lcp).substr(0, std::to_string(mean_lcp).find('.') + 2));

		run_construction_benchmark(num_words, word_length, mean_lcp, num_repetitions);
		run_search_benchmark(num_words, word_length, mean_lcp, num_repetitions);

		timer.print();
	}
}

/**
 * @brief Runs benchmarks with varying word lengths.
 *
 * @details Runs construction and search benchmarks for different word lengths,
 * starting from 16 and doubling until reaching the maximum word length. This helps
 * analyze how the performance of different trie implementations varies with the word length.
 *
 * @param num_words The number of words to use in the benchmark.
 * @param max_word_length The maximum word length to test.
 * @param mean_lcp The average Longest Common Prefix length between words.
 * @param num_repetitions The number of repetitions for each benchmark (default: 1000).
 */
void run_variable_word_length_benchmarks(size_t num_words, size_t max_word_length, double mean_lcp, unsigned num_repetitions = 1000)
{
	WallTimer timer;

	for (size_t word_length = 4; word_length <= max_word_length; word_length *= 2)
	{
		timer.start("Running benchmarks for word length " + std::to_string(word_length));

		run_construction_benchmark(num_words, word_length, mean_lcp, num_repetitions);
		run_search_benchmark(num_words, word_length, mean_lcp, num_repetitions);

		timer.print();
	}
}

/**
 * @brief Runs benchmarks with varying numbers of words.
 *
 * @details Runs construction and search benchmarks for different numbers of words,
 * starting from 4 and doubling until reaching the maximum number of words. This helps
 * analyze how the performance of different trie implementations scales with the number
 * of elements in the trie.
 *
 * @param max_num_words The maximum number of words to test.
 * @param word_length The length of each word.
 * @param mean_lcp The average Longest Common Prefix length between words.
 * @param num_repetitions The number of repetitions for each benchmark (default: 1000).
 */
void run_variable_num_words_benchmarks(size_t max_num_words, size_t word_length, double mean_lcp, unsigned num_repetitions = 1000)
{
	WallTimer timer;

	for (size_t num_words = 4; num_words <= max_num_words; num_words *= 2)
	{
		timer.start("Running benchmarks for num words " + std::to_string(num_words));

		run_construction_benchmark(num_words, word_length, mean_lcp, num_repetitions);
		run_search_benchmark(num_words, word_length, mean_lcp, num_repetitions * 100);

		timer.print();
	}
}

/**
 * @brief Main function that runs the synthetic benchmarks.
 *
 * @details Runs construction and search benchmarks on various trie implementations
 * using synthetic data with controlled properties. The benchmark parameters (number of words,
 * word length, mean LCP length, and number of repetitions) can be specified as command-line
 * arguments.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line argument strings.
 * @return int Returns 0 on successful execution, 1 on invalid arguments.
 */
int main(int argc, char* argv[])
{
	if (argc != 5) {
		std::cout << "Usage: " << argv[0] << " <num_words> <word_length> <mean_lcp_length> <num_repetitions>\n";
		return 1;
	}

	size_t num_words = std::stoi(argv[1]);
	size_t word_length = std::stoi(argv[2]);
	double mean_lcp = std::stod(argv[3]);
	unsigned num_repetitions = std::stoi(argv[4]);

	WallTimer timer;

	setlocale(LC_NUMERIC, "");
	printf("\tRunning benchmarks for num words %'zu, word length %'zu, mean LCP length %'.1f:\n", num_words, word_length, mean_lcp);

	timer.start("\t\tConstruction");

	run_construction_benchmark(num_words, word_length, mean_lcp, num_repetitions);

	timer.print();

	timer.start("\t\tSearch");

	run_search_benchmark(num_words, word_length, mean_lcp, num_repetitions * 10);

	timer.print();

	return 0;
}
