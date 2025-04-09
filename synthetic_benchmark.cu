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

BitString<char> from_string(const std::string& str)
{
	return BitString<char>(str);
}

struct DataPair
{
	BitString<char> bit_string;
	LongString long_string;
};

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
