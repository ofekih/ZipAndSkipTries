#include "src/BitString.hpp"
#include "src/SkipTrie.hpp"
#include "src/synthetic.hpp"
#include "src/data.hpp"

#include "ctriepp/ctriepp/CTriePP.hpp"
#include "ctriepp/ctriepp/LongString.hpp"

#include <chrono>
#include <algorithm>
#include <random>
#include <string>
#include <vector>

#include <iostream>

using namespace ctriepp;

BitString<char> from_string(const std::string& str)
{
	return BitString<char>(str);
}

int main(int argc, char* argv[])
{
	if (argc != 5) {
		std::cout << "Usage: " << argv[0] << " <num_words> <word_length> <num_trials> <mean_lcp_length>\n";
		return 1;
	}

	size_t num_words = std::stoi(argv[1]);
	size_t word_length = std::stoi(argv[2]);
	unsigned num_trials = std::stoi(argv[3]);
	double mean_lcp_length = std::stod(argv[4]);

	CPUTimer timer;

	timer.start("Generating random words");

	const auto& words = get_random_words(word_length, num_words, mean_lcp_length);

	timer.print();

	timer.start("Generating BitStrings from random words");

	std::vector<BitString<char>> bit_strings;
	std::transform(words.begin(), words.end(), std::back_inserter(bit_strings), from_string);

	timer.print();

	timer.start("Generating LongStrings from random words");

	std::vector<LongString> long_strings;
	std::transform(words.begin(), words.end(), std::back_inserter(long_strings), [](const std::string& str) { return LongString(&str); });

	timer.print();

	timer.start("Creating the skip-trie " + std::to_string(num_trials) + " times");

	for (unsigned i = 0; i < num_trials; ++i)
	{
		SkipTrie<char> trie;

		for (const auto& word : bit_strings)
		{
			trie.insert(&word);
		}
	}

	timer.print();

	timer.start("Creating the ctrie++ " + std::to_string(num_trials) + " times");

	for (unsigned i = 0; i < num_trials; ++i)
	{
		CTriePP<char, false> ctriepp;

		for (const auto& word : long_strings)
		{
			ctriepp.insert(word, true);
		}
	}

	timer.print();

	return 0;
}
