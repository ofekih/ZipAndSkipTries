#include "src/BitString.cuh"
#include "src/SkipTrie.hpp"
#include "src/ParallelSkipTrie.cuh"
#include "src/synthetic.hpp"
// #include "src/ZipZipTrie.hpp"
#include "src/ZipTrie.hpp"
#include "src/ParallelZipTrie.cuh"

#include <chrono>
#include <algorithm>
#include <random>
#include <string>
#include <vector>

#include <iostream>

const std::vector<std::string> WORDS = {
	"ANT", "APPLE", "APPLY", "APT", "APTLY", "AQUA", "PAL", "PEACE", "PEACH", "PEAK", "PENGUIN", "PALA", "AAAAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAVB"
};

BitString<char, 64> from_string(const std::string& str)
{
	return BitString<char, 64>(str);
}

// ZipZipTrie<char, false> create_test_zzt()
// {
// 	static std::unordered_map<std::string, BitString<char>> word_to_bs;
// 	static std::unordered_map<std::string, unsigned> word_to_index;
	
// 	if (word_to_bs.empty())
// 	{
// 		for (unsigned i = 0; i < WORDS.size(); ++i)
// 		{
// 			word_to_bs[WORDS[i]] = from_string(WORDS[i]);
// 			word_to_index[WORDS[i]] = i;
// 		}
// 	}

// 	unsigned longest_word_length = std::max_element(WORDS.begin(), WORDS.end(), [](const std::string& a, const std::string& b) { return a.size() < b.size(); })->size();

// 	ZipZipTrie<char, false> trie(WORDS.size(), longest_word_length);
// 	const unsigned NULLPTR = std::numeric_limits<unsigned>::max();

// 	trie.set_root_index(word_to_index["ANT"]);

// 	trie.set(&word_to_bs["ANT"], { 3, 52 }, NULLPTR, word_to_index["APPLY"], 0, 0);
// 	trie.set(&word_to_bs["APPLE"], { 2, 12 }, NULLPTR, NULLPTR, 4, 1);
// 	trie.set(&word_to_bs["APPLY"], { 3, 21 }, word_to_index["APPLE"], word_to_index["APT"], 1, 0);
// 	trie.set(&word_to_bs["APT"], { 3, 3 }, NULLPTR, word_to_index["PEAK"], 2, 0);
// 	trie.set(&word_to_bs["APTLY"], { 0, 45 }, NULLPTR, NULLPTR, 1, 3);
// 	trie.set(&word_to_bs["AQUA"], { 2, 0 }, word_to_index["APTLY"], word_to_index["PEACE"], 0, 1);
// 	trie.set(&word_to_bs["PAL"], { 0, 33 }, NULLPTR, NULLPTR, 1, 0);
// 	trie.set(&word_to_bs["PEACE"], { 1, 8 }, word_to_index["PAL"], word_to_index["PEACH"], 0, 3);
// 	trie.set(&word_to_bs["PEACH"], { 1, 2 }, NULLPTR, NULLPTR, 4, 3);
// 	trie.set(&word_to_bs["PEAK"], { 2, 9 }, word_to_index["AQUA"], word_to_index["PENGUIN"], 0, 0);
// 	trie.set(&word_to_bs["PENGUIN"], { 0, 1 }, NULLPTR, NULLPTR, 2, 0);

// 	return trie;
// }

int main()
{
	// auto word = get_random_word(100);
	// std::cout << word << std::endl;
	std::vector<BitString<char, 64>> bit_strings;
	std::transform(WORDS.begin(), WORDS.end(), std::back_inserter(bit_strings), from_string);

	unsigned longest_word_length = std::max_element(WORDS.begin(), WORDS.end(), [](const std::string& a, const std::string& b) { return a.size() < b.size(); })->size();

	ZipTrie<char, true, GeometricRank, 64> trie(WORDS.size(), longest_word_length);
	// ParallelZipTrie<char, false, GeometricRank, 64> trie(WORDS.size(), longest_word_length);
	// ParallelSkipTrie<char, 64> trie(longest_word_length);
	// ZipZipTrie<char, unsigned> trie = create_test_zzt();

	// if (!trie.contains(&bit_strings[4]))
	// {
	// 	std::cout << "Word " << bit_strings[4].to_string() << " not found in trie" << std::endl;
	// }


	trie.insert(&bit_strings[0]);
	trie.insert(&bit_strings[1]);
	trie.insert(&bit_strings[2]);
	trie.insert(&bit_strings[3]);
	trie.insert(&bit_strings[4]);
	trie.insert(&bit_strings[5]);
	trie.insert(&bit_strings[6]);
	trie.insert(&bit_strings[7]);
	trie.insert(&bit_strings[8]);
	trie.insert(&bit_strings[9]);
	trie.insert(&bit_strings[10]);
	trie.insert(&bit_strings[11]);
	// trie.insert(&bit_strings[12]);
	trie.insert(&bit_strings[13]);

	for (const auto& word : bit_strings)
	{
		if (!trie.contains(&word))
		{
			// should printf something like: <word> not found in trie, but shares an LCP length of <lcp_length>
			// use trie.lcp(&word) to get the LCP length
			printf("%s not found in trie, but shares an LCP length of %u\n", word.to_string().c_str(), trie.lcp(&word));
		}
	}

	trie.to_dot("test.dot");

	// use system to run dot -Tpng test.dot -o test.png
	int result = system("dot -Tpng test.dot -o test.png");
	if (result != 0) {
		std::cerr << "Error: Failed to generate PNG from DOT file. Return code: " << result << std::endl;
		return 1;
	}

	// SkipTrie<char> trie;
	// for (const auto& bit_string : bit_strings)
	// {
	// 	printf("Inserting %s\n", bit_string.to_string().c_str());
	// 	trie.insert(&bit_string);
	// 	std::cout << trie << std::endl;

	// 	printf("\n\n");
	// }

	return 0;
}
