#include "src/BitString.hpp"
#include "src/SkipTrie.hpp"
#include "src/synthetic.hpp"

#include <chrono>
#include <algorithm>
#include <random>
#include <string>
#include <vector>

#include <iostream>

const std::vector<std::string> WORDS = {
	"ANT", "APPLE",
	"APPLY", "APTLY", "AMP",
	"AQUA", "LIME", "PAT", "PEA", "PEAK",
	"PEACH"
};

BitString<char> from_string(const std::string& str)
{
	return BitString<char>(str);
}

int main()
{
	auto word = get_random_word(100);
	std::cout << word << std::endl;
	// std::vector<BitString<char>> bit_strings;
	// std::transform(WORDS.begin(), WORDS.end(), std::back_inserter(bit_strings), from_string);

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
