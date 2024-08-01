#include "synthetic.hpp"

#include <algorithm>
#include <random>
#include <string>
#include <vector>

static const std::string ALPHABET = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!\"#$%&'()*+,-./:;<=>?@[\\]^_`{|}~";

std::string get_random_word(size_t length) noexcept
{
	static std::random_device rd;
	static std::mt19937 gen(rd());
	
	std::uniform_int_distribution<size_t> dis(0, ALPHABET.size() - 1);

	std::string str(length, '\0');
	std::generate_n(str.begin(), length, [&] { return ALPHABET[dis(gen)]; });

	return str;
}

std::vector<std::string> get_random_words(size_t length, size_t num_words, double mean_lcp_length) noexcept
{
	// poisson distribution
	static std::random_device rd;
	static std::mt19937 gen(rd());
	
	std::poisson_distribution<size_t> dis(mean_lcp_length);
	
	std::vector<std::string> words(num_words, get_random_word(length));

	for (auto& word : words)
	{
		size_t difference = dis(gen);
		if (difference >= length)
		{
			continue;
		}

		std::string random_word = get_random_word(length - difference);
		std::copy(random_word.begin(), random_word.end(), word.begin() + difference);
	}

	return words;
}
