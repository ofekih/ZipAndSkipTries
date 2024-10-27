#include "synthetic.hpp"

#include <algorithm>
#include <random>
#include <string>
#include <vector>

static const std::string ALPHABET = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!\"#$%&'()*+,-./:;<=>?@[\\]^_`{|}~";

std::string get_random_word(size_t length) noexcept
{
	static std::random_device rd;
	static std::mt19937 gen(0);
	// static std::mt19937 gen(rd());
	
	std::uniform_int_distribution<size_t> dis(0, ALPHABET.size() - 1);

	std::string str(length, '\0');
	std::generate_n(str.begin(), length, [&] { return ALPHABET[dis(gen)]; });

	return str;
}

std::vector<std::string> get_random_words(size_t length, size_t num_words, double mean_lcp_length) noexcept
{
	// poisson distribution
	static std::random_device rd;
	static std::mt19937 gen(0);
	// static std::mt19937 gen(rd());
	
	std::poisson_distribution<size_t> dis(mean_lcp_length);
	
	std::vector<std::string> words(num_words, get_random_word(length));

	// at the end of each word, I want to put a unique signature, to avoid duplicates
	// the length of this signature must be such that alphabet_length ^ signature_length > num_words
	// so that I can generate num_words unique signatures
	// this can be calculated using logs
	unsigned signature_length = std::ceil(std::log(num_words) / std::log(ALPHABET.size()));

	for (size_t i = 0; i < num_words; ++i)
	{
		auto& word = words[i];

		size_t i_copy = i;

		for (size_t j = 0; j < signature_length; ++j)
		{
			word[j + length - signature_length] = ALPHABET[i_copy % ALPHABET.size()];
			i_copy /= ALPHABET.size();
		}

		size_t difference = dis(gen);
		if (difference >= length - signature_length)
		{
			continue;
		}

		std::string random_word = get_random_word(length - difference - signature_length);
		std::copy(random_word.begin(), random_word.end(), word.begin() + difference);
	}

	return words;
}
