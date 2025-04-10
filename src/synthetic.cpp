// Implementation of synthetic data generation utilities
// See synthetic.hpp for comprehensive documentation

#include "synthetic.hpp"

#include <algorithm>
#include <random>
#include <string>
#include <vector>

// The alphabet used for generating random strings
// Contains lowercase and uppercase letters, digits, and special characters
static const std::string ALPHABET = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!\"#$%&'()*+,-./:;<=>?@[\\]^_`{|}~";

// Implementation of get_random_word
// See synthetic.hpp for detailed documentation
std::string get_random_word(size_t length) noexcept
{
	// Use a static random device for initialization of the generator
	static std::random_device rd;
	// Uncomment for deterministic generation with fixed seed:
	// static std::mt19937 gen(0);
	// Use random seed for non-deterministic generation
	static std::mt19937 gen(rd());
	
	// Create uniform distribution for selecting characters from the alphabet
	std::uniform_int_distribution<size_t> dis(0, ALPHABET.size() - 1);

	// Pre-allocate string of requested length
	std::string str(length, '\0');
	// Fill the string with randomly selected characters from the alphabet
	std::generate_n(str.begin(), length, [&] { return ALPHABET[dis(gen)]; });

	return str;
}

// Implementation of get_random_words
// See synthetic.hpp for detailed documentation
std::vector<std::string> get_random_words(size_t length, size_t num_words, double mean_lcp_length) noexcept
{
	// Use Poisson distribution to generate varying LCP lengths with specified mean
	static std::random_device rd;
	// Uncomment for deterministic generation with fixed seed:
	// static std::mt19937 gen(0);
	// Use random seed for non-deterministic generation
	static std::mt19937 gen(rd());
	
	// Create distribution for common prefix lengths
	std::poisson_distribution<size_t> dis(mean_lcp_length);
	
	// Initialize vector with completely random words of specified length
	std::vector<std::string> words(num_words, get_random_word(length));

	// Calculate required signature length to ensure uniqueness
	// The length must satisfy: alphabet_size^signature_length > num_words
	// Using logarithm properties: signature_length > log(num_words)/log(alphabet_size)
	unsigned signature_length = std::ceil(std::log(num_words) / std::log(ALPHABET.size()));

	for (size_t i = 0; i < num_words; ++i)
	{
		auto& word = words[i];

		// Create a unique signature at the end of each word
		// This encodes the index i as a base-ALPHABET.size() number
		size_t i_copy = i;
		for (size_t j = 0; j < signature_length; ++j)
		{
			word[j + length - signature_length] = ALPHABET[i_copy % ALPHABET.size()];
			i_copy /= ALPHABET.size();
		}

		// Generate a random LCP length from the Poisson distribution
		size_t difference = dis(gen);
		
		// Skip modification if the difference would affect the unique signature
		if (difference >= length - signature_length)
		{
			continue;
		}

		// Generate a new random substring to replace part of the word
		// This creates varying common prefix lengths between words
		std::string random_word = get_random_word(length - difference - signature_length);
		std::copy(random_word.begin(), random_word.end(), word.begin() + difference);
	}

	return words;
}
