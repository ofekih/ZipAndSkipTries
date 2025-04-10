/**
 * @file synthetic.cpp
 * @brief Implementation of synthetic data generation utilities for benchmark testing.
 *
 * @details This file implements the functions declared in synthetic.hpp for generating
 * random strings with controlled properties. The implementation provides two key functions:
 *
 * 1. get_random_word: Generates a single random string of specified length using
 *    characters from a diverse alphabet.
 *
 * 2. get_random_words: Generates a collection of random strings with controlled
 *    Longest Common Prefix (LCP) lengths following a Poisson distribution.
 *
 * The implementation uses C++ standard library random number generators (std::mt19937)
 * and distributions (std::uniform_int_distribution, std::poisson_distribution) to
 * create randomized but statistically controlled datasets. Each generated string
 * includes a unique signature to ensure distinctness while maintaining the desired
 * LCP distribution properties.
 *
 * These synthetic datasets are essential for benchmarking trie data structures under
 * controlled conditions, allowing for systematic performance analysis with varying
 * input characteristics.
 *
 * @see synthetic.hpp
 * @see get_random_word
 * @see get_random_words
 */

#include "synthetic.hpp"

#include <algorithm>
#include <random>
#include <string>
#include <vector>

/**
 * @brief The alphabet used for generating random strings.
 * @details Contains 94 characters including lowercase letters (a-z), uppercase letters (A-Z),
 * digits (0-9), and special characters. This diverse alphabet provides a large character space
 * for generating random strings, which helps create realistic test data with high entropy.
 * The size of this alphabet also affects the calculation of unique signature length in
 * get_random_words().
 */
static const std::string ALPHABET = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!\"#$%&'()*+,-./:;<=>?@[\\]^_`{|}~";

/**
 * @brief Implementation of get_random_word function.
 * @details Generates a random string of specified length using characters from the ALPHABET.
 * The implementation uses a static random number generator (std::mt19937) initialized with
 * either a random seed from std::random_device or an optional fixed seed for reproducible results.
 *
 * Each character in the string is selected independently using a uniform distribution over
 * the entire alphabet, ensuring equal probability for all characters.
 *
 * @param length The desired length of the random string.
 * @return std::string A newly generated random string of the specified length.
 * @see synthetic.hpp for the function declaration and additional documentation.
 */
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

/**
 * @brief Implementation of get_random_words function.
 * @details Generates a collection of random strings with controlled Longest Common Prefix (LCP)
 * properties. The implementation follows these steps:
 *
 * 1. Generate num_words random strings of the specified length.
 * 2. Calculate the minimum signature length needed to ensure uniqueness.
 * 3. Add a unique signature at the end of each string by encoding its index.
 * 4. For each string, generate a random LCP length from a Poisson distribution.
 * 5. Replace part of each string to create the desired LCP pattern while preserving the unique signature.
 *
 * The Poisson distribution ensures that the LCP lengths follow a realistic distribution with
 * the specified mean, while the unique signatures guarantee that all strings remain distinct
 * regardless of the LCP modifications.
 *
 * @param length The length of each string to generate.
 * @param num_words The number of strings to generate.
 * @param mean_lcp_length The mean length of the Longest Common Prefix between strings.
 * @return std::vector<std::string> A vector containing the generated strings.
 * @see synthetic.hpp for the function declaration and additional documentation.
 */
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
