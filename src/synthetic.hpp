/**
 * @file synthetic.hpp
 * @brief Provides utilities for generating synthetic string datasets for benchmark testing.
 * @details This file contains functions for generating random strings with controllable 
 * properties such as length, quantity, and shared prefix length. These synthetic datasets
 * are used for benchmarking various trie implementations under controlled conditions.
 * 
 * The implementation uses a diverse alphabet containing lowercase letters, uppercase letters, 
 * digits, and special characters to generate strings with randomized content.
 * 
 * @see BitString.cuh
 * @see SkipTrie.hpp
 * @see ParallelSkipTrie.cuh
 */

#pragma once

#include <string>
#include <vector>

/**
 * @brief Generates a random string of specified length.
 * @details Creates a string using characters from a predefined alphabet containing 
 * lowercase letters, uppercase letters, digits, and special characters. Each character
 * is selected independently using a uniform random distribution.
 * 
 * The implementation uses C++'s standard random number facilities (std::mt19937) 
 * with either a random seed or an optional fixed seed for reproducible results.
 *
 * @param length The desired length of the random string.
 * @return std::string A newly generated random string.
 * @see get_random_words()
 */
std::string get_random_word(size_t length) noexcept;

/**
 * @brief Generates a collection of random strings with controlled commonality.
 * @details Creates a set of random strings where each string has a specified length,
 * and pairs of strings are likely to share a common prefix of a length drawn from
 * a Poisson distribution with the specified mean. To ensure uniqueness, each string
 * includes a unique signature at the end.
 *
 * The generated dataset has the following specific properties:
 *  1. Each string has the exact same length (specified by 'length' parameter)
 *  2. Each string has a unique signature at the end to ensure uniqueness
 *  3. Strings share common prefixes with lengths following a Poisson distribution
 *
 * The algorithm first generates random words, then adds unique signatures, and finally
 * creates controlled common prefixes between the strings. The unique signatures are
 * created by encoding the string's index as a base-N number, where N is the alphabet size.
 * The length of this signature is calculated to ensure all strings can have a unique value.
 * 
 * This function is particularly useful for benchmarking trie data structures with
 * controlled LCP (Longest Common Prefix) distributions.
 *
 * @param length The length of each string to generate.
 * @param num_words The number of strings to generate.
 * @param mean_lcp_length The mean length of the Longest Common Prefix (LCP) between strings.
 * @return std::vector<std::string> A vector containing the generated strings.
 * @see get_random_word()
 */
std::vector<std::string> get_random_words(size_t length, size_t num_words, double mean_lcp_length) noexcept;
