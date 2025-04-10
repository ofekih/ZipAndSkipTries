/**
 * @file nucleotide.hpp
 * @brief Defines nucleotide representation and related utility functions.
 * @details This file provides a type-safe representation of DNA nucleotides and
 * utility functions for converting between character and nucleotide formats.
 * It supports case-insensitive conversion and implements streaming operators for
 * convenient output. This is used primarily for genetic sequence processing in
 * the trie data structures.
 * 
 * @see genetics.hpp
 * @see BitString.cuh
 */

#pragma once

#include <array>
#include <cstdint> // for uint8_t
#include <stdexcept> // for std::invalid_argument
#include <fstream>

/**
 * @enum Nucleotide
 * @brief Enumerates the four nucleotides of DNA.
 * 
 * This enumeration provides a type-safe way of representing the four nucleotides
 * in DNA: Adenine (A), Cytosine (C), Guanine (G), and Thymine (T). It also includes
 * an INVALID value for error handling.
 */
enum class Nucleotide : uint8_t
{
	A = 0, ///< Adenine
	C = 1, ///< Cytosine
	G = 2, ///< Guanine
	T = 3,  ///< Thymine
	INVALID ///< Invalid nucleotide
};

/**
 * @brief Lookup table for converting characters to nucleotides.
 * @details This constexpr array provides O(1) conversion from ASCII characters to Nucleotide enum values.
 * It supports both uppercase and lowercase letters (A/a, C/c, G/g, T/t) and marks all other
 * characters as INVALID. Using a lookup table avoids branching in critical code paths.
 */
static constexpr std::array<Nucleotide, 256> char_to_nucleotide = []() {
	std::array<Nucleotide, 256> result;
	result.fill(Nucleotide::INVALID);

	// Initialize valid nucleotide mappings
	result['A'] = Nucleotide::A;
	result['a'] = Nucleotide::A;
	result['C'] = Nucleotide::C;
	result['c'] = Nucleotide::C;
	result['G'] = Nucleotide::G;
	result['g'] = Nucleotide::G;
	result['T'] = Nucleotide::T;
	result['t'] = Nucleotide::T;

	return result;
}();

/**
 * @brief Lookup table for converting nucleotides to characters.
 * @details This array provides O(1) conversion from Nucleotide enum values to their
 * corresponding ASCII character representation. The array is indexed by the underlying
 * uint8_t value of the Nucleotide enum, so Nucleotide::A (0) maps to 'A', etc.
 * Note that this array doesn't include a mapping for Nucleotide::INVALID.
 * @see Nucleotide
 */
static constexpr std::array<char, 4> nucleotide_to_char = {
	'A', 'C', 'G', 'T'
};

/**
 * @brief Checks if a character represents a valid nucleotide.
 * @details Determines whether the provided character is a valid nucleotide representation
 * (A, C, G, T in either uppercase or lowercase).
 * 
 * @param nucleotide The character to check.
 * @return true if the character is a valid nucleotide representation, false otherwise.
 * @see char_to_nucleotide
 */
bool is_valid_nucleotide(char nucleotide);

/**
 * @brief Inserts a textual representation of a nucleotide into an output stream.
 * @details Converts the Nucleotide enum value to its character representation and
 * inserts it into the provided output stream. Uses a switch statement to handle
 * each nucleotide type explicitly and throws an std::invalid_argument exception
 * for Nucleotide::INVALID rather than silently outputting a placeholder.
 * This allows Nucleotide objects to be used directly with stream operators.
 * 
 * @param os The output stream to insert into.
 * @param nucleotide The nucleotide to insert.
 * @return Reference to the output stream for chaining.
 * @throws std::invalid_argument if nucleotide is Nucleotide::INVALID.
 * @see nucleotide_to_char
 */
std::ostream& operator<<(std::ostream& os, const Nucleotide& nucleotide);
