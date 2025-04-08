#pragma once

#include <array>
#include <cstdint> // for uint8_t
#include <stdexcept> // for std::invalid_argument

#include <fstream>

/**
 * @enum NucleotideEnum
 * @brief Enumerates the four nucleotides of DNA.
 * 
 * This enumeration provides a type-safe way of representing the four nucleotides
 * in DNA: Adenine (A), Cytosine (C), Guanine (G), and Thymine (T).
 */
enum class Nucleotide : uint8_t
{
	A = 0, ///< Adenine
	C = 1, ///< Cytosine
	G = 2, ///< Guanine
	T = 3,  ///< Thymine
	INVALID ///< Invalid nucleotide
};

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

static constexpr std::array<char, 4> nucleotide_to_char = {
	'A', 'C', 'G', 'T'
};

bool is_valid_nucleotide(char nucleotide);

std::ostream& operator<<(std::ostream& os, const Nucleotide& nucleotide);
