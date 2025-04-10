/**
 * @file nucleotide.cpp
 * @brief Implementation of nucleotide utility functions.
 * @details This file implements the functions declared in nucleotide.hpp for
 * nucleotide validation and stream output operations. It provides efficient
 * implementations for handling nucleotide conversions and validation.
 * 
 * @see nucleotide.hpp
 */
#include "nucleotide.hpp"

#include <fstream>

bool is_valid_nucleotide(char nucleotide)
{
	// Use the lookup table approach for efficient character validation
	// This avoids multiple comparisons by directly indexing into the table
	return char_to_nucleotide[nucleotide] != Nucleotide::INVALID;
}

std::ostream& operator<<(std::ostream& os, const Nucleotide& nucleotide)
{
	// Handle each nucleotide case explicitly in the switch statement
	// Note: We could use the nucleotide_to_char array directly for valid nucleotides,
	// but the switch approach provides clearer error handling for invalid values
	switch (nucleotide)
	{
	case Nucleotide::A:
		os << 'A';
		break;
	case Nucleotide::C:
		os << 'C';
		break;
	case Nucleotide::G:
		os << 'G';
		break;
	case Nucleotide::T:
		os << 'T';
		break;
	default:
		// Throw an exception for invalid nucleotides rather than silently outputting a placeholder
		// This helps catch programming errors early
		throw std::invalid_argument("Invalid nucleotide");
	}

	return os;
}
