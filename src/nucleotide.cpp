#include "nucleotide.hpp"

#include <fstream>

bool is_valid_nucleotide(char nucleotide)
{
	return char_to_nucleotide[nucleotide] != Nucleotide::INVALID;
}

std::ostream& operator<<(std::ostream& os, const Nucleotide& nucleotide)
{
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
		throw std::invalid_argument("Invalid nucleotide");
	}

	return os;
}
