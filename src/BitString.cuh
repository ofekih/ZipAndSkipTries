#pragma once

#include <array>
#include <compare> // for spaceship <=> operator
#include <cstddef> // for size_t
#include <cstdint> // for uintmax_t
#include <numeric> // for std::countl_zero()
#include <vector>

// For file I/O
#include <fstream>
#include <string>

#include <iostream>
#include <stdio.h>

/**
 * @brief A class to store strings compactly into integers, packing as many characters that fit into a single word.
 * 
 * This class allows for compact storage of character strings by packing them into integer arrays, with a focus
 * on efficiency in storage and comparison operations. It leverages bitwise operations for comparisons,
 * optimizing the process by identifying the longest common prefix between two strings.
 * 
 * @tparam CHAR_T Character type to be stored.
 * @tparam CHAR_SIZE_BITS Size of relevant data of CHAR_T in bits.
 */
template<typename CHAR_T, unsigned CHAR_SIZE_BITS = sizeof(CHAR_T) * 8>
class BitString
{
public:
	static constexpr unsigned WORD_SIZE_BITS = sizeof(uintmax_t) * 8; ///< Size of storage word in bits.
	static constexpr unsigned ALPHA = WORD_SIZE_BITS / CHAR_SIZE_BITS; ///< Number of characters per word.

	/**
	 * @brief Constructs a BitString object from a raw character array.
	 * 
	 * @param data Pointer to the beginning of the character array.
	 * @param size Length of the character array.
	 */
	explicit BitString(const CHAR_T* data, size_t size);

	/**
	 * @brief Constructs a BitString object from a container of characters.
	 * 
	 * The container should support size() and operator[].
	 * 
	 * @tparam Container Type of the container.
	 * @param data The container holding characters.
	 */
	template<typename Container>
	explicit BitString(const Container& data);

	/**
	 * @brief Constructs a BitString object from a given string.
	 *
	 * @param data The string used to initialize the BitString object.
	 */
	explicit BitString(const std::string& data)
		: BitString(data.c_str(), data.size())
	{
	}

	/**
	 * @brief Represents a BitString.
	 *
	 * This class represents a BitString, which is a sequence of bits.
	 * The default constructor is only intended to be used by the `from_file` function.
	 */
	BitString() = default;

	void push_back(CHAR_T c) noexcept;

	/**
	 * @brief Accesses a character at a given index.
	 * 
	 * @param index The zero-based index of the character to access.
	 * @return CHAR_T The character at the specified index.
	 */
	CHAR_T operator[](size_t index) const;

	size_t size() const noexcept { return m_size; } ///< Returns the size of the stored string.
	size_t word_count() const noexcept { return m_data.size(); } ///< Returns the number of words used to store the string.

	void clear() noexcept
	{
		m_data.clear();
		m_size = 0;
	}

	struct ResultLCP
	{
		std::strong_ordering result;
		size_t lcp;
	};

	ResultLCP compare(const BitString& other, size_t lcp, size_t max_compare) const noexcept;
	ResultLCP compare(const BitString& other, size_t lcp) const noexcept
	{
		return compare(other, lcp, size());
	}

	const uintmax_t* data() const noexcept { return m_data.data(); }

private:
	std::vector<uintmax_t> m_data; ///< Storage for the compacted string.
	size_t m_size = 0; ///< Size of the string in characters.

	static constexpr std::array<unsigned, ALPHA> GET_SHIFTS(); ///< Generates bit shifts for character extraction.
	static constexpr std::array<unsigned, ALPHA> m_shifts = GET_SHIFTS(); ///< Holds the bit shifts for character extraction.
	static constexpr std::array<uintmax_t, ALPHA> GET_MASKS(); ///< Generates bit masks for character extraction.
	static constexpr std::array<uintmax_t, ALPHA> m_masks = GET_MASKS(); ///< Holds the bit masks for character extraction.
	
	static constexpr CHAR_T GET_CHAR(uintmax_t word, unsigned bit_index); ///< Extracts a character from a word.

public:
	/**
	 * @brief Overloaded stream insertion operator for BitString objects.
	 * 
	 * This function allows BitString objects to be printed to an output stream using the << operator.
	 * 
	 * @tparam CHAR_T The character type used in the BitString.
	 * @tparam CHAR_SIZE_BITS The number of bits in each character.
	 * @param os The output stream to write to.
	 * @param bs The BitString object to be printed.
	 * @return A reference to the output stream after the BitString has been written.
	 */
	friend std::ostream& operator<< <CHAR_T, CHAR_SIZE_BITS>(std::ostream& os, const BitString& bs);

	/**
	 * Prints the bytes of the BitString to the specified output stream.
	 *
	 * @param os The output stream to print the bytes to.
	 */
	void print_bytes(std::ostream& os) const;

	/**
	 * Writes the contents of the BitString to a file in binary format.
	 *
	 * @param filename The name of the file to write the BitString to.
	 */
	void to_file(const std::string& filename) const;

	void to_file(std::ofstream& file) const;

	/**
	 * Reads a BitString from a file in binary format.
	 *
	 * @param filename The name of the file to read from.
	 * @return The BitString read from the file.
	 */
	static BitString from_file(const std::string& filename);

	static BitString from_file(std::ifstream& file);

	std::string to_string() const noexcept;
public:
	class Iterator
	{
	public:
		Iterator(const BitString& bs, size_t index)
			: m_bs(bs), m_index(index)
		{
		}

		CHAR_T operator*() const noexcept
		{
			return m_bs[m_index];
		}

		Iterator& operator++() noexcept
		{
			++m_index;
			return *this;
		}

		bool operator==(const Iterator& other) const noexcept
		{
			return m_index == other.m_index;
		}

	private:
		const BitString& m_bs;
		size_t m_index;
	};

	Iterator begin() const noexcept
	{
		return Iterator(*this, 0);
	}

	Iterator end() const noexcept
	{
		return Iterator(*this, size());
	}
};

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
constexpr std::array<unsigned, BitString<CHAR_T, CHAR_SIZE_BITS>::ALPHA> BitString<CHAR_T, CHAR_SIZE_BITS>::GET_SHIFTS()
{
	std::array<unsigned, ALPHA> shifts;
	for (unsigned i = 0; i < ALPHA; ++i)
	{
		shifts[ALPHA - i - 1] = i * CHAR_SIZE_BITS;
	}

	return shifts;
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
constexpr std::array<uintmax_t, BitString<CHAR_T, CHAR_SIZE_BITS>::ALPHA> BitString<CHAR_T, CHAR_SIZE_BITS>::GET_MASKS()
{
	uintmax_t mask = 0;
	for (unsigned i = 0; i < CHAR_SIZE_BITS; ++i)
	{
		mask |= static_cast<uintmax_t>(1) << i;
	}

	std::array<uintmax_t, ALPHA> masks;
	for (unsigned i = 0; i < ALPHA; ++i)
	{
		masks[i] = mask << m_shifts[i];
	}

	return masks;
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
constexpr CHAR_T BitString<CHAR_T, CHAR_SIZE_BITS>::GET_CHAR(uintmax_t word, unsigned bit_index)
{
	const uintmax_t masked_word = word & m_masks[bit_index];
	const uintmax_t shifted_word = masked_word >> m_shifts[bit_index];

	return static_cast<CHAR_T>(shifted_word);
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
BitString<CHAR_T, CHAR_SIZE_BITS>::BitString(const CHAR_T* data, size_t data_size)
	: m_size(data_size)
{
	const size_t word_count = (size() + ALPHA - 1) / ALPHA;
	m_data.resize(word_count, 0);

	for (size_t i = 0; i < size(); ++i)
	{
		const size_t word_index = i / ALPHA;
		const size_t bit_index = i % ALPHA;
		const uintmax_t mask = static_cast<uintmax_t>(data[i]) << m_shifts[bit_index];
		m_data[word_index] |= mask;
	}
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
template<typename Container>
BitString<CHAR_T, CHAR_SIZE_BITS>::BitString(const Container& data)
	: m_size(data.size())
{
	const size_t word_count = (size() + ALPHA - 1) / ALPHA;
	m_data.resize(word_count, 0);

	for (size_t i = 0; i < size(); ++i)
	{
		const size_t word_index = i / ALPHA;
		const size_t bit_index = i % ALPHA;
		const uintmax_t mask = static_cast<uintmax_t>(data[i]) << m_shifts[bit_index];
		m_data[word_index] |= mask;
	}
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
void BitString<CHAR_T, CHAR_SIZE_BITS>::push_back(CHAR_T c) noexcept
{
	const size_t word_index = size() / ALPHA;
	const size_t bit_index = size() % ALPHA;
	const uintmax_t mask = static_cast<uintmax_t>(c) << m_shifts[bit_index];

	if (word_index >= m_data.size())
	{
		m_data.push_back(0);
	}

	m_data[word_index] |= mask;
	++m_size;
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
CHAR_T BitString<CHAR_T, CHAR_SIZE_BITS>::operator[](size_t index) const
{
	return GET_CHAR(m_data[index / ALPHA], index % ALPHA);
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
auto BitString<CHAR_T, CHAR_SIZE_BITS>::compare(const BitString& other, size_t lcp, size_t max_compare) const noexcept -> ResultLCP
{
	const size_t min_size = std::min(size(), other.size());
	const size_t max_compare_until = lcp + max_compare;
	size_t word_index = lcp / ALPHA;

	lcp -= lcp % ALPHA; // Reset to word boundary

	while (lcp < min_size && lcp < max_compare_until)
	{
		const uintmax_t XOR = m_data[word_index] ^ other.m_data[word_index];

		if (XOR == 0)
		{
			lcp += ALPHA;
			++word_index;
			continue;
		}

		const size_t l_zero = std::countl_zero(XOR);

		lcp += l_zero / CHAR_SIZE_BITS;

		const size_t bit_index = lcp % ALPHA;

		const CHAR_T lhs = GET_CHAR(m_data[word_index], bit_index);
		const CHAR_T rhs = GET_CHAR(other.m_data[word_index], bit_index);

		return { lhs <=> rhs, lcp };
	}

	if (lcp >= min_size)
	{
		return { size() <=> other.size(), min_size };
	}

	return { std::strong_ordering::equal, lcp };
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
void BitString<CHAR_T, CHAR_SIZE_BITS>::print_bytes(std::ostream& os) const
{
	static constexpr unsigned shift = WORD_SIZE_BITS - 8;

	for (size_t i = 0, a = 0; i + ALPHA <= size(); i += ALPHA, ++a)
	{
		uintmax_t word = m_data[a];
		
		for (unsigned c = 0; c < WORD_SIZE_BITS / 8; c++, word <<= 8)
		{
			os << static_cast<uint8_t>(word >> shift);
		}
	}

	uintmax_t word = m_data.back();
	unsigned remaining_bits = (size() % ALPHA) * CHAR_SIZE_BITS;

	for (; remaining_bits > 0; remaining_bits -= 8, word <<= 8)
	{
		os << static_cast<uint8_t>(word >> shift);
	}
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
std::ostream& operator<<(std::ostream& os, const BitString<CHAR_T, CHAR_SIZE_BITS>& bs)
{
	for (size_t i = 0; i < bs.size(); ++i)
	{
		os << bs[i];
	}

	return os;
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
void BitString<CHAR_T, CHAR_SIZE_BITS>::to_file(std::ofstream& file) const
{
	file.write(reinterpret_cast<const char*>(&m_size), sizeof(size_t));
	file.write(reinterpret_cast<const char*>(m_data.data()), m_data.size() * sizeof(uintmax_t));
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
void BitString<CHAR_T, CHAR_SIZE_BITS>::to_file(const std::string& filename) const
{
	std::ofstream file(filename, std::ios::binary);
	if (!file)
	{
		throw std::runtime_error("Failed to open file for writing.");
	}

	to_file(file);
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
BitString<CHAR_T, CHAR_SIZE_BITS> BitString<CHAR_T, CHAR_SIZE_BITS>::from_file(std::ifstream& file)
{
	BitString<CHAR_T, CHAR_SIZE_BITS> bs;

	file.read(reinterpret_cast<char*>(&bs.m_size), sizeof(size_t));

	const size_t word_count = (bs.size() + ALPHA - 1) / ALPHA;
	bs.m_data.resize(word_count, 0);
	file.read(reinterpret_cast<char*>(bs.m_data.data()), bs.m_data.size() * sizeof(uintmax_t));

	return bs;
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
BitString<CHAR_T, CHAR_SIZE_BITS> BitString<CHAR_T, CHAR_SIZE_BITS>::from_file(const std::string& filename)
{
	std::ifstream file(filename, std::ios::binary);
	if (!file)
	{
		throw std::runtime_error("Failed to open file for reading.");
	}

	return from_file(file);
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
std::string BitString<CHAR_T, CHAR_SIZE_BITS>::to_string() const noexcept
{
	std::string str;
	str.reserve(size());

	for (size_t i = 0; i < size(); ++i)
	{
		str.push_back((*this)[i]);
	}

	return str;
}