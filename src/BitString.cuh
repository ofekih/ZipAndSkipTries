/**
 * @file BitString.cuh
 * @brief Defines the BitString class for compact character string storage and comparison.
 *
 * This file contains the definition of the BitString template class, which provides
 * mechanisms for storing sequences of characters efficiently by packing them into
 * integer words. It also includes methods for comparing BitString objects, including
 * sequential and parallel (CUDA-based) comparison routines.
 */
#pragma once

#include <array>
#include <bit>       // for std::countl_zero, std::bit_ceil
#include <compare>   // for std::strong_ordering and <=> operator
#include <cstddef>   // for size_t
#include <cstdint>   // for uintmax_t
#include <numeric>   // Provides std::min (for two arguments)
#include <vector>
#include <algorithm> // Include for std::min/max

// For file I/O
#include <fstream>
#include <stdexcept> // for std::runtime_error
#include <string>

#include <iostream>

// Include CUDA headers for parallel operations
#include "msw.cuh"
#include "utility.cuh"

/**
 * @brief Minimum number of words required to trigger parallel comparison.
 * @note If the number of words to compare is less than or equal to this value,
 * sequential comparison (`seq_k_compare`) will be used instead of
 * parallel comparison (`par_k_compare`). Set to 1 to always prefer parallel
 * comparison when applicable.
 * @warning Should not be set to 0.
 */
static constexpr size_t MIN_PAR_COMPARE_WORD_SIZE = 1;

/**
 * @brief A class to store strings compactly into integers, packing multiple characters into a single word.
 *
 * This class allows for space-efficient storage of character strings by packing them into an array
 * of `uintmax_t` integers. It focuses on optimizing storage and comparison operations.
 * Comparisons leverage bitwise operations and can identify the longest common prefix (LCP)
 * efficiently. Parallel comparison using CUDA is supported for large strings.
 *
 * @tparam CHAR_T The character type to be stored (e.g., `char`, `uint8_t`, `wchar_t`).
 * @tparam CHAR_SIZE_BITS The number of bits representing relevant data within `CHAR_T`.
 * Defaults to the full size of `CHAR_T` in bits. This allows handling
 * types where only a subset of bits are meaningful (e.g., 7-bit ASCII in an 8-bit char,
 * or 2-bit nucleotides for DNA strings).
 */
template<typename CHAR_T, unsigned CHAR_SIZE_BITS = sizeof(CHAR_T) * 8>
class BitString
{
public:
	/** @brief Size of the underlying storage word (`uintmax_t`) in bits. */
	static constexpr unsigned WORD_SIZE_BITS = sizeof(uintmax_t) * 8;
	/** @brief Number of characters that can be packed into a single storage word. */
	static constexpr unsigned ALPHA = WORD_SIZE_BITS / CHAR_SIZE_BITS;
	/** @brief Minimum number of characters required to trigger parallel comparison. Derived from `MIN_PAR_COMPARE_WORD_SIZE`. */
	static constexpr size_t MIN_PAR_COMPARE_CHAR_SIZE = MIN_PAR_COMPARE_WORD_SIZE * ALPHA;

	/**
	 * @brief Constructs a BitString object from a raw character array.
	 * @details Packs the input characters into the internal `uintmax_t` vector.
	 * @param data Pointer to the beginning of the character array. Must not be null if size > 0.
	 * @param size The number of characters in the array.
	 */
	explicit BitString(const CHAR_T* data, size_t size);

	/**
	 * @brief Constructs a BitString object from a container of characters.
	 * @details Packs characters from the container into the internal vector.
	 *
	 * The container must support `size()` and random access via `operator[]`.
	 * Examples include `std::vector<CHAR_T>`, `std::basic_string<CHAR_T>`, `std::array<CHAR_T, N>`.
	 *
	 * @tparam Container Type of the container holding the character data.
	 * @param data The container instance.
	 */
	template<typename Container>
	explicit BitString(const Container& data);

	/**
	 * @brief Constructs a BitString object from a `std::string`.
	 * @note This constructor assumes `CHAR_T` is compatible with `char`.
	 * @param data The `std::string` used to initialize the BitString object.
	 */
	explicit BitString(const std::string& data)
		: BitString(reinterpret_cast<const CHAR_T*>(data.c_str()), data.size()) // Reinterpret cast assumes CHAR_T layout matches char
	{
		// Body intentionally empty
	}

	/**
	 * @brief Default constructor.
	 * @note Creates an empty BitString. Primarily intended for use by the `from_file` static method.
	 */
	BitString() = default;

	/**
	 * @brief Appends a single character to the end of the BitString.
	 * @details Adds the character to the last word, potentially allocating a new word if needed.
	 * @param c The character to append.
	 * @note This operation may reallocate the underlying storage if capacity is exceeded.
	 */
	void push_back(CHAR_T c) noexcept;

	/**
	 * @brief Accesses the character at a specific index.
	 * @details Calculates the word and sub-index, then uses GET_CHAR for extraction. Behavior undefined if index is out of bounds.
	 *
	 * Behavior is undefined for out-of-bounds access.
	 *
	 * @param index The zero-based index of the character to access. Must be less than `size()`.
	 * @return CHAR_T The character at the specified index.
	 * @note This operation involves bitwise extraction and may be slower than direct array access.
	 */
	CHAR_T operator[](size_t index) const;

	/**
	 * @brief Returns the number of characters stored in the BitString.
	 * @return size_t The total number of characters.
	 */
	size_t size() const noexcept { return m_size; }

	/**
	 * @brief Returns the number of `uintmax_t` words used for storage.
	 * @return size_t The number of words in the internal data vector.
	 */
	size_t word_count() const noexcept { return m_data.size(); }

	/**
	 * @brief Clears the contents of the BitString, making it empty.
	 * @note Resets size to 0 and clears the internal data vector.
	 */
	void clear() noexcept
	{
		m_data.clear();
		m_size = 0;
	}

	/**
	 * @brief Structure to hold the result of a comparison, including the ordering and the length of the common prefix.
	 */
	struct ResultLCP
	{
		std::strong_ordering result; ///< The comparison result (less, greater, equal).
		size_t lcp;                  ///< The length of the longest common prefix in characters.
	};

	/**
	 * @brief Performs sequential comparison with another BitString, starting from a given LCP offset, up to a maximum length.
	 * @details Iterates word by word, finds the first mismatch using XOR and countl_zero.
	 *
	 * Compares this BitString with `other`, assuming the first `lcp` characters are already known to be equal.
	 * The comparison proceeds character by character (packed within words) sequentially, but stops strictly
	 * at `lcp + max_compare` characters or the end of the shorter string, whichever comes first.
	 *
	 * @param other The BitString to compare against.
	 * @param lcp The starting offset (longest common prefix length) for the comparison, in characters.
	 * @param max_compare The maximum number of additional characters to compare beyond the initial `lcp`.
	 * The comparison will not proceed beyond `lcp + max_compare` total characters from the start.
	 * @return ResultLCP A struct containing the comparison result (`std::strong_ordering`) of the substrings
	 * up to the comparison limit and the final calculated LCP length (which will not exceed the limit).
	 * Returns `std::strong_ordering::equal` if the substrings match up to the limit, even if the
	 * full strings might differ later.
	 * @note This comparison is performed entirely on the CPU. Uses `std::countl_zero`.
	 */
	ResultLCP seq_k_compare(const BitString& other, size_t lcp, size_t max_compare) const noexcept;

	/**
	 * @brief Performs sequential comparison with another BitString, starting from a given LCP offset.
	 *
	 * Compares this BitString with `other`, assuming the first `lcp` characters are equal.
	 * Compares up to the end of the shorter string.
	 *
	 * @param other The BitString to compare against.
	 * @param lcp The starting offset (longest common prefix length) for the comparison, in characters.
	 * @return ResultLCP A struct containing the comparison result and the final LCP length.
	 * @note This comparison is performed entirely on the CPU. Equivalent to `seq_k_compare(other, lcp, size())`.
	 */
	ResultLCP seq_k_compare(const BitString& other, size_t lcp) const noexcept
	{
		return seq_k_compare(other, lcp, size());
	}

	/**
	 * @brief Performs parallel comparison (potentially using CUDA) with another BitString, starting from a given LCP offset, up to a maximum length.
	 * @details Checks size, copies data to GPU if necessary, calls parallel mismatch kernel, then refines LCP. Falls back to sequential if too small.
	 *
	 * Compares this BitString with `other`, assuming the first `lcp` characters are equal.
	 * If the number of words to compare exceeds `MIN_PAR_COMPARE_WORD_SIZE`, this function attempts
	 * to use a parallel CUDA kernel (`par_find_mismatch_s`) for faster mismatch detection.
	 * Otherwise, it falls back to `seq_k_compare`. The comparison stops strictly at `lcp + max_compare`
	 * characters or the end of the shorter string, whichever comes first.
	 *
	 * @param other The BitString to compare against.
	 * @param lcp The starting offset (longest common prefix length) for the comparison, in characters.
	 * @param max_compare The maximum number of additional characters to compare beyond the initial `lcp`.
	 * The comparison will not proceed beyond `lcp + max_compare` total characters from the start.
	 * @param d_a Pointer to device (GPU) memory allocated for this BitString's data. Data might be copied here if needed.
	 * @param d_largeblock Pointer to auxiliary device (GPU) memory required by the parallel kernel.
	 * @param [in,out] max_copied Tracks the amount of data already copied to `d_a` to avoid redundant transfers. Updated by this function.
	 * @return ResultLCP A struct containing the comparison result (`std::strong_ordering`) of the substrings
	 * up to the comparison limit and the final calculated LCP length (which will not exceed the limit).
	 * Returns `std::strong_ordering::equal` if the substrings match up to the limit, even if the
	 * full strings might differ later.
	 * @see seq_k_compare
	 * @see MIN_PAR_COMPARE_WORD_SIZE
	 * @note Requires properly allocated CUDA device memory (`d_a`, `d_largeblock`) and a valid CUDA context.
	 * @note The `other` BitString's data is assumed to be accessible directly (e.g., host pinned memory or already on device if applicable to `par_find_mismatch_s`).
	 * @note Uses `std::bit_ceil` and `std::countl_zero`.
	 */
	ResultLCP par_k_compare(const BitString& other, size_t lcp, size_t max_compare, uintmax_t* d_a, uintmax_t* d_largeblock, size_t& max_copied) const noexcept;

	/**
	 * @brief Performs parallel comparison (potentially using CUDA) with another BitString, starting from a given LCP offset.
	 *
	 * Compares this BitString with `other`, assuming the first `lcp` characters are equal.
	 * Compares up to the end of the shorter string using parallel methods if applicable.
	 *
	 * @param other The BitString to compare against.
	 * @param lcp The starting offset (longest common prefix length) for the comparison, in characters.
	 * @param d_a Pointer to device (GPU) memory for this BitString's data.
	 * @param d_largeblock Pointer to auxiliary device (GPU) memory.
	 * @param [in,out] max_copied Tracks the amount of data copied to `d_a`.
	 * @return ResultLCP A struct containing the comparison result and the final LCP length.
	 * @note Equivalent to `par_k_compare(other, lcp, size(), d_a, d_largeblock, max_copied)`.
	 * @see par_k_compare(const BitString&, size_t, size_t, uintmax_t*, uintmax_t*, size_t&)
	 */
	ResultLCP par_k_compare(const BitString& other, size_t lcp, uintmax_t* d_a, uintmax_t* d_largeblock, size_t& max_copied) const noexcept
	{
		return par_k_compare(other, lcp, size(), d_a, d_largeblock, max_copied);
	}

	/**
	 * @brief Provides direct access to the underlying packed data array.
	 * @return const uintmax_t* A pointer to the beginning of the internal `uintmax_t` data array. Returns `nullptr` if the BitString is empty.
	 * @warning Modifying the data through this pointer may invalidate the BitString's state. Use with caution.
	 */
	const uintmax_t* data() const noexcept { return m_data.data(); }

private:
	/** @brief Internal storage vector holding the packed character data. */
	std::vector<uintmax_t> m_data;
	/** @brief The number of characters stored in the BitString. */
	size_t m_size = 0;

	/**
	 * @brief Generates the bit shift amounts needed for character extraction/insertion at compile time.
	 * @details Calculates the left shift amount needed for each character position within a word.
	 * @tparam CHAR_T Character type.
	 * @tparam CHAR_SIZE_BITS Bits per character.
	 * @return std::array<unsigned, ALPHA> An array where each element `i` contains the left-shift amount
	 * required to position the character at sub-word index `i` correctly within a `uintmax_t`.
	 */
	static constexpr std::array<unsigned, ALPHA> GET_SHIFTS();
	/** @brief Compile-time generated array of bit shift amounts. */
	static constexpr std::array<unsigned, ALPHA> m_shifts = GET_SHIFTS();

	/**
	 * @brief Generates the bit masks needed for character extraction at compile time.
	 * @details Creates a mask for each character position within a word.
	 * @tparam CHAR_T Character type.
	 * @tparam CHAR_SIZE_BITS Bits per character.
	 * @return std::array<uintmax_t, ALPHA> An array where each element `i` is a mask that isolates
	 * the bits corresponding to the character at sub-word index `i`.
	 */
	static constexpr std::array<uintmax_t, ALPHA> GET_MASKS();
	/** @brief Compile-time generated array of bit masks. */
	static constexpr std::array<uintmax_t, ALPHA> m_masks = GET_MASKS();

	/**
	 * @brief Extracts a single character from a given word at a specific bit index (sub-word position).
	 * @details Applies the appropriate mask and shift to isolate and retrieve the character.
	 * @tparam CHAR_T Character type.
	 * @tparam CHAR_SIZE_BITS Bits per character.
	 * @param word The `uintmax_t` word containing the packed characters.
	 * @param bit_index The index within the word (0 to ALPHA-1) corresponding to the desired character.
	 * @return CHAR_T The extracted character.
	 */
	static constexpr CHAR_T GET_CHAR(uintmax_t word, unsigned bit_index);

public:
	// Forward declaration for the friend function template
	template<typename CT, unsigned CSB>
	friend std::ostream& operator<<(std::ostream& os, const BitString<CT, CSB>& bs);

	/**
	 * @brief Prints the raw byte representation of the BitString to an output stream.
	 * @details Iterates through words and extracts bytes based on `WORD_SIZE_BITS`. Handles the last partial word correctly.
	 *
	 * Writes the packed `uintmax_t` data to the stream byte by byte, respecting the actual
	 * number of characters (`size()`). Handles partial final words correctly.
	 *
	 * @param os The output stream (e.g., `std::cout`, `std::ofstream`) to write to.
	 */
	void print_bytes(std::ostream& os) const;

	/**
	 * @brief Writes the BitString contents to a file in binary format.
	 * @details Opens the file, calls the stream version of `to_file`, and closes the file.
	 *
	 * Serializes the size and the packed data array to the specified file.
	 *
	 * @param filename The path to the output file.
	 * @throw std::runtime_error If the file cannot be opened for writing.
	 */
	void to_file(const std::string& filename) const;

	/**
	 * @brief Writes the BitString contents to an already open output file stream.
	 * @details Writes the size, then the raw `uintmax_t` data. Does not perform error checking on writes.
	 *
	 * Serializes the size and the packed data array to the stream.
	 *
	 * @param file The `std::ofstream` object, opened in binary mode.
	 */
	void to_file(std::ofstream& file) const;

	/**
	 * @brief Reads a BitString from a file in binary format.
	 * @details Opens the file, calls the stream version of `from_file`, and closes the file.
	 *
	 * Deserializes the size and packed data from the specified file.
	 *
	 * @param filename The path to the input file.
	 * @return BitString The BitString object reconstructed from the file data.
	 * @throw std::runtime_error If the file cannot be opened for reading.
	 */
	static BitString from_file(const std::string& filename);

	/**
	 * @brief Reads a BitString from an already open input file stream.
	 * @details Reads the size, calculates words needed, resizes, then reads raw data. Does not perform error checking on reads.
	 *
	 * Deserializes the size and packed data from the stream.
	 *
	 * @param file The `std::ifstream` object, opened in binary mode.
	 * @return BitString The BitString object reconstructed from the stream data.
	 */
	static BitString from_file(std::ifstream& file);

	/**
	 * @brief Converts the BitString back to a standard `std::string`.
	 * @details Iterates through characters using iterators and appends them to a string. Assumes `CHAR_T` is convertible to `char`.
	 * @return std::string A string containing the characters stored in the BitString.
	 * @note This assumes `CHAR_T` is convertible to `char`.
	 */
	std::string to_string() const noexcept;

public:
	/**
	 * @brief An iterator class for traversing the characters within a BitString.
	 *
	 * Provides read-only, forward iteration over the characters stored in the BitString.
	 * Dereferencing yields the character value, and incrementing moves to the next character.
	 */
	class Iterator
	{
	public:
		using iterator_category = std::forward_iterator_tag;
		using value_type = CHAR_T;
		using difference_type = std::ptrdiff_t;
		using pointer = const CHAR_T*;   // Or void if value is returned by copy
		using reference = CHAR_T; // Return by value

		/**
		 * @brief Constructs an iterator pointing to a specific character index.
		 * @param bs The BitString to iterate over.
		 * @param index The initial character index.
		 */
		Iterator(const BitString& bs, size_t index)
			: m_bs(bs), m_index(index)
		{
			// Body intentionally empty
		}

		/**
		 * @brief Dereferences the iterator to get the character at the current position.
		 * @return CHAR_T The character value.
		 */
		CHAR_T operator*() const noexcept
		{
			return m_bs[m_index];
		}

		/**
		 * @brief Pre-increments the iterator to point to the next character.
		 * @return Iterator& A reference to the incremented iterator.
		 */
		Iterator& operator++() noexcept
		{
			++m_index;
			return *this;
		}

		/**
		 * @brief Post-increments the iterator.
		 * @return Iterator A copy of the iterator before incrementing.
		 * @note Prefer pre-increment (`++it`) when possible for efficiency.
		 */
		Iterator operator++(int) noexcept
		{
			Iterator temp = *this;
			++(*this);
			return temp;
		}


		/**
		 * @brief Compares this iterator with another for equality.
		 * @param other The iterator to compare against.
		 * @return bool True if both iterators point to the same index within the same BitString (implicitly).
		 */
		bool operator==(const Iterator& other) const noexcept
		{
			return m_index == other.m_index;
		}

	private:
		/** @brief Reference to the BitString being iterated over. */
		const BitString& m_bs;
		/** @brief Current character index. */
		size_t m_index;
	};

	/**
	 * @brief Returns an iterator pointing to the beginning of the BitString.
	 * @return Iterator An iterator to the first character.
	 */
	Iterator begin() const noexcept
	{
		return Iterator(*this, 0);
	}

	/**
	 * @brief Returns an iterator pointing past the end of the BitString.
	 * @return Iterator An iterator representing the end sentinel.
	 */
	Iterator end() const noexcept
	{
		return Iterator(*this, size());
	}
};

// == Implementation of Template Methods ==

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
	: m_size(data.size()) // Get size from the container
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
auto BitString<CHAR_T, CHAR_SIZE_BITS>::seq_k_compare(const BitString& other, size_t lcp, size_t max_compare) const noexcept -> ResultLCP
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

		// Mismatch found within limit
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
auto BitString<CHAR_T, CHAR_SIZE_BITS>::par_k_compare(const BitString& other, size_t lcp, size_t max_compare, uintmax_t* d_a, uintmax_t* d_largeblock, size_t& max_copied) const noexcept -> ResultLCP
{
	// Fallback for small comparisons
	if (max_compare <= MIN_PAR_COMPARE_CHAR_SIZE)
	{
		return seq_k_compare(other, lcp, max_compare);
	}

	const size_t min_size = std::min(size(), other.size());
	const size_t min_word_count = std::min(word_count(), other.word_count());
	size_t word_index = lcp / ALPHA;

	// Calculate how many full words remain to potentially compare
	const size_t num_remaining_words = (min_word_count > word_index) ? (min_word_count - word_index) : 0;

	// Calculate max words to compare based on max_compare characters
	// Need to compare up to the word containing character (lcp + max_compare - 1)
	const size_t end_char_index = lcp + max_compare;
	const size_t end_word_index = (end_char_index + ALPHA - 1) / ALPHA; // Word index containing the last char to check
	const size_t max_compare_words = (end_word_index > word_index) ? (end_word_index - word_index) : 0;


	const size_t actual_max_compare_words = std::min(max_compare_words, num_remaining_words);

	// Fallback if the actual number of words to compare is too small
	if (actual_max_compare_words <= MIN_PAR_COMPARE_WORD_SIZE)
	{
		return seq_k_compare(other, lcp, max_compare);
	}

	// Ensure data is on the device
	if (max_copied < word_index + actual_max_compare_words)
	{
		size_t end_copy_index = word_index + actual_max_compare_words;
		
		// get next power of 2 greater than or equal to end_copy_index
		size_t end_copy = std::min(std::bit_ceil(end_copy_index), word_count());
		if (end_copy > max_copied) {
			copy_to_device(d_a + max_copied, data() + max_copied, end_copy - max_copied);
			max_copied = end_copy;
		}
	}

	lcp -= lcp % ALPHA; // Reset LCP to the start of the current word boundary for parallel comparison

	// Perform parallel mismatch search
	auto msw = par_find_mismatch_s(d_a + word_index, other.data() + word_index, d_largeblock, actual_max_compare_words);

	// Calculate LCP based on mismatch word index (msw)
	lcp += msw * ALPHA;

	// If no mismatch found within the compared words
	if (msw == actual_max_compare_words)
	{
		// Check if the comparison reached the end of either string
		if (lcp >= min_size)
		{
			return { size() <=> other.size(), min_size };
		}

		return { std::strong_ordering::equal, lcp };
	}

	// Mismatch found by parallel kernel. Refine LCP within the mismatch word.
	word_index = lcp / ALPHA; // Update word_index to the mismatch word

	// Find exact mismatch character within the word
	const uintmax_t XOR = m_data[word_index] ^ other.m_data[word_index];
	const size_t l_zero = std::countl_zero(XOR);

	lcp += l_zero / CHAR_SIZE_BITS;

	const size_t bit_index = lcp % ALPHA;

	// Mismatch is within the limit
	const CHAR_T lhs = GET_CHAR(m_data[word_index], bit_index);
	const CHAR_T rhs = GET_CHAR(other.m_data[word_index], bit_index);

	return { lhs <=> rhs, lcp };
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

	// Handle last partial word if it exists
	if (size() % ALPHA != 0 && !m_data.empty())
	{
		uintmax_t word = m_data.back();
		unsigned remaining_chars = size() % ALPHA;
		unsigned remaining_bits = remaining_chars * CHAR_SIZE_BITS;
		unsigned bytes_to_print = (remaining_bits + 7) / 8; // Round up bytes

		// Process the required bytes from the most significant side
		for (unsigned byte_idx = 0; byte_idx < bytes_to_print; ++byte_idx)
		{
			os << static_cast<uint8_t>(word >> (shift - byte_idx * 8));
		}
	}
}


/**
 * @brief Overloaded stream insertion operator for printing BitString content.
 * @details Iterates through the BitString using its iterators and prints each character.
 * @tparam CHAR_T Character type.
 * @tparam CHAR_SIZE_BITS Bits per character.
 * @param os The output stream (e.g., `std::cout`).
 * @param bs The BitString object to print.
 * @return std::ostream& Reference to the output stream.
 */
template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
std::ostream& operator<<(std::ostream& os, const BitString<CHAR_T, CHAR_SIZE_BITS>& bs)
{
	for (CHAR_T c : bs)
	{
		os << c;
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
		throw std::runtime_error("Failed to open file for writing: " + filename);
	}
	to_file(file);
	// file automatically closed by destructor
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
BitString<CHAR_T, CHAR_SIZE_BITS> BitString<CHAR_T, CHAR_SIZE_BITS>::from_file(std::ifstream& file)
{
	BitString<CHAR_T, CHAR_SIZE_BITS> bs;

	file.read(reinterpret_cast<char*>(&bs.m_size), sizeof(size_t));
	if (!file) { // Check read success
		throw std::runtime_error("Failed to read size from file stream.");
	}

	const size_t word_count = (bs.size() + ALPHA - 1) / ALPHA;
	bs.m_data.resize(word_count, 0); // Resize before reading into it

	if (word_count > 0) { // Only read if there's data expected
		file.read(reinterpret_cast<char*>(bs.m_data.data()), word_count * sizeof(uintmax_t));
		if (!file) { // Check read success
			throw std::runtime_error("Failed to read data from file stream.");
		}
	}


	return bs;
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
BitString<CHAR_T, CHAR_SIZE_BITS> BitString<CHAR_T, CHAR_SIZE_BITS>::from_file(const std::string& filename)
{
	std::ifstream file(filename, std::ios::binary);
	if (!file)
	{
		throw std::runtime_error("Failed to open file for reading: " + filename);
	}
	return from_file(file);
	// file automatically closed by destructor
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
std::string BitString<CHAR_T, CHAR_SIZE_BITS>::to_string() const noexcept
{
	std::string str;
	str.reserve(size());

	for (CHAR_T c : *this)
	{
		// Assumes CHAR_T is convertible to char, potential truncation if wider.
		str.push_back(static_cast<char>(c));
	}

	return str;
}
