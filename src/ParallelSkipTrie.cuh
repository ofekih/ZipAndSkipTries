/**
 * @file ParallelSkipTrie.hpp
 * @brief Definition of the ParallelSkipTrie class for efficient string storage and retrieval.
 *
 * The ParallelSkipTrie is a templated class designed to store strings in a manner that allows for fast insertion, deletion,
 * and lookup operations. It is particularly optimized for strings, utilizing a bit-level representation to minimize
 * comparison operations. The implementation uses a skip list data structure, where each node has links that skip over
 * multiple elements, allowing for logarithmic average time complexities for its primary operations.
 */

#pragma once

#include <cmath> // for log2
#include <cstddef> // for size_t
#include <memory> // for smart pointers
#include <random> // for random height generation
#include <vector> // for dynamic arrays

#include <unordered_map>

#include <fstream> // for file operations

#include "BitString.cuh" // for string representation
#include "SkipTrie.hpp" // for Direction

#include "cuda_utils.cuh"
#include "msw.cuh"


/**
 * @class ParallelSkipTrie
 * @brief Implements a skip list specifically optimized for string storage.
 *
 * @tparam CHAR_T The character type of the strings to be stored. Defaults to char.
 * @tparam CHAR_SIZE_BITS The size in bits of the character type. Defaults to the bit-size of CHAR_T.
 */
template<typename CHAR_T, unsigned CHAR_SIZE_BITS = sizeof(CHAR_T) * 8>
class ParallelSkipTrie : public SkipTrie<CHAR_T, CHAR_SIZE_BITS>
{
public:
	using ST = SkipTrie<CHAR_T, CHAR_SIZE_BITS>;

	using KEY_T = BitString<CHAR_T, CHAR_SIZE_BITS>;

	/**
	 * @brief Constructs a new String Skip List object.
	 */
	ParallelSkipTrie(size_t max_size);

	~ParallelSkipTrie();

	using ST::get_random_height;

	bool insert(const KEY_T* key) noexcept;

	/**
	 * @brief Inserts a new key into the skip list with a specified height.
	 * @param key A pointer to the BitString representing the key to insert.
	 * @param height The height at which to insert the key.
	 */
	bool insert(const KEY_T* key, size_t height) noexcept;

protected:
	size_t m_max_size;
	uintmax_t* d_a;
	uintmax_t* d_largeblock;
	mutable size_t m_comparison_size;
	mutable size_t m_max_copied;

	using EqualOrSuccessor = ST::EqualOrSuccessor;

	// NOTE: RESETS LCP
	EqualOrSuccessor find_equal_or_successor(const KEY_T* key, const bool require_level0) const noexcept;

	using ResultLCP = ST::ResultLCP;

	ResultLCP compare(const KEY_T* key1, const KEY_T* key2, size_t lcp) const noexcept;
};

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
ParallelSkipTrie<CHAR_T, CHAR_SIZE_BITS>::ParallelSkipTrie(size_t max_size)
	: ST(), m_max_size(max_size)
{
	size_t max_size_words = (max_size + KEY_T::ALPHA - 1) / KEY_T::ALPHA;
	d_a = alloc_to_device<uintmax_t>(max_size_words);
	d_largeblock = alloc_large_block_to_device_s(max_size_words);
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
ParallelSkipTrie<CHAR_T, CHAR_SIZE_BITS>::~ParallelSkipTrie()
{
	device_free(d_a);
	device_free(d_largeblock);
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
auto ParallelSkipTrie<CHAR_T, CHAR_SIZE_BITS>::compare(const KEY_T* key1, const KEY_T* key2, size_t lcp) const noexcept -> ResultLCP
{
	ResultLCP result = { std::strong_ordering::equal, lcp };

	while (true)
	{
		auto [comparison, next_lcp] = key1->par_k_compare(*key2, lcp, m_comparison_size, d_a, d_largeblock, m_max_copied);

		result.result = comparison;
		result.lcp = next_lcp;

		if (comparison == std::strong_ordering::equal)
		{
			if (result.lcp == key1->size() && result.lcp == key2->size())
			{
				return result;
			}

			m_comparison_size *= 2;
			lcp = result.lcp;
		}
		else
		{
			break;
		}
	}

	m_comparison_size = std::max(KEY_T::MIN_PAR_COMPARE_CHAR_SIZE, m_comparison_size / 2);

	return result;
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
bool ParallelSkipTrie<CHAR_T, CHAR_SIZE_BITS>::insert(const KEY_T* key) noexcept
{
	return insert(key, get_random_height());
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
bool ParallelSkipTrie<CHAR_T, CHAR_SIZE_BITS>::insert(const KEY_T* key, size_t height) noexcept
{
	m_comparison_size = KEY_T::MIN_PAR_COMPARE_CHAR_SIZE;
	m_max_copied = 0;

	return ST::insert(key, height);
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
typename ParallelSkipTrie<CHAR_T, CHAR_SIZE_BITS>::EqualOrSuccessor ParallelSkipTrie<CHAR_T, CHAR_SIZE_BITS>::find_equal_or_successor(const KEY_T* key, const bool require_level0) const noexcept
{
	m_comparison_size = KEY_T::MIN_PAR_COMPARE_CHAR_SIZE;
	m_max_copied = 0;

	return ST::find_equal_or_successor(key, require_level0);
}
