/**
 * @file ParallelSkipTrie.cuh
 * @brief Defines the ParallelSkipTrie class, a GPU-accelerated variant of SkipTrie for efficient string operations.
 * @details This file contains the definition of the ParallelSkipTrie template class, which extends the SkipTrie
 * data structure with parallel GPU-accelerated comparison operations. It leverages CUDA for fast string
 * comparisons, particularly beneficial for long strings or large datasets. The class maintains the skip list
 * structure of its parent class while enhancing performance through parallel processing capabilities.
 * @see SkipTrie
 * @see BitString
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
 * @brief Implements a GPU-accelerated skip list data structure optimized for efficient string operations.
 * @details This class extends SkipTrie by implementing parallel comparison operations using CUDA.
 * The primary performance enhancement comes from parallel key comparison (`par_k_compare`), which
 * offloads string comparison work to the GPU. This is particularly effective for long strings
 * or large datasets where sequential comparisons would otherwise become a bottleneck.
 * The class maintains GPU memory buffers for performing these parallel operations efficiently.
 *
 * @tparam CHAR_T The underlying character type used in the keys (`BitString`).
 * @tparam CHAR_SIZE_BITS The number of significant bits per character in `CHAR_T`. Defaults to `sizeof(CHAR_T) * 8`.
 * @see SkipTrie
 * @see BitString
 */
template<typename CHAR_T, unsigned CHAR_SIZE_BITS = sizeof(CHAR_T) * 8>
class ParallelSkipTrie : public SkipTrie<CHAR_T, CHAR_SIZE_BITS>
{
public:
	using ST = SkipTrie<CHAR_T, CHAR_SIZE_BITS>;

	using KEY_T = BitString<CHAR_T, CHAR_SIZE_BITS>;

	/**
	 * @brief Constructs a new ParallelSkipTrie with allocated GPU memory for parallel operations.
	 * @details Initializes the skip list structure and allocates GPU memory buffers for parallel string
	 * comparisons. The size of the allocated memory is determined by the `max_size` parameter.
	 * @param max_size The maximum total size (in characters) of keys that will be stored in the trie,
	 * used to allocate sufficient GPU memory for parallel comparison operations.
	 */
	ParallelSkipTrie(size_t max_size);

	/**
	 * @brief Destructor for the ParallelSkipTrie class.
	 * @details Frees all GPU memory allocated for parallel comparison operations,
	 * including the device arrays used for string comparison.
	 */
	~ParallelSkipTrie();

	using ST::get_random_height;

	/**
	 * @brief Inserts a new key into the skip list with a randomly generated height.
	 * @details Calls `get_random_height()` to determine the height and then calls the height-specific insert method.
	 * Initializes the GPU comparison parameters before delegating to the parent class implementation.
	 * @param key A pointer to the `BitString` key to insert. The ParallelSkipTrie stores this pointer directly.
	 * @return bool True if the key was inserted successfully, false if the key already exists.
	 * @warning Key Lifetime: The lifetime of the object pointed to by `key` must exceed the lifetime of the ParallelSkipTrie.
	 * @see insert(const KEY_T*, size_t)
	 * @see get_random_height()
	 */
	bool insert(const KEY_T* key) noexcept;

	/**
	 * @brief Inserts a new key into the skip list with a specified height.
	 * @details Initializes the GPU comparison parameters (m_comparison_size, m_max_copied) before 
	 * calling the parent class's implementation. Uses parallel comparison operations on the GPU
	 * for efficient string comparisons during the insertion process.
	 * @param key A pointer to the `BitString` key to insert. The ParallelSkipTrie stores this pointer directly.
	 * @param height The height (number of levels above the base) for the new node.
	 * @return bool True if the key was inserted successfully, false if the key already exists.
	 * @warning Key Lifetime: The lifetime of the object pointed to by `key` must exceed the lifetime of the ParallelSkipTrie.
	 * @see SkipTrie::insert
	 */
	bool insert(const KEY_T* key, size_t height) noexcept;

protected:
	/** @brief Maximum total size of keys that can be stored, used for GPU memory allocation. */
	size_t m_max_size;
	
	/** @brief Device (GPU) memory for parallel string comparison operations. */
	uintmax_t* d_a;
	
	/** @brief Large block of device (GPU) memory for more complex parallel operations. */
	uintmax_t* d_largeblock;
	
	/** @brief Current size (in characters) used for parallel comparisons, adjusted dynamically. */
	mutable size_t m_comparison_size;
	
	/** @brief Tracks the maximum number of characters copied to the GPU during comparison operations. */
	mutable size_t m_max_copied;

	/** @brief Alias for the parent class's EqualOrSuccessor struct, used for search results. */
	using EqualOrSuccessor = ST::EqualOrSuccessor;

	/**
	 * @brief Finds the node containing the exact key or its immediate successor, using parallel comparison.
	 * @details Overrides the parent class's method to initialize GPU comparison parameters before
	 * performing the search. Resets the comparison state (m_comparison_size and m_max_copied) 
	 * before delegating to the parent class implementation.
	 * @param key A pointer to the `BitString` key to search for.
	 * @param require_level0 If true, ensures the search descends fully to level 0 before returning the node.
	 * @return EqualOrSuccessor A struct containing the found node (match or successor), an equality flag, and the search path LCP.
	 * @note This method resets the LCP calculation state as part of its operation.
	 * @see SkipTrie::find_equal_or_successor
	 */
	EqualOrSuccessor find_equal_or_successor(const KEY_T* key, const bool require_level0) const noexcept;

	/** @brief Alias for the parent class's ResultLCP struct, used for comparison results. */
	using ResultLCP = ST::ResultLCP;

	/**
	 * @brief Compares two keys using GPU-accelerated parallel comparison.
	 * @details Overrides the parent class's method to use GPU-accelerated string comparison.
	 * Uses `BitString::par_k_compare` to perform the comparison on the GPU, significantly improving
	 * performance for long strings. Adaptively adjusts the comparison size (m_comparison_size)
	 * based on comparison results, doubling it when more characters need to be compared and
	 * halving it (but not below minimum) when the comparison is complete.
	 * 
	 * @param key1 Pointer to the first `BitString` key.
	 * @param key2 Pointer to the second `BitString` key.
	 * @param lcp The length of the known common prefix (in characters) from which to start the comparison.
	 * @return ResultLCP A struct containing the comparison result (`std::strong_ordering`) and the total LCP length found.
	 * @see SkipTrie::compare
	 * @see BitString::par_k_compare
	 */
	ResultLCP compare(const KEY_T* key1, const KEY_T* key2, size_t lcp) const noexcept;
};

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
ParallelSkipTrie<CHAR_T, CHAR_SIZE_BITS>::ParallelSkipTrie(size_t max_size)
	: ST(), m_max_size(max_size)
{
	// Calculate number of words needed to store max_size characters
	// Use ceiling division (max_size / ALPHA + potential remainder)
	size_t max_size_words = (max_size + KEY_T::ALPHA - 1) / KEY_T::ALPHA;
	
	// Allocate device memory for string comparison operations
	d_a = alloc_to_device<uintmax_t>(max_size_words);
	
	// Allocate larger block for more complex parallel operations
	d_largeblock = alloc_large_block_to_device_s(max_size_words);
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
ParallelSkipTrie<CHAR_T, CHAR_SIZE_BITS>::~ParallelSkipTrie()
{
	// Free GPU memory allocated for string comparison
	device_free(d_a);
	
	// Free GPU memory allocated for large block operations
	device_free(d_largeblock);
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
auto ParallelSkipTrie<CHAR_T, CHAR_SIZE_BITS>::compare(const KEY_T* key1, const KEY_T* key2, size_t lcp) const noexcept -> ResultLCP
{
	// Initialize result with equality and starting LCP value
	ResultLCP result = { std::strong_ordering::equal, lcp };

	while (true)
	{
		// Perform parallel comparison on GPU using the current comparison size
		// This is the core GPU-accelerated operation that makes ParallelSkipTrie faster
		auto [comparison, next_lcp] = key1->par_k_compare(*key2, lcp, m_comparison_size, d_a, d_largeblock, m_max_copied);

		// Update the result with comparison outcome and new LCP length
		result.result = comparison;
		result.lcp = next_lcp;

		if (comparison == std::strong_ordering::equal)
		{
			// If both keys match completely, we're done
			if (result.lcp == key1->size() && result.lcp == key2->size())
			{
				return result;
			}

			// Keys match so far but are longer - double comparison size for next iteration
			// This adaptive sizing improves performance by processing more characters at once
			m_comparison_size *= 2;
			lcp = result.lcp;
		}
		else
		{
			// Found a difference between keys, exit the loop
			break;
		}
	}

	// Halve comparison size for future operations, but don't go below minimum
	// This adapts the comparison granularity for optimal performance
	m_comparison_size = std::max(KEY_T::MIN_PAR_COMPARE_CHAR_SIZE, m_comparison_size / 2);

	return result;
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
bool ParallelSkipTrie<CHAR_T, CHAR_SIZE_BITS>::insert(const KEY_T* key) noexcept
{
	// Call the height-specific insert with a random height
	// This follows skip list standard practice of randomized height assignment
	return insert(key, get_random_height());
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
bool ParallelSkipTrie<CHAR_T, CHAR_SIZE_BITS>::insert(const KEY_T* key, size_t height) noexcept
{
	// Reset the GPU comparison parameters to their initial values
	// This ensures optimal starting point for parallel comparison operations
	m_comparison_size = KEY_T::MIN_PAR_COMPARE_CHAR_SIZE;
	m_max_copied = 0;

	// Delegate to the parent class implementation after initializing GPU parameters
	return ST::insert(key, height);
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
typename ParallelSkipTrie<CHAR_T, CHAR_SIZE_BITS>::EqualOrSuccessor ParallelSkipTrie<CHAR_T, CHAR_SIZE_BITS>::find_equal_or_successor(const KEY_T* key, const bool require_level0) const noexcept
{
	// Reset GPU comparison parameters to their initial values
	// This ensures consistent behavior across multiple search operations
	m_comparison_size = KEY_T::MIN_PAR_COMPARE_CHAR_SIZE;
	m_max_copied = 0;

	// Delegate to parent class implementation after initializing GPU parameters
	// The actual comparisons will use the overridden compare() method which uses GPU acceleration
	return ST::find_equal_or_successor(key, require_level0);
}
