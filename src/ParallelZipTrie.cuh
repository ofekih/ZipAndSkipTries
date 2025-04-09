/**
 * @file ParallelZipTrie.cuh
 * @brief Defines the ParallelZipTrie class, extending ZipTrie with parallel key comparison capabilities using CUDA.
 */
#pragma once

#include "BitString.cuh"

#include <compare>   // For std::strong_ordering

#include "ZipTrie.hpp"

#include "utility.cuh" // CUDA utility functions like alloc_to_device, device_free

/**
 * @class ParallelZipTrie
 * @brief Extends the ZipTrie class to incorporate parallel key comparisons leveraging GPU acceleration (CUDA).
 *
 * This class inherits the core structure and logic of `ZipTrie` but overrides key comparison-related methods
 * (`k_compare`, `search`, `insert`) to utilize a parallel comparison function (`par_k_compare` from `BitString`).
 * It manages CUDA device memory buffers (`d_a`, `d_largeblock`) required for these parallel comparisons.
 *
 * @tparam CHAR_T The underlying character type used in the keys (`BitString`).
 * @tparam MEMORY_EFFICIENT If true, uses `MemoryEfficientLCP`; otherwise, uses `unsigned`.
 * @tparam RANK_T The type used for node ranks (default: `GeometricRank`).
 * @tparam CHAR_SIZE_BITS The number of significant bits per character in `CHAR_T`.
 */
template <typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T = GeometricRank, unsigned CHAR_SIZE_BITS = sizeof(CHAR_T) * 8>
class ParallelZipTrie : public ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>
{
public:
	/** @brief Alias for the base ZipTrie class template instantiation. */
	using ZT = ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>;

	/** @brief Alias for the LCP type used, inherited from the base class. */
	using LCP_T = typename ZT::LCP_T;
	/** @brief Alias for the Key type used, inherited from the base class. */
	using KEY_T = typename ZT::KEY_T;

	/**
	 * @brief Constructs a ParallelZipTrie.
	 * Initializes the base ZipTrie and allocates necessary CUDA device memory buffers
	 * for parallel key comparisons based on the maximum LCP length.
	 * @param max_size Hint for reserving memory for the expected maximum number of nodes (passed to base class).
	 * @param max_lcp_length The maximum possible LCP length, used for device memory allocation size calculation.
	 */
	ParallelZipTrie(unsigned max_size, unsigned max_lcp_length);

	/**
	 * @brief Destructor for ParallelZipTrie.
	 * Frees the allocated CUDA device memory buffers (`d_a`, `d_largeblock`).
	 */
	~ParallelZipTrie();

	/** @brief Inherits the `size()` method directly from the base `ZipTrie` class. */
	using ZT::size;

	/**
	 * @brief Inserts a key into the parallel zip trie using parallel comparison logic.
	 *
	 * This method overrides the base class `insert`. It adds the key with a random rank
	 * to the internal storage and then calls the recursive insertion helper (`insert_recursive`),
	 * providing a comparison function that utilizes the parallel `k_compare` implementation.
	 *
	 * @warning Assumes keys are unique. Inserting duplicate keys leads to undefined behavior
	 * (specifically, the duplicate key might remain in the `_buckets` vector but not be part of the tree).
	 * @warning The lifetime of the object pointed to by `key` must exceed the lifetime of the ParallelZipTrie.
	 *
	 * @param key Pointer to the new key to insert. Must be non-null and point to a valid `KEY_T` object.
	 */
	void insert(const KEY_T* key) noexcept;

	/**
	 * @brief Removes a node with a given key from the zip tree. (Not implemented)
	 * @param key key of node to remove
	 * @return true if a node was removed, false otherwise
	 * @note This functionality is declared but not implemented.
	 */
	// bool remove(const KEY_T& key) noexcept;

	/** @brief Alias for the SearchResults struct defined in the base class. */
	using SearchResults = typename ZT::SearchResults;

	/**
	 * @brief Performs a search for a key within the ParallelZipTrie using parallel comparison logic.
	 * Overrides the base class `search` method. It initializes state for parallel comparison
	 * and calls the recursive search helper (`search_recursive`) with a comparison function
	 * that uses the parallel `k_compare` implementation.
	 * @param key Pointer to the key to search for.
	 * @return SearchResults A struct containing the search outcome (found status, LCP, depth).
	 */
	SearchResults search(const KEY_T* key) const noexcept;

protected:
	// --- Inherited members from ZipTrie ---
	/** @brief Using declaration for the root index pointer from the base class. */
	using ZT::_root_index;
	/** @brief Using declaration for the maximum LCP length storage from the base class. */
	using ZT::_max_lcp_length;
	/** @brief Using declaration for the log max size storage from the base class. */
	using ZT::_log_max_size;
	/** @brief Using declaration for the node storage vector from the base class. */
	using ZT::_buckets;

	/** @brief Using declaration for the null pointer constant from the base class. */
	using ZT::NULLPTR;

	// --- Inherited types/structs from ZipTrie ---
	/** @brief Alias for the AncestorLCPs struct from the base class. */
	using AncestorLCPs = typename ZT::AncestorLCPs;
	/** @brief Alias for the Bucket (node) struct from the base class. */
	using Bucket = typename ZT::Bucket;
	/** @brief Alias for the ComparisonResult struct from the base class. */
	using ComparisonResult = typename ZT::ComparisonResult;

	// --- Inherited methods from ZipTrie ---
	/** @brief Using declaration for the LCP conversion helper from the base class. */
	using ZT::get_memory_efficient_lcp; // Note: This might need adjustment if parallel LCP calculation differs
	/** @brief Using declaration for the recursive search helper from the base class. */
	using ZT::search_recursive;
	/** @brief Using declaration for the recursive insert helper from the base class. */
	using ZT::insert_recursive;
	/** @brief Using declaration for the LCP prefix check helper from the base class. */
	using ZT::k_compare_prefix_check;
	/** @brief Using declaration for the LCP value conversion helper from the base class. */
	using ZT::convert_lcp;

private:
	/** @brief Pointer to device memory buffer 'a', used for parallel key comparisons (e.g., via `par_k_compare`). */
	uintmax_t* d_a;
	/** @brief Pointer to another device memory buffer ('largeblock'), used as temporary storage or workspace during parallel key comparisons. */
	uintmax_t* d_largeblock;

	/**
	 * @brief Performs key comparison using parallel logic, overriding the base class `k_compare`.
	 *
	 * This method implements the core parallel comparison strategy. It first uses the base class's
	 * `k_compare_prefix_check` to leverage ancestor LCP optimizations. If a full comparison is needed,
	 * it iteratively calls the parallel comparison function `par_k_compare`.
	 * The `comparison_size` parameter controls the chunk size for parallel comparison and is
	 * adapted dynamically based on whether a mismatch is found or the end of the keys is reached.
	 *
	 * @param x Pointer to the key being searched for or inserted.
	 * @param v Constant reference to the current node (`Bucket`) being compared against.
	 * @param ancestor_lcps The LCPs computed along the path from the root down to `v`'s parent.
	 * @param[in, out] comparison_size Controls the block size or granularity of the parallel comparison. May be adjusted by this function.
	 * @param[in, out] max_copied Tracks state related to data copied to/from device memory, potentially for optimization.
	 * @return ComparisonResult A struct containing the relative ordering (`comparison`) of `x` vs `v`
	 * and the calculated Longest Common Prefix (`lcp`) between them.
	 * @note This method directly interacts with the `par_k_compare` function and relies on the device
	 * memory buffers `d_a` and `d_largeblock`.
	 */
	inline ComparisonResult k_compare(const KEY_T* x, const Bucket& v, const AncestorLCPs& ancestor_lcps, size_t& comparison_size, size_t& max_copied) const noexcept;
};

/**
 * @brief Constructor implementation for ParallelZipTrie.
 * Initializes the base class and allocates CUDA device memory.
 */
template<typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
ParallelZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::ParallelZipTrie(unsigned max_size, unsigned max_lcp_length):
	ZT(max_size, max_lcp_length) // Call base class constructor
{
	// Calculate required device buffer size in words based on max LCP length.
	size_t max_lcp_length_words = (max_lcp_length + KEY_T::ALPHA - 1) / KEY_T::ALPHA;

	d_a = alloc_to_device<uintmax_t>(max_lcp_length_words);
	d_largeblock = alloc_large_block_to_device_s(max_lcp_length_words);
}

/**
 * @brief Destructor implementation for ParallelZipTrie.
 * Frees allocated CUDA device memory.
 */
template<typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
ParallelZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::~ParallelZipTrie()
{
	device_free(d_a);
	device_free(d_largeblock);
}

/**
 * @brief Implementation of the parallel key comparison logic.
 * Uses prefix check optimization and iterative calls to `par_k_compare`.
 */
template<typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
typename ParallelZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::ComparisonResult ParallelZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::k_compare(const KEY_T* x, const Bucket& v, const AncestorLCPs& ancestor_lcps, size_t& comparison_size, size_t& max_copied) const noexcept
{
	size_t start_lcp_val = 0; // Initialize LCP offset for comparison.

	// Use the base class's prefix check for initial optimization.
	if (auto prefix_result = k_compare_prefix_check(v, ancestor_lcps, start_lcp_val))
	{
		return *prefix_result; // Return early if comparison decided by LCPs.
	}

	// Initialize comparison result and final LCP with the value from prefix check.
	std::strong_ordering final_comparison = std::strong_ordering::equal;
	size_t final_lcp_val = start_lcp_val;

	// Iteratively call parallel compare until a mismatch is found or keys are fully compared.
	while (true)
	{
		// Call the parallel comparison function provided by the KEY_T (BitString) class.
		// Pass the keys, current LCP offset, comparison block size, device buffers, and copy tracker.
		auto [current_comparison, next_lcp] = x->par_k_compare(
			*v.key, final_lcp_val, comparison_size,
			d_a, d_largeblock, max_copied);

		// Update the overall comparison result and the total LCP found so far.
		final_comparison = current_comparison;
		final_lcp_val = next_lcp;

		// If the current block compared equal...
		if (final_comparison == std::strong_ordering::equal)
		{
			// Check if we have compared the entire length of both keys.
			// If keys are fully compared and equal, break the loop.
			if (final_lcp_val == x->size() && final_lcp_val == v.key->size())
			{
				break;
			}
			// Otherwise, keys matched so far but aren't finished.
			// Increase the comparison size (e.g., double it) for the next iteration to potentially speed up.
			comparison_size *= 2;
			// Continue to the next iteration to compare the next block.
			continue;
		}

		// Mismatch found (comparison is less or greater).
		// Reduce the comparison size for future comparisons (adaptive strategy).
		// Ensure it doesn't go below a minimum threshold.
		comparison_size = std::max(BitString<CHAR_T, CHAR_SIZE_BITS>::MIN_PAR_COMPARE_CHAR_SIZE, comparison_size / 2);
		// Break the loop as the final comparison result is determined.
		break;
	}

	// Convert the final exact LCP length to the appropriate LCP_T format and return the result.
	return { final_comparison, convert_lcp(final_lcp_val) };
}

/**
 * @brief Implementation of the parallel search method.
 * Sets up and calls the recursive search helper with a parallel comparison lambda.
 */
template<typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
typename ParallelZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::SearchResults ParallelZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::search(const KEY_T* key) const noexcept
{
	// Initialize state variables for the parallel comparison process.
	size_t comparison_size = KEY_T::MIN_PAR_COMPARE_CHAR_SIZE;
	size_t max_copied = 0; // Tracks device memory copy status/optimization.

	// Define a lambda function that captures the current state and calls the parallel k_compare method.
	// This lambda will be passed to the base class's search_recursive helper.
	auto parallel_compare = [&](const KEY_T* k, const Bucket& v, const AncestorLCPs& ancestor_lcps)
		{
			// Note: comparison_size and max_copied are captured by reference, allowing k_compare to modify them.
			return k_compare(k, v, ancestor_lcps, comparison_size, max_copied);
		};

	// Call the inherited recursive search function, passing the key and the parallel comparison lambda.
	return search_recursive(key, parallel_compare);
}

/**
 * @brief Implementation of the parallel insert method.
 * Sets up and calls the recursive insert helper with a parallel comparison lambda.
 */
template <typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
void ParallelZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::insert(const KEY_T* key) noexcept
{
	// Initialize state variables for the parallel comparison process.
	size_t comparison_size = KEY_T::MIN_PAR_COMPARE_CHAR_SIZE;
	size_t max_copied = 0;

	// Add a new bucket for the key with a random rank to the storage vector.
	_buckets.push_back({ key, RANK_T::get_random() });
	unsigned new_node_index = _buckets.size() - 1;

	// Define the parallel comparison lambda, capturing state by reference.
	auto parallel_compare = [&](const KEY_T* k, const Bucket& v, const AncestorLCPs& ancestor_lcps)
		{
			return k_compare(k, v, ancestor_lcps, comparison_size, max_copied);
		};

	// Call the inherited recursive insert function to place the node and perform rotations.
	// Pass the new node details, initial empty ancestor LCPs, and the parallel comparison lambda.
	// Update the root index, as it might change due to rotations.
	_root_index = insert_recursive(&_buckets[new_node_index], new_node_index, _root_index, {}, parallel_compare);
}
