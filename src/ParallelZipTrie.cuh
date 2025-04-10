/**
 * @file ParallelZipTrie.cuh
 * @brief Defines the ParallelZipTrie class, extending ZipTrie with parallel key comparison capabilities using CUDA.
 * @details This class inherits from `ZipTrie` and overrides key methods (`insert`, `search`, internal `k_compare`)
 * to leverage parallel comparison algorithms (specifically `BitString::par_k_compare`) executed on the GPU via CUDA.
 * It manages the necessary device memory buffers for these parallel operations. The goal is to accelerate
 * trie navigation and updates for scenarios involving long keys or high throughput requirements where
 * parallel comparison offers a performance benefit over sequential CPU comparison.
 * @see ZipTrie
 * @see BitString::par_k_compare
 */
#pragma once

#include "BitString.cuh" // Requires BitString definition with par_k_compare

#include <compare>   // For std::strong_ordering

#include "ZipTrie.hpp" // Base class definition

#include "cuda_utils.cuh" // CUDA utility functions like alloc_to_device, device_free

/**
 * @class ParallelZipTrie
 * @brief Extends the `ZipTrie` class to incorporate parallel key comparisons leveraging GPU acceleration (CUDA).
 *
 * @details This class inherits the core structure and logic of `ZipTrie` but overrides key comparison-related methods
 * (`k_compare`, `search`, `insert`) to utilize a parallel comparison function (`par_k_compare` from `BitString`).
 * It manages CUDA device memory buffers (`d_a`, `d_largeblock`) required for these parallel comparisons,
 * allocating them in the constructor and freeing them in the destructor. This allows for potentially faster
 * trie operations when comparing long keys, at the cost of requiring CUDA-enabled hardware and managing GPU memory.
 * The LCP-based prefix check optimization (`k_compare_prefix_check`) from the base class is still utilized
 * to avoid unnecessary parallel comparisons when possible.
 *
 * @tparam CHAR_T The underlying character type used in the keys (`BitString`).
 * @tparam MEMORY_EFFICIENT If true, uses `MemoryEfficientLCP`; otherwise, uses `unsigned` for LCP storage.
 * @tparam RANK_T The type used for node ranks (default: `GeometricRank`).
 * @tparam CHAR_SIZE_BITS The number of significant bits per character in `CHAR_T`.
 * @see ZipTrie
 * @see BitString::par_k_compare
 */
template <typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T = GeometricRank, unsigned CHAR_SIZE_BITS = sizeof(CHAR_T) * 8>
class ParallelZipTrie : public ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>
{
public:
	/** @brief Alias for the base ZipTrie class template instantiation for brevity. */
	using ZT = ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>;

	/** @brief Alias for the LCP type used, inherited from the base class (`ZipTrie::LCP_T`). */
	using LCP_T = typename ZT::LCP_T;
	/** @brief Alias for the Key type used, inherited from the base class (`ZipTrie::KEY_T`). */
	using KEY_T = typename ZT::KEY_T;

	/**
	 * @brief Constructs a ParallelZipTrie.
	 * @details Initializes the base `ZipTrie` and allocates necessary CUDA device memory buffers
	 * (`d_a`, `d_largeblock`) required for the parallel key comparisons (`BitString::par_k_compare`).
	 * The size of the allocated buffers is determined by the `max_lcp_length` parameter.
	 * @param max_size Hint for reserving memory for the expected maximum number of nodes (passed to base class constructor).
	 * @param max_lcp_length The maximum possible LCP length (in characters/bits), used for device memory allocation size calculation.
	 * @see alloc_to_device
	 * @see alloc_large_block_to_device_s
	 */
	ParallelZipTrie(unsigned max_size, unsigned max_lcp_length);

	/**
	 * @brief Destructor for ParallelZipTrie.
	 * @details Frees the CUDA device memory buffers (`d_a`, `d_largeblock`) allocated in the constructor.
	 * @see device_free
	 */
	~ParallelZipTrie();

	/**
	 * @brief Returns the number of nodes (keys) currently stored in the trie.
	 * @details Inherited directly from the base `ZipTrie` class via the `using` declaration.
	 * @return unsigned The number of nodes.
	 */
	using ZT::size;

	/**
	 * @brief Inserts a key into the parallel zip trie using parallel comparison logic.
	 *
	 * @details This method overrides the base class `insert`. It adds the key with a random rank
	 * to the internal storage (`_buckets`) and then calls the recursive insertion helper (`insert_recursive`),
	 * providing a comparison function lambda that utilizes the parallel `k_compare` implementation of this class.
	 * The parallel `k_compare` leverages `BitString::par_k_compare` for GPU acceleration when beneficial.
	 *
	 * @param key Pointer to the new key to insert. Must be non-null and point to a valid `KEY_T` object.
	 * @warning **Key Lifetime:** The lifetime of the object pointed to by `key` must exceed the lifetime of the ParallelZipTrie.
	 * @warning **Uniqueness:** Assumes keys are unique. Inserting duplicate keys results in the insertion being skipped,
	 * but the initially added bucket might remain unused in storage.
	 * @see insert_recursive
	 * @see k_compare(const KEY_T*, const Bucket&, const AncestorLCPs&, size_t&, size_t&) const
	 * @see RANK_T::get_random
	 */
	void insert(const KEY_T* key) noexcept;

	/**
	 * @brief Removes a node with a given key from the zip tree. (Not implemented)
	 * @param key key of node to remove
	 * @return true if a node was removed, false otherwise
	 * @note This functionality is declared but **not implemented**. Requires implementing 'unzip' operations.
	 */
	// bool remove(const KEY_T& key) noexcept;

	/** @brief Alias for the `SearchResults` struct defined in the base `ZipTrie` class. */
	using SearchResults = typename ZT::SearchResults;

	/**
	 * @brief Performs a search for a key within the ParallelZipTrie using parallel comparison logic.
	 * @details Overrides the base class `search` method. It initializes state required for parallel comparison
	 * (`comparison_size`, `max_copied`) and calls the recursive search helper (`search_recursive`) inherited
	 * from the base class. Crucially, it passes a comparison function lambda that invokes the parallel
	 * `k_compare` implementation defined within this `ParallelZipTrie` class.
	 * @param key Pointer to the key to search for.
	 * @return SearchResults A struct containing the search outcome (found status, LCP, depth).
	 * @see search_recursive
	 * @see k_compare(const KEY_T*, const Bucket&, const AncestorLCPs&, size_t&, size_t&) const
	 */
	SearchResults search(const KEY_T* key) const noexcept override; // Mark as override

protected:
	// --- Inherited members from ZipTrie ---
	/** @brief Using declaration for the root index from the base class `ZipTrie`. */
	using ZT::_root_index;
	/** @brief Using declaration for the maximum LCP length storage from the base class `ZipTrie`. */
	using ZT::_max_lcp_length;
	/** @brief Using declaration for the log max size storage from the base class `ZipTrie`. */
	using ZT::_log_max_size;
	/** @brief Using declaration for the node storage vector (`std::vector<Bucket>`) from the base class `ZipTrie`. */
	using ZT::_buckets;

	/** @brief Using declaration for the null pointer constant (`unsigned`) from the base class `ZipTrie`. */
	using ZT::NULLPTR;

	// --- Inherited types/structs from ZipTrie ---
	/** @brief Alias for the `AncestorLCPs` struct from the base class `ZipTrie`. */
	using AncestorLCPs = typename ZT::AncestorLCPs;
	/** @brief Alias for the `Bucket` (node) struct from the base class `ZipTrie`. */
	using Bucket = typename ZT::Bucket;
	/** @brief Alias for the `ComparisonResult` struct from the base class `ZipTrie`. */
	using ComparisonResult = typename ZT::ComparisonResult;

	// --- Inherited methods from ZipTrie ---
	/** @brief Using declaration for the LCP approximation helper (if `MEMORY_EFFICIENT` is true) from `ZipTrie`. */
	using ZT::get_memory_efficient_lcp;
	/** @brief Using declaration for the recursive search helper from the base class `ZipTrie`. */
	using ZT::search_recursive;
	/** @brief Using declaration for the recursive insert helper from the base class `ZipTrie`. */
	using ZT::insert_recursive;
	/** @brief Using declaration for the LCP prefix check optimization helper from `ZipTrie`. */
	using ZT::k_compare_prefix_check;
	/** @brief Using declaration for the LCP value conversion helper from `ZipTrie`. */
	using ZT::convert_lcp;

private:
	/** @brief Pointer to device memory buffer 'a', used as primary input/output for parallel key comparisons (e.g., via `BitString::par_k_compare`). */
	uintmax_t* d_a;
	/** @brief Pointer to auxiliary device memory buffer ('largeblock'), used as temporary storage or workspace during parallel key comparisons by kernels like `par_find_mismatch_s`. */
	uintmax_t* d_largeblock;

	/**
	 * @brief Performs key comparison using parallel logic, overriding the base class `k_compare`.
	 *
	 * @details This method implements the core parallel comparison strategy for `ParallelZipTrie`.
	 * It first leverages the inherited `k_compare_prefix_check` optimization. If a full comparison
	 * is necessary, it iteratively calls the `BitString::par_k_compare` function, which performs
	 * the comparison on the GPU. The `comparison_size` parameter controls the chunk size (number of characters)
	 * compared in each parallel step and is adapted dynamically: it increases if blocks match to compare
	 * larger chunks next time, and decreases if a mismatch is found to potentially refine the LCP more quickly
	 * in subsequent comparisons within the same search/insert operation.
	 *
	 * @param x Pointer to the key being searched for or inserted.
	 * @param v Constant reference to the current node (`Bucket`) being compared against.
	 * @param ancestor_lcps The LCPs computed along the path from the root down to `v`'s parent.
	 * @param[in, out] comparison_size Controls the block size (number of characters) for the parallel comparison step.
	 * Its value is adapted based on the comparison outcome.
	 * @param[in, out] max_copied Tracks state related to data copied to device memory (`d_a`) by `par_k_compare`
	 * to potentially optimize subsequent copies within the same search/insert operation.
	 * @return ComparisonResult A struct containing the relative ordering (`comparison`) of `x` vs `v`
	 * and the calculated Longest Common Prefix (`lcp`) between them (converted to `LCP_T`).
	 * @note This method directly interacts with the `BitString::par_k_compare` function and relies on the device
	 * memory buffers `d_a` and `d_largeblock` managed by this class.
	 * @see ZT::k_compare_prefix_check
	 * @see BitString::par_k_compare
	 * @see convert_lcp
	 */
	inline ComparisonResult k_compare(const KEY_T* x, const Bucket& v, const AncestorLCPs& ancestor_lcps, size_t& comparison_size, size_t& max_copied) const noexcept;
};

//-----------------------------------------------------------------------------
// ParallelZipTrie Method Implementations
//-----------------------------------------------------------------------------

template<typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
ParallelZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::ParallelZipTrie(unsigned max_size, unsigned max_lcp_length):
	ZT(max_size, max_lcp_length) // Call base class constructor
{
	// Calculate required device buffer size in words based on max LCP length.
	size_t max_lcp_length_words = (max_lcp_length + KEY_T::ALPHA - 1) / KEY_T::ALPHA;

	// Allocate the primary device buffer.
	d_a = alloc_to_device<uintmax_t>(max_lcp_length_words);
	// Allocate the auxiliary large block buffer (size calculation might depend on the parallel algorithm used).
	// Assuming alloc_large_block_to_device_s is defined appropriately, potentially in cuda_utils.cuh.
	d_largeblock = alloc_large_block_to_device_s(max_lcp_length_words);
}

template<typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
ParallelZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::~ParallelZipTrie()
{
	// Free the allocated device memory buffers.
	device_free(d_a);
	device_free(d_largeblock);
}

template<typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
typename ParallelZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::ComparisonResult ParallelZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::k_compare(const KEY_T* x, const Bucket& v, const AncestorLCPs& ancestor_lcps, size_t& comparison_size, size_t& max_copied) const noexcept
{
	size_t start_lcp_val = 0; // Initialize LCP offset for comparison.

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
		comparison_size = std::max(KEY_T::MIN_PAR_COMPARE_CHAR_SIZE, comparison_size / 2);
		// Break the loop as the final comparison result is determined.
		break;
	}

	// Convert the final exact LCP length to the appropriate LCP_T format and return the result.
	return { final_comparison, convert_lcp(final_lcp_val) };
}

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

	// Call the inherited recursive insert function to place the node and perform zip operations.
	// Pass the new node details, initial empty ancestor LCPs, and the parallel comparison lambda.
	// Update the root index, as it might change due to zip operations.
	_root_index = insert_recursive(&_buckets[new_node_index], new_node_index, _root_index, {}, parallel_compare);
}
