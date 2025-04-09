#pragma once

#include "BitString.cuh"

#include <compare>
#include <limits>
#include <vector>

#include <fstream>
#include <string>

#include "ZipTrie.hpp" // for MemoryEfficientLCP and GeometricRank

#include "utility.cuh"

#include <iostream>


template <typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T = GeometricRank, unsigned CHAR_SIZE_BITS = sizeof(CHAR_T) * 8>
class ParallelZipTrie : public ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>
{
public:
	using ZT = ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>;

	using LCP_T = typename ZT::LCP_T;
	using KEY_T = typename ZT::KEY_T;

	ParallelZipTrie(unsigned max_size, unsigned max_lcp_length);

	~ParallelZipTrie();

	using ZT::size;

	/**
	 * Inserts a key, value pair into the zip tree. Note that inserting there is
	 * no validation that the keys don't already exist. Add only unique keys to
	 * avoid undefined behavior.
	 *
	 * @param key new node key
	 * @param val new node value
	 */
	void insert(const KEY_T* key) noexcept;

	/**
	 * Removes a node with a given key from the zip tree.
	 *
	 * @param  key key of node to remove
	 * @return     true if a node was removed, false otherwise
	 */
	// bool remove(const KEY_T& key) noexcept;

	using SearchResults = typename ZT::SearchResults;

	SearchResults search(const KEY_T* key) const noexcept;

protected:
	using ZT::_root_index;
	using ZT::_max_lcp_length;
	using ZT::_log_max_size;
	using ZT::_buckets;
	
	using ZT::NULLPTR;
	
	using AncestorLCPs = typename ZT::AncestorLCPs;
	using Bucket = typename ZT::Bucket;
	using ComparisonResult = typename ZT::ComparisonResult;

	using ZT::get_memory_efficient_lcp;
	using ZT::search_recursive;
	using ZT::insert_recursive;

	using ZT::k_compare_prefix_check;
	using ZT::convert_lcp;

private:
	uintmax_t* d_a;
	uintmax_t* d_largeblock;

	inline ComparisonResult k_compare(const KEY_T* x, const Bucket& v, const AncestorLCPs& ancestor_lcps, size_t& comparison_size, size_t& max_copied) const noexcept;
};

template<typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
ParallelZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::ParallelZipTrie(unsigned max_size, unsigned max_lcp_length): 
	ZT(max_size, max_lcp_length)
{
	size_t max_lcp_length_words = (max_lcp_length + KEY_T::ALPHA - 1) / KEY_T::ALPHA;
	d_a = alloc_to_device<uintmax_t>(max_lcp_length_words);
	d_largeblock = alloc_large_block_to_device_s(max_lcp_length_words);
}

template<typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
ParallelZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::~ParallelZipTrie()
{
	device_free(d_a);
	device_free(d_largeblock);
}

template<typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
typename ParallelZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::ComparisonResult ParallelZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::k_compare(const KEY_T* x, const Bucket& v, const AncestorLCPs& ancestor_lcps, size_t& comparison_size, size_t& max_copied) const noexcept
{
	size_t start_lcp_val = 0;

	if (auto prefix_result = k_compare_prefix_check(v, ancestor_lcps, start_lcp_val))
	{
		return *prefix_result; // Return early result
	}

	std::strong_ordering final_comparison = std::strong_ordering::equal;
	size_t final_lcp_val = start_lcp_val; // Start with LCP from prefix check

	while (true)
	{
		auto [current_comparison, next_lcp] = x->par_k_compare(
			*v.key, final_lcp_val, comparison_size,
			d_a, d_largeblock, max_copied);

		final_comparison = current_comparison;
		final_lcp_val = next_lcp; // Update LCP for potential next iteration

		if (final_comparison == std::strong_ordering::equal)
		{
			// Check completion (ensure KEY_T has size() or similar method)
			if (final_lcp_val == x->size() && final_lcp_val == v.key->size())
			{
				break;
			}
			comparison_size *= 2;
			continue;
		}

		// Mismatch found
		comparison_size = std::max(BitString<CHAR_T, CHAR_SIZE_BITS>::MIN_PAR_COMPARE_CHAR_SIZE, comparison_size / 2);
		break;
	}

	return { final_comparison, convert_lcp(final_lcp_val) };
}

template<typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
typename ParallelZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::SearchResults ParallelZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::search(const KEY_T* key) const noexcept
{
	// Initialize parallel comparison state
	size_t comparison_size = KEY_T::MIN_PAR_COMPARE_CHAR_SIZE;
	size_t max_copied = 0;

	// Define a lambda that calls the parallel k_compare logic
	auto parallel_compare = [&](const KEY_T* k, const Bucket& v, const AncestorLCPs& ancestor_lcps)
		{
			return k_compare(k, v, ancestor_lcps, comparison_size, max_copied);
		};

	return search_recursive(key, parallel_compare);
}

template <typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
void ParallelZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::insert(const KEY_T* key) noexcept
{
	// Initialize parallel comparison state
	size_t comparison_size = KEY_T::MIN_PAR_COMPARE_CHAR_SIZE;
	size_t max_copied = 0;

	_buckets.push_back({ key, RANK_T::get_random() });
	unsigned new_node_index = _buckets.size() - 1;

	auto parallel_compare = [&](const KEY_T* k, const Bucket& v, const AncestorLCPs& ancestor_lcps)
		{
			return k_compare(k, v, ancestor_lcps, comparison_size, max_copied);
		};

	_root_index = insert_recursive(&_buckets[new_node_index], new_node_index, _root_index, {}, parallel_compare);
}
