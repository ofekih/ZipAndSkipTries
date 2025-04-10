#include <bit> // bit_width
#include <cstddef> // size_t
#include <cstdint> // uintmax_t
#include <iostream>

#include <stdio.h>

#include "cuda_utils.cuh"

/**
 * @brief CUDA kernel that implements the first phase of the square root decomposition.
 * @details This kernel performs two critical operations:
 *   1. Computes the XOR (as inequality check) between corresponding elements of two arrays
 *   2. Maps each mismatch to its respective sqrt(n)-sized chunk in the signature array
 *
 * Each thread processes multiple elements based on its thread ID. When a mismatch
 * is found, the corresponding chunk in sig_out is marked with a 1.
 */
__global__ void par_sig_sqrt(const uintmax_t * const p1, uintmax_t *p2, uintmax_t *sig_out, size_t n) {
	// Calculate number of chunks, each of size sqrt(n)
	size_t num_sqrt = std::sqrt(n - 1) + 1;

	// Each thread processes elements at positions (threadIdx + k*blockDim*gridDim)
	for (auto i = get_tid(); i < n; i += get_num_threads()) {
		// Store XOR result in p2 - non-zero values indicate mismatches
		p2[i] = p1[i] != p2[i];

		// If there's a mismatch, mark its chunk in the signature array
		if (p2[i] != 0)
			sig_out[i / num_sqrt] = 1;
	}
}

/**
 * @brief Parallel kernel to find the leftmost non-zero element in an array.
 * @details Uses a tournament-style approach where each pair of elements is compared,
 * and if a left element (lower index) is non-zero, it "eliminates" the right element
 * by setting it to zero. After all comparisons, only the leftmost non-zero element remains.
 * This kernel achieves O(n²) work in O(1) span when sufficient threads are available.
 *
 * The algorithm distributes all n-choose-2 pairwise comparisons across available threads.
 */
__global__ void naive_leftmost_prisoner(uintmax_t *p, size_t n)
{
	// Calculate total number of pairwise comparisons needed
	size_t n_choose_2 = n * (n - 1) / 2;
	
	// Distribute comparisons evenly across threads
	size_t num_computations = (n_choose_2 - 1) / get_num_threads() + 1;

	// Calculate starting comparison index for this thread
	auto i = get_tid() * num_computations;
	
	// Convert linear index i to pair (a,b) using inverse of triangular number formula
	// This maps from a 1D index to a 2D comparison between elements at positions a and b
	size_t b = (1 + std::sqrt(1 + 8 * i)) / 2;  // Row index calculation
	size_t a = i - b * (b - 1) / 2;             // Column index calculation
	
	// Process this thread's allocated comparisons
	while (num_computations-- && b < n)
	{
		// If the left element is non-zero, eliminate the right element
		if (p[a]) p[b] = 0;

		// Move to next comparison pair
		++a;

		// If we reach the end of a row, move to the next row
		if (a == b)
		{
			a = 0;  // Start at the leftmost element
			++b;    // Move to the next row
		}
	}
}

/**
 * @brief Parallel kernel to identify and record the index of the first non-zero element.
 * @details After the naive_leftmost_prisoner kernel has eliminated all non-leftmost elements,
 * this kernel scans the array to find the remaining non-zero element and records its index.
 * Multiple threads may write to the same location, but they'll all write the same value
 * since only one element should remain non-zero after the elimination process.
 * 
 * This kernel is used in both phases of the square root algorithm:
 * 1. To find which sqrt(n)-sized chunk contains the first mismatch
 * 2. To find the exact position within that chunk
 */
__global__ void par_get_section_diff(uintmax_t *tree, uintmax_t *section, size_t n) 
{
	// Each thread checks a subset of elements
	for (auto i = get_tid(); i < n; i += get_num_threads())
		// If we find a non-zero element, it must be the leftmost one after elimination
		if (tree[i] != 0)
			*section = i;  // Record its index as the result
}

size_t seq_find_mismatch(const uintmax_t * const arr1, const uintmax_t * const arr2, size_t size) 
{
	for (size_t i = 0; i < size; ++i)
		if (arr1[i] != arr2[i])
			return i;

	return size;
}

// assume d_a is already populated on the device
// assume large block has size big enough to hold b and signature array
// that size should be ~n + sqrt(n) + 2
size_t par_find_mismatch(const uintmax_t * const d_a, const uintmax_t * const b, uintmax_t *d_large_block, size_t n)
{
	// Copy the second array to device memory using pre-allocated buffer
	uintmax_t *d_b = copy_to_device(d_large_block, b, n);

	// Set up memory pointers within the large block:
	// - d_section: single value to store result
	// - d_sqrt: array of size sqrt(n) to track which chunks have mismatches
	uintmax_t *d_section = d_b + n;
	uintmax_t *d_sqrt = d_section + 1;
	size_t num_sqrt = std::sqrt(n - 1) + 1;  // Calculate number of sqrt-sized chunks

	// Initialize the section result to max value (no match found yet)
	memset_within_device(d_section, 1, 0xFF);

	// Initialize the sqrt array to zeros (no mismatches found yet)
	memset_within_device(d_sqrt, num_sqrt, 0);

	// PHASE 1: Find which sqrt(n)-sized chunk has the first mismatch
	// For each element, compute XOR of the two arrays and mark its chunk in d_sqrt if mismatch
	par_sig_sqrt<<<BLOCKS, THREADS>>>(d_a, d_b, d_sqrt, n);

	// Find the first non-zero chunk (chunk with a mismatch) in parallel
	naive_leftmost_prisoner<<<BLOCKS, THREADS>>>(d_sqrt, num_sqrt);

	// Extract the index of the first chunk with a mismatch
	par_get_section_diff<<<BLOCKS, THREADS>>>(d_sqrt, d_section, num_sqrt);

	// Copy the result back to host memory
	uintmax_t *h_section = copy_from_device(d_section, 1);
	size_t section = *h_section;
	
	// If no mismatch found, return n (indicating arrays are identical)
	if (section == std::numeric_limits<uintmax_t>::max()) {
		delete[] h_section;
		return n;
	}

	// Reset the section result for the second phase
	memset_within_device(d_section, 1, 0xFF);

	// PHASE 2: Find the exact mismatch position within the identified chunk
	// Calculate the offset and size of the identified chunk
	size_t section_offset = section * num_sqrt;
	size_t section_end = std::min(n, section_offset + num_sqrt);
	uintmax_t *d_sig_section = d_b + section_offset;
	size_t num_words = section_end - section_offset;

	// Find the first non-zero word within the chunk (containing XOR results)
	naive_leftmost_prisoner<<<BLOCKS, THREADS>>>(d_sig_section, num_words);

	// Extract the index of the first mismatching word within the chunk
	par_get_section_diff<<<BLOCKS, THREADS>>>(d_sig_section, d_section, num_words);

	// Copy the result back to host memory
	h_section = copy_from_device(d_section, 1);
	size_t result = *h_section;

	delete[] h_section;

	// Return the global index of the first mismatch (chunk offset + position within chunk)
	return section_offset + result;
}

uintmax_t* alloc_large_block_to_device(size_t n)
{
	size_t large_block_size = n + std::sqrt(n - 1) + 2;

	return alloc_to_device<uintmax_t>(large_block_size);
}

size_t par_find_mismatch(const uintmax_t * const a, const uintmax_t * const b, size_t n) 
{
	uintmax_t *d_a = copy_to_device(a, n);
	uintmax_t *d_large_block = alloc_large_block_to_device(n);

	size_t result = par_find_mismatch(d_a, b, d_large_block, n);

	device_free(d_a);
	device_free(d_large_block);

	return result;
}
