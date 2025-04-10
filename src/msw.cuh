/**
 * @file msw.cuh
 * @brief CUDA-accelerated Most Significant Word finding algorithms.
 * @details Provides functions for finding the most significant word (MSW) and
 * the first mismatching word between two arrays of words. The implementation uses
 * GPU-accelerated algorithms with a square root decomposition approach for efficient
 * parallel processing.
 */

#pragma once

#include <cstddef> // size_t
#include <cstdint> // uintmax_t

/**
 * @brief Sequential implementation of finding the first mismatching word between two arrays.
 * @details Compares two arrays element by element until finding the first mismatch.
 * This is the baseline CPU implementation.
 *
 * @param a First input array (host memory)
 * @param b Second input array (host memory)
 * @param n Number of elements in the arrays
 * @return Index of the first mismatching word, or n if arrays are identical
 */
size_t seq_find_mismatch(const uintmax_t * const a, const uintmax_t * const b, size_t n);

/**
 * @brief GPU-accelerated implementation of finding the first mismatching word.
 * @details Allocates device memory, copies data, and uses parallel algorithms to 
 * find the first mismatch. Uses a square root decomposition approach that divides 
 * the arrays into √n chunks for efficient parallel processing.
 *
 * @param a First input array (host memory)
 * @param b Second input array (host memory)
 * @param n Number of elements in the arrays
 * @return Index of the first mismatching word, or n if arrays are identical
 */
size_t par_find_mismatch(const uintmax_t * const a, const uintmax_t * const b, size_t n);

/**
 * @brief GPU-accelerated mismatch finding with pre-allocated device memory.
 * @details Uses a square root decomposition approach: first divides the input into √n chunks,
 * identifies the first chunk with a mismatch in O(n) time and constant span,
 * then finds the exact mismatch within that chunk, also in O(n) time and constant span.
 * 
 * @param d_a First input array (already in device memory)
 * @param b Second input array (host memory, will be copied to device)
 * @param d_large_block Pre-allocated device memory block for temporary storage
 * @param n Number of elements in the arrays
 * @return Index of the first mismatching word, or n if arrays are identical
 */
size_t par_find_mismatch(const uintmax_t * const d_a, const uintmax_t * const b, uintmax_t *d_large_block, size_t n);

/**
 * @brief Allocates a large block of device memory for parallel mismatch operations.
 * @details Allocates enough memory for the square root decomposition algorithm,
 * which requires space for the second array, section information, and √n additional elements.
 *
 * @param max_n Maximum number of elements to be processed
 * @return Pointer to allocated device memory
 */
uintmax_t* alloc_large_block_to_device(size_t max_n);
