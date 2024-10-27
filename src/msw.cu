#include <bit> // bit_width
#include <cstddef> // size_t
#include <cstdint> // uintmax_t
#include <iostream>

#include <stdio.h>

#include "utility.cuh"

// TODO: Naming conventions in general suck, should fix

/*
__host__ inline size_t host_fast_log2(size_t n) 
{
	return std::bit_width(n) - 1;
}

__host__ inline size_t host_fast_log2_ceil(size_t n) 
{
	return std::bit_width(n - 1);
}

__device__ __forceinline__ size_t device_fast_log2(size_t n) 
{
	return std::numeric_limits<size_t>::digits - 1 - __clzll(n);
}

__device__ __forceinline__ size_t device_fast_log2_ceil(size_t n) 
{
	return std::numeric_limits<size_t>::digits - __clzll(n - 1);
}
*/

__device__ __forceinline__ size_t device_get_ancestor(size_t i, size_t height) 
{
	return ((i + 1) >> height) - 1;
}

__global__ void par_xor_sig(const uintmax_t * const p1, uintmax_t *p2, uintmax_t *sig_out, size_t n) 
{
	for (auto i = get_tid(); i < n; i += get_num_threads()) 
	{
		p2[i] = p1[i] ^ p2[i];

		if (p2[i] != 0)
			sig_out[i / device_fast_log2(n)] = 1;
	}
}

__global__ void par_poptree(uintmax_t *tree, size_t n) 
{
	auto num_ancestors = device_fast_log2_ceil(n);
	auto leaf_start = (1uLL << num_ancestors) - 1; // there are 2^num_ancestors - 1 internal nodes

	for (auto i = get_tid(); i < n * num_ancestors; i += get_num_threads()) 
	{
		auto leaf_index = i / num_ancestors;
		auto leaf_node = leaf_start + leaf_index;

		if (tree[leaf_node] == 0)
			continue;

		auto height_from_bottom = (i % num_ancestors) + 1;

		tree[device_get_ancestor(leaf_node, height_from_bottom)] = 1;
	}
}

__global__ void par_mark_msw(uintmax_t *tree, uintmax_t n) 
{
	auto num_ancestors = device_fast_log2_ceil(n);
	auto leaf_start = (1uLL << num_ancestors) - 1; // there are 2^num_ancestors - 1 internal nodes

	for (auto i = get_tid(); i < n * num_ancestors; i += get_num_threads()) 
	{
		auto leaf_index = i / num_ancestors;
		auto leaf_node = leaf_start + leaf_index;

		if (tree[leaf_node] == 0)
			continue;

		auto height_from_bottom = (i % num_ancestors) + 1;
		auto ancestor = device_get_ancestor(leaf_node, height_from_bottom);
		auto ancestor_left_child = (ancestor << 1) + 1;

		if (tree[ancestor_left_child] == 0)
			continue;

		auto one_lower_ancestor = device_get_ancestor(leaf_node, height_from_bottom - 1);

		if (one_lower_ancestor == ancestor_left_child)
			continue;

		tree[leaf_node] = 0;
	}
}

__global__ void par_get_section_diff(uintmax_t *tree, uintmax_t *section, size_t n) 
{
	for (auto i = get_tid(); i < n; i += get_num_threads())
		if (tree[i] != 0)
			*section = i;
}

size_t seq_find_mismatch(const uintmax_t * const arr1, const uintmax_t * const arr2, size_t size) 
{
	for (size_t i = 0; i < size; ++i)
		if (arr1[i] != arr2[i])
			return i;

	return size;
}

size_t seq_find_msw(const uintmax_t * const words, size_t n)
{
	for (size_t i = 0; i < n; ++i)
		if (words[i] != 0)
			return i;

	return n;
}

size_t par_find_msw(uintmax_t *d_tree, uintmax_t *d_section, size_t num_leaves)
{
	size_t padded_num_leaves = 1 << host_fast_log2_ceil(num_leaves);
	uintmax_t *leaves = d_tree + padded_num_leaves - 1;

	memset_within_device(d_section, 1, 0xFF);

	par_poptree<<<BLOCKS, THREADS>>>(d_tree, num_leaves);

	par_mark_msw<<<BLOCKS, THREADS>>>(d_tree, num_leaves);

	par_get_section_diff<<<BLOCKS, THREADS>>>(leaves, d_section, num_leaves);

	uintmax_t *h_section = copy_from_device(d_section, 1);
	size_t section = *h_section;

	delete[] h_section;

	if (section == std::numeric_limits<uintmax_t>::max())
		return num_leaves;

	return section;
}

// assume d_a and d_b are already populated
// assume d_tree and d_section are already allocated
size_t par_find_mismatch(const uintmax_t * const d_a, uintmax_t *d_b, uintmax_t *d_tree, uintmax_t *d_section, size_t n)
{
	std::cout << "par_find_mismatch" << std::endl;
	std::cout << n << std::endl;
	size_t num_leaves_per_section = host_fast_log2(n);
	size_t num_leaves = (n + num_leaves_per_section - 1) / num_leaves_per_section;
	size_t padded_num_leaves = 1 << host_fast_log2_ceil(num_leaves);
	uintmax_t *leaves = d_tree + padded_num_leaves - 1;

	memset_within_device(d_tree, num_leaves + padded_num_leaves - 1, 0);

	std::cout << "par_find_mismatch 1" << std::endl;

	// stores XOR in d_b
	par_xor_sig<<<BLOCKS, THREADS>>>(d_a, d_b, leaves, n);

	uintmax_t *h_leaves = copy_from_device(leaves, padded_num_leaves);
	std::cout << "leaves: ";
	for (size_t i = 0; i < padded_num_leaves; ++i)
		std::cout << h_leaves[i] << " ";

	std::cout << std::endl;

	std::cout << "par_find_mismatch 2" << std::endl;

	size_t h_section = par_find_msw(d_tree, d_section, num_leaves);
	if (h_section == num_leaves)
		return n;

	std::cout << "par_find_mismatch 3" << std::endl;

	size_t xor_offset = h_section * num_leaves_per_section;
	size_t xor_end = std::min(n, xor_offset + num_leaves_per_section);
	uintmax_t *xor_section = d_b + xor_offset;
	size_t num_words = xor_end - xor_offset;

	// print all the above variables, very pretty, fancy message

	std::cout << "h_section: " << h_section << std::endl;
	std::cout << "num_leaves_per_section: " << num_leaves_per_section << std::endl;
	std::cout << "xor_offset: " << xor_offset << std::endl;
	std::cout << "xor_end: " << xor_end << std::endl;
	std::cout << "num_words: " << num_words << std::endl;


	uintmax_t *h_xor_section = copy_from_device(xor_section, num_words);

	std::cout << "par_find_mismatch 4" << std::endl;

	size_t result = seq_find_msw(h_xor_section, num_words);
	delete[] h_xor_section;

	std::cout << "par_find_mismatch end" << std::endl;

	return xor_offset + result;
}

// assume d_a is already populated
// assume large block has size big enough to hold b, tree, and section (1)
// that size should be ~n + 2 * (n / log2(n)) + 1
size_t par_find_mismatch(const uintmax_t * const d_a, const uintmax_t * const b, uintmax_t *d_large_block, size_t n)
{
	uintmax_t *d_b = copy_to_device(d_large_block, b, n);
	uintmax_t *d_section = d_b + n;
	uintmax_t *d_tree = d_section + 1;

	return par_find_mismatch(d_a, d_b, d_tree, d_section, n);
}

uintmax_t* alloc_large_block_to_device(size_t n)
{
	size_t num_leaves_per_section = host_fast_log2(n);
	size_t num_leaves = (n + num_leaves_per_section - 1) / num_leaves_per_section;
	size_t padded_num_leaves = 1 << host_fast_log2_ceil(num_leaves);

	size_t large_block_size = n + num_leaves + padded_num_leaves;

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

__global__ void par_sig_sqrt(const uintmax_t * const p1, uintmax_t *p2, uintmax_t *sig_out, size_t n) {
	size_t num_sqrt = std::sqrt(n - 1) + 1;

	for (auto i = get_tid(); i < n; i += get_num_threads()) {
		p2[i] = p1[i] != p2[i];

		if (p2[i] != 0)
			sig_out[i / num_sqrt] = 1;
	}
}

size_t par_find_mismatch_sa(const uintmax_t * const d_a, const uintmax_t * const b, uintmax_t *d_large_block, size_t n)
{
	uintmax_t *d_b = copy_to_device(d_large_block, b, n);

	uintmax_t *d_section = d_b + n;
	uintmax_t *d_sqrt = d_section + 1;

	par_sig_sqrt<<<BLOCKS, THREADS>>>(d_a, d_b, d_sqrt, n);

	size_t num_sqrt = std::sqrt(n - 1) + 1;

	uintmax_t* h_sqrt = copy_from_device(d_sqrt, num_sqrt);

	size_t section = seq_find_msw(h_sqrt, num_sqrt);

	delete[] h_sqrt;

	if (section == num_sqrt)
		return n;

	size_t section_offset = section * num_sqrt;
	size_t section_end = std::min(n, section_offset + num_sqrt);
	uintmax_t *d_sig_section = d_b + section_offset;
	size_t num_words = section_end - section_offset;

	uintmax_t* h_sig_section = copy_from_device(d_sig_section, num_words);

	size_t result = seq_find_msw(h_sig_section, num_words);

	delete[] h_sig_section;

	return section_offset + result;
}

__global__ void naive_leftmost_prisoner(uintmax_t *p, size_t n)
{
	size_t n_choose_2 = n * (n - 1) / 2;
	size_t num_computations = (n_choose_2 - 1) / get_num_threads() + 1;

	auto i = get_tid() * num_computations;
	size_t b = (1 + std::sqrt(1 + 8 * i)) / 2;
	size_t a = i - b * (b - 1) / 2;

	while (num_computations-- && b < n)
	{
		if (p[a]) p[b] = 0;

		++a;

		if (a == b)
		{
			a = 0;
			++b;
		}
	}
}

size_t par_find_mismatch_s(const uintmax_t * const d_a, const uintmax_t * const b, uintmax_t *d_large_block, size_t n)
{
	uintmax_t *d_b = copy_to_device(d_large_block, b, n);

	uintmax_t *d_section = d_b + n;
	uintmax_t *d_sqrt = d_section + 1;
	size_t num_sqrt = std::sqrt(n - 1) + 1;

	memset_within_device(d_section, 1, 0xFF);

	memset_within_device(d_sqrt, num_sqrt, 0);

	par_sig_sqrt<<<BLOCKS, THREADS>>>(d_a, d_b, d_sqrt, n);

	naive_leftmost_prisoner<<<BLOCKS, THREADS>>>(d_sqrt, num_sqrt);

	par_get_section_diff<<<BLOCKS, THREADS>>>(d_sqrt, d_section, num_sqrt);

	uintmax_t *h_section = copy_from_device(d_section, 1);
	size_t section = *h_section;
	
	if (section == std::numeric_limits<uintmax_t>::max()) {
		delete[] h_section;
		return n;
	}

	memset_within_device(d_section, 1, 0xFF);

	size_t section_offset = section * num_sqrt;
	size_t section_end = std::min(n, section_offset + num_sqrt);
	uintmax_t *d_sig_section = d_b + section_offset;
	size_t num_words = section_end - section_offset;

	naive_leftmost_prisoner<<<BLOCKS, THREADS>>>(d_sig_section, num_words);

	par_get_section_diff<<<BLOCKS, THREADS>>>(d_sig_section, d_section, num_words);

	h_section = copy_from_device(d_section, 1);
	size_t result = *h_section;

	delete[] h_section;

	return section_offset + result;
}

uintmax_t* alloc_large_block_to_device_s(size_t n)
{
	size_t large_block_size = n + std::sqrt(n - 1) + 2;

	return alloc_to_device<uintmax_t>(large_block_size);
}
