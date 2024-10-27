#pragma once

#include <cstddef> // size_t
#include <cstdint> // uintmax_t

__global__ void par_xor_sig(const uintmax_t * const p1, uintmax_t *p2, uintmax_t *sig_out, size_t n);
__global__ void par_poptree(uintmax_t *tree, size_t n);
size_t seq_find_mismatch(const uintmax_t * const a, const uintmax_t * const b, size_t n);
size_t par_find_mismatch(const uintmax_t * const a, const uintmax_t * const b, size_t n);
size_t par_find_mismatch(const uintmax_t * const d_a, const uintmax_t * const b, uintmax_t *d_large_block, size_t n);
size_t par_find_mismatch(const uintmax_t * const d_a, uintmax_t *d_b, uintmax_t *d_tree, uintmax_t *d_section, size_t n);

size_t par_find_msw(uintmax_t *d_tree, uintmax_t *d_section, size_t num_leaves);
size_t seq_find_msw(const uintmax_t * const words, size_t n);

uintmax_t* alloc_large_block_to_device(size_t max_n);

uintmax_t* alloc_large_block_to_device_s(size_t max_n);
size_t par_find_mismatch_s(const uintmax_t * const d_a, const uintmax_t * const b, uintmax_t *d_large_block, size_t n);

