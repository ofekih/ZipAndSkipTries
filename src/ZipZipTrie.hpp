#pragma once

#include "GeneralizedZipTrie.hpp"

#include <limits>
#include <random>

struct GeometricRank
{
	uint8_t geometric_rank;
	uint8_t uniform_rank;

	auto operator<=>(const GeometricRank&) const = default;
};;

template <typename CHAR_T, bool MEMORY_EFFICIENT, unsigned CHAR_SIZE_BITS = sizeof(CHAR_T) * 8>
class ZipZipTrie : public GeneralizedZipTrie<CHAR_T, GeometricRank, MEMORY_EFFICIENT, CHAR_SIZE_BITS>
{
public:
	ZipZipTrie(unsigned max_size, unsigned max_lcp_length) : GeneralizedZipTrie<CHAR_T, GeometricRank, MEMORY_EFFICIENT, CHAR_SIZE_BITS>(max_size, max_lcp_length) {}

protected:
	GeometricRank get_random_rank() const noexcept override
	{
		static std::random_device rd;
		// static std::default_random_engine generator(rd());
		static std::default_random_engine generator(1);
		static std::geometric_distribution<uint8_t> g_dist(0.5);
		static std::uniform_int_distribution<uint8_t> u_dist(0, std::numeric_limits<uint8_t>::max());

		return { g_dist(generator), u_dist(generator) };
	}
};
