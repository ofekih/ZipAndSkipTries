#pragma once

#include "BitString.hpp"

#include <compare>
#include <limits>
#include <vector>

#include <fstream>
#include <string>

struct MemoryEfficientLCP
{
	uint8_t exp_of_2 = 0;
	uint8_t multiple = 0;

	auto operator<=>(const MemoryEfficientLCP&) const = default;

	unsigned value() const noexcept
	{
		return std::pow(2, exp_of_2) * multiple;
	}

	std::ostream& print(std::ostream& os) const noexcept;
};

struct GeometricRank
{
	uint8_t geometric_rank;
	uint8_t uniform_rank;

	auto operator<=>(const GeometricRank&) const = default;

	static GeometricRank get_random()
	{
		static std::random_device rd;
		static std::default_random_engine generator(rd());
		// static std::default_random_engine generator(1);
		static std::geometric_distribution<uint8_t> g_dist(0.5);
		static std::uniform_int_distribution<uint8_t> u_dist(0, std::numeric_limits<uint8_t>::max());

		return { g_dist(generator), u_dist(generator) };
	}
};

template <typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T = GeometricRank, unsigned CHAR_SIZE_BITS = sizeof(CHAR_T) * 8>
class ZipTrie
{
public:
	using LCP_T = std::conditional_t<MEMORY_EFFICIENT, MemoryEfficientLCP, unsigned>;
	using KEY_T = BitString<CHAR_T, CHAR_SIZE_BITS>;

	ZipTrie(unsigned max_size, unsigned max_lcp_length);

	int get_depth(const KEY_T* key) const noexcept;
	int height() const noexcept;
	double get_average_height() const noexcept;
	unsigned size() const noexcept;
	bool contains(const KEY_T* key) const noexcept;
	LCP_T lcp(const KEY_T* key) const noexcept;

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

	const RANK_T& get_root_rank() const noexcept
	{
		return _buckets[_root_index].rank;
	}

	// used for testing only
	void set(const KEY_T* key, RANK_T rank, unsigned left, unsigned right, LCP_T predecessor_lcp, LCP_T successor_lcp) noexcept
	{
		_buckets.push_back({ key, rank, left, right, predecessor_lcp, successor_lcp });
	}

	// used for testing only
	void set_root_index(unsigned root_index) noexcept
	{
		_root_index = root_index;
	}

	struct SearchResults
	{
		bool contains;
		LCP_T max_lcp;
		int depth;
	};

	SearchResults search(const KEY_T* key) const noexcept;

	void to_dot(const std::string& file_path) const noexcept;

protected:
	unsigned _root_index;
	unsigned _max_lcp_length;
	uint8_t _log_max_size;

	static constexpr unsigned NULLPTR = std::numeric_limits<unsigned>::max();

	struct AncestorLCPs
	{
		LCP_T predecessor = {};
		LCP_T successor = {};
	};

	struct Bucket
	{
		const KEY_T* key;
		RANK_T rank;
		unsigned left = NULLPTR, right = NULLPTR;
		AncestorLCPs ancestor_lcps = {};
	};

	std::vector<Bucket> _buckets;

private:
	struct ComparisonResult
	{
		std::strong_ordering comparison;
		LCP_T lcp;
	};

	inline ComparisonResult compare(const KEY_T* x_key, const Bucket& v, const AncestorLCPs& ancestor_lcps) const noexcept;

	int height(unsigned node_index) const noexcept;
	uint64_t get_total_depth(unsigned node_index, uint64_t depth) const noexcept;

	
	unsigned insert_recursive(Bucket* x, unsigned x_index, unsigned node_index, AncestorLCPs ancestor_lcps = {}) noexcept;

	void to_dot_recursive(std::ofstream& os, unsigned node_index) const noexcept;

	inline MemoryEfficientLCP get_memory_efficient_lcp(size_t num) const noexcept;
};

template<typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::ZipTrie(unsigned max_size, unsigned max_lcp_length): _root_index(NULLPTR), _max_lcp_length(max_lcp_length), _log_max_size(std::log2(max_size))
{
	_buckets.reserve(max_size);
}

template<typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
bool ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::contains(const BitString<CHAR_T, CHAR_SIZE_BITS>* key) const noexcept
{
	return search(key).contains;
}

template<typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::LCP_T ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::lcp(const BitString<CHAR_T, CHAR_SIZE_BITS>* key) const noexcept
{
	return search(key).max_lcp;
}

template<typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
typename ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::ComparisonResult ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::compare(const KEY_T* x, const Bucket& v, const AncestorLCPs& ancestor_lcps) const noexcept
{
	auto predecessor_lcp = ancestor_lcps.predecessor;
	auto successor_lcp = ancestor_lcps.successor;

	auto x_max_lcp = std::max(predecessor_lcp, successor_lcp);
	auto corr_v_lcp = predecessor_lcp > successor_lcp ? v.ancestor_lcps.predecessor  : v.ancestor_lcps.successor;

	if (x_max_lcp != corr_v_lcp)
	{
		auto ordering = (x_max_lcp > corr_v_lcp) == (predecessor_lcp > successor_lcp) ? std::strong_ordering::less : std::strong_ordering::greater;
		return { ordering, std::min(x_max_lcp, corr_v_lcp) };
	}

	if constexpr (MEMORY_EFFICIENT)
	{
		auto [comparison, lcp] = x->compare(*v.key, x_max_lcp.value());
		return { comparison, get_memory_efficient_lcp(lcp) };
	}
	else
	{
		auto [comparison, lcp] = x->compare(*v.key, x_max_lcp);
		return { comparison, lcp };
	}
}

template<typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
typename ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::SearchResults ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::search(const BitString<CHAR_T, CHAR_SIZE_BITS>* key) const noexcept
{
	if (_buckets.empty())
	{
		return { false, 0, -1 };
	}

	unsigned v_index = _root_index;
	AncestorLCPs ancestor_lcps = {};
	int depth = 0;

	while (v_index != NULLPTR)
	{
		auto [ comparison, lcp ] = compare(key, _buckets[v_index], ancestor_lcps);

		if (comparison == std::strong_ordering::equal)
		{
			return { true, lcp, depth };
		}

		if (comparison == std::strong_ordering::less)
		{
			ancestor_lcps.successor = std::max(ancestor_lcps.successor, lcp);
			v_index = _buckets[v_index].left;
		}
		else
		{
			ancestor_lcps.predecessor = std::max(ancestor_lcps.predecessor, lcp);
			v_index = _buckets[v_index].right;
		}

		++depth;
	}

	return { false, std::max(ancestor_lcps.predecessor, ancestor_lcps.successor), -1 };
}

template <typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
void ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::insert(const BitString<CHAR_T, CHAR_SIZE_BITS>* key) noexcept
{
	_buckets.push_back({ key, RANK_T::get_random() });
	_root_index = insert_recursive(&_buckets.back(), size() - 1, _root_index);
}

template<typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
unsigned ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::insert_recursive(Bucket* x, unsigned x_index, unsigned v_index, AncestorLCPs ancestor_lcps) noexcept
{
	if (v_index == NULLPTR)
	{
		x->ancestor_lcps = ancestor_lcps;
		return x_index;
	}

	auto [ comparison, lcp ] = compare(x->key, _buckets[v_index], ancestor_lcps);

	if (comparison == std::strong_ordering::equal)
	{
		return v_index;
	}

	if (comparison == std::strong_ordering::less)
	{
		unsigned subroot_index = insert_recursive(x, x_index, _buckets[v_index].left, { ancestor_lcps.predecessor, std::max(ancestor_lcps.successor, lcp) });

		if (subroot_index == x_index && x->rank >= _buckets[v_index].rank)
		{
			_buckets[v_index].left = x->right;
			x->right = v_index;

			_buckets[v_index].ancestor_lcps.predecessor = lcp;
			x->ancestor_lcps = ancestor_lcps;

			return x_index;
		}
		else
		{
			_buckets[v_index].left = subroot_index;
		}
	}
	else
	{
		unsigned subroot_index = insert_recursive(x, x_index, _buckets[v_index].right, { std::max(ancestor_lcps.predecessor, lcp), ancestor_lcps.successor });

		if (subroot_index == x_index && x->rank >= _buckets[v_index].rank)
		{
			_buckets[v_index].right = x->left;
			x->left = v_index;

			_buckets[v_index].ancestor_lcps.successor = lcp;
			x->ancestor_lcps = ancestor_lcps;

			return x_index;
		}
		else
		{
			_buckets[v_index].right = subroot_index;
		}
	}

	return v_index;
}

template<typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
unsigned ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::size() const noexcept
{
	return _buckets.size();
}

template <typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
int ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::height() const noexcept
{
	return height(_root_index);
}

template <typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
int ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::height(unsigned node_index) const noexcept
{
	if (node_index == NULLPTR)
	{
		return -1;
	}

	return std::max(height(_buckets[node_index].left), height(_buckets[node_index].right)) + 1;
}

template <typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
int ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::get_depth(const BitString<CHAR_T, CHAR_SIZE_BITS>* key) const noexcept
{
	return search(key).depth;
}

template <typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
double ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::get_average_height() const noexcept
{
	return static_cast<double>(get_total_depth(_root_index, 0)) / size();
}

template <typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
uint64_t ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::get_total_depth(unsigned node_index, uint64_t depth) const noexcept
{
	if (node_index == NULLPTR)
	{
		return 0;
	}

	return get_total_depth(_buckets[node_index].left, depth + 1) + get_total_depth(_buckets[node_index].right, depth + 1) + depth;
}

template <typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
MemoryEfficientLCP ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::get_memory_efficient_lcp(size_t num) const noexcept
{
	if (num < _log_max_size)
	{
		return { 0, static_cast<uint8_t>(num) };
	}

	MemoryEfficientLCP lcp;
	lcp.exp_of_2 = std::log2(num / _log_max_size);
	lcp.multiple = num / std::pow(2, lcp.exp_of_2);

	return lcp;
}

template <typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
void ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::to_dot(const std::string& file_path) const noexcept
{
	std::ofstream os(file_path);
	os << "digraph G {\n";
	to_dot_recursive(os, _root_index);
	os << "}\n";
}

std::ostream& MemoryEfficientLCP::print(std::ostream& os) const noexcept
{
	return os << "2^" << static_cast<unsigned>(exp_of_2) << " * " << static_cast<unsigned>(multiple);
}

std::ostream& operator<<(std::ostream& os, const MemoryEfficientLCP& lcp)
{
	return lcp.print(os);
}

template <typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
void ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::to_dot_recursive(std::ofstream& os, unsigned node_index) const noexcept
{
	if (node_index == NULLPTR)
	{
		return;
	}

	const Bucket& node = _buckets[node_index];

	os << *node.key << " [label=\"" << *node.key << "\\n(" << node.ancestor_lcps.predecessor << ", " << node.ancestor_lcps.successor << ")\"];\n";

	if (node.left != NULLPTR)
	{
		const Bucket& left = _buckets[node.left];

		os << *node.key << " -> " << *left.key << ";\n";
		to_dot_recursive(os, node.left);
	}
	else
	{
		os << *node.key << "_left [label=\"Ø\"];\n";
		os << *node.key << " -> " << *node.key << "_left;\n";
	}

	if (node.right != NULLPTR)
	{
		const Bucket& right = _buckets[node.right];

		os << *node.key << " -> " << *right.key << ";\n";
		to_dot_recursive(os, node.right);
	}
	else
	{
		os << *node.key << "_right [label=\"Ø\"];\n";
		os << *node.key << " -> " << *node.key << "_right;\n";
	}
}
