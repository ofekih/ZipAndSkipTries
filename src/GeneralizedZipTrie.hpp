#pragma once

#include "BitString.hpp"

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

template <typename CHAR_T, typename RANK_T, bool MEMORY_EFFICIENT, unsigned CHAR_SIZE_BITS = sizeof(CHAR_T) * 8>
class GeneralizedZipTrie
{
public:
	using LCP_T = std::conditional_t<MEMORY_EFFICIENT, MemoryEfficientLCP, unsigned>;

	GeneralizedZipTrie(unsigned max_size, unsigned max_lcp_length);

	int get_depth(const BitString<CHAR_T, CHAR_SIZE_BITS>* key) const noexcept;
	int height() const noexcept;
	double get_average_height() const noexcept;
	unsigned size() const noexcept;
	bool contains(const BitString<CHAR_T, CHAR_SIZE_BITS>* key) const noexcept;
	LCP_T lcp(const BitString<CHAR_T, CHAR_SIZE_BITS>* key) const noexcept;

	/**
	 * Inserts a key, value pair into the zip tree. Note that inserting there is
	 * no validation that the keys don't already exist. Add only unique keys to
	 * avoid undefined behavior.
	 *
	 * @param key new node key
	 * @param val new node value
	 */
	void insert(const BitString<CHAR_T, CHAR_SIZE_BITS>* key) noexcept;

	/**
	 * Removes a node with a given key from the zip tree.
	 *
	 * @param  key key of node to remove
	 * @return     true if a node was removed, false otherwise
	 */
	// bool remove(const BitString<CHAR_T, CHAR_SIZE_BITS>& key) noexcept;

	const RANK_T& get_root_rank() const noexcept
	{
		return _buckets[_root_index].rank;
	}

	// used for testing only
	void set(const BitString<CHAR_T, CHAR_SIZE_BITS>* key, RANK_T rank, unsigned left, unsigned right, LCP_T parent_lcp, LCP_T spine_lcp) noexcept
	{
		_buckets.push_back({ key, rank, left, right, parent_lcp, spine_lcp });
	}

	// used for testing only
	void set_root_index(unsigned root_index) noexcept
	{
		_root_index = root_index;
	}

	struct ContainsLCP
	{
		bool contains;
		LCP_T max_lcp;
	};

	ContainsLCP contains_lcp(const BitString<CHAR_T, CHAR_SIZE_BITS>* key) const noexcept;

	void to_dot(const std::string& file_path) const noexcept;

protected:
	unsigned _root_index;
	unsigned _max_lcp_length;
	uint8_t _log_max_size;

	static constexpr unsigned NULLPTR = std::numeric_limits<unsigned>::max();

	struct Bucket
	{
		const BitString<CHAR_T, CHAR_SIZE_BITS>* key;
		RANK_T rank;
		unsigned left = NULLPTR, right = NULLPTR;
		LCP_T parent_lcp = {}, spine_lcp = {};
	};

	std::vector<Bucket> _buckets;

	virtual RANK_T get_random_rank() const noexcept = 0;

private:
	struct TraversalKnowledge
	{
		unsigned index;
		LCP_T delta = {};
		LCP_T lcp_shared = {};
		bool lcp_shared_with_parent = false;
		std::strong_ordering comparison = std::strong_ordering::less;
	};

	inline TraversalKnowledge compare(const BitString<CHAR_T, CHAR_SIZE_BITS>* key, TraversalKnowledge tk) const noexcept;

	int height(unsigned node_index) const noexcept;
	uint64_t get_total_depth(unsigned node_index, uint64_t depth) const noexcept;

	struct Ancestors
	{
		unsigned parent = NULLPTR;
		unsigned spine = NULLPTR;
	};

	struct Memory
	{
		Ancestors left_ancestors, right_ancestors;
		LCP_T left_lcp_shared_with_x, right_lcp_shared_with_x;
	};
	
	unsigned insert_recursive(Bucket* x, unsigned x_index, const TraversalKnowledge& tk, const Ancestors& _ancestors, Memory& memory, unsigned& spine_index) noexcept;

	void to_dot_recursive(std::ofstream& os, unsigned node_index) const noexcept;

	inline MemoryEfficientLCP get_memory_efficient_lcp(size_t num) const noexcept;
};

template<typename CHAR_T, typename RANK_T, bool MEMORY_EFFICIENT, unsigned CHAR_SIZE_BITS>
GeneralizedZipTrie<CHAR_T, RANK_T, MEMORY_EFFICIENT, CHAR_SIZE_BITS>::GeneralizedZipTrie(unsigned max_size, unsigned max_lcp_length): _root_index(NULLPTR), _max_lcp_length(max_lcp_length), _log_max_size(std::log2(max_size))
{
	_buckets.reserve(max_size);
}

template<typename CHAR_T, typename RANK_T, bool MEMORY_EFFICIENT, unsigned CHAR_SIZE_BITS>
bool GeneralizedZipTrie<CHAR_T, RANK_T, MEMORY_EFFICIENT, CHAR_SIZE_BITS>::contains(const BitString<CHAR_T, CHAR_SIZE_BITS>* key) const noexcept
{
	return contains_lcp(key).contains;
}

template<typename CHAR_T, typename RANK_T, bool MEMORY_EFFICIENT, unsigned CHAR_SIZE_BITS>
GeneralizedZipTrie<CHAR_T, RANK_T, MEMORY_EFFICIENT, CHAR_SIZE_BITS>::LCP_T GeneralizedZipTrie<CHAR_T, RANK_T, MEMORY_EFFICIENT, CHAR_SIZE_BITS>::lcp(const BitString<CHAR_T, CHAR_SIZE_BITS>* key) const noexcept
{
	return contains_lcp(key).max_lcp;
}

template<typename CHAR_T, typename RANK_T, bool MEMORY_EFFICIENT, unsigned CHAR_SIZE_BITS>
typename GeneralizedZipTrie<CHAR_T, RANK_T, MEMORY_EFFICIENT, CHAR_SIZE_BITS>::TraversalKnowledge GeneralizedZipTrie<CHAR_T, RANK_T, MEMORY_EFFICIENT, CHAR_SIZE_BITS>::compare(const BitString<CHAR_T, CHAR_SIZE_BITS>* key, TraversalKnowledge tk) const noexcept
{
	const Bucket& curr = _buckets[tk.index];
	LCP_T relevant_lcp = tk.lcp_shared_with_parent ? curr.parent_lcp : curr.spine_lcp;

	if (tk.delta == relevant_lcp) // we must actually compare when equal
	{
		tk.lcp_shared_with_parent = true;

		if constexpr (MEMORY_EFFICIENT)
		{
			auto [comparison, delta] = key->compare(*curr.key, tk.delta.value());
			tk.delta = get_memory_efficient_lcp(delta);
			tk.comparison = comparison;
		}
		else
		{
			auto [comparison, delta] = key->compare(*curr.key, tk.delta);
			tk.delta = delta;
			tk.comparison = comparison;
		}

		tk.lcp_shared = tk.delta;

		if (tk.comparison == std::strong_ordering::equal)
		{
			return tk;
		}
	}
	else
	{
		tk.comparison = ((tk.comparison == std::strong_ordering::less) == (tk.lcp_shared_with_parent == tk.delta < relevant_lcp)) ? std::strong_ordering::less : std::strong_ordering::greater;
		tk.lcp_shared_with_parent = tk.delta < relevant_lcp;
		tk.lcp_shared = std::min(tk.delta, relevant_lcp);
	}

	tk.index = tk.comparison == std::strong_ordering::less ? curr.left : curr.right;

	return tk;
}

template<typename CHAR_T, typename RANK_T, bool MEMORY_EFFICIENT, unsigned CHAR_SIZE_BITS>
typename GeneralizedZipTrie<CHAR_T, RANK_T, MEMORY_EFFICIENT, CHAR_SIZE_BITS>::ContainsLCP GeneralizedZipTrie<CHAR_T, RANK_T, MEMORY_EFFICIENT, CHAR_SIZE_BITS>::contains_lcp(const BitString<CHAR_T, CHAR_SIZE_BITS>* key) const noexcept
{
	if (_buckets.empty())
	{
		return { false, 0 };
	}

	TraversalKnowledge tk = { _root_index };

	while (tk.index != NULLPTR)
	{
		tk = compare(key, tk);

		if (tk.comparison == std::strong_ordering::equal)
		{
			return { true, tk.delta };
		}
	}

	return { false, tk.delta };
}

template <typename CHAR_T, typename RANK_T, bool MEMORY_EFFICIENT, unsigned CHAR_SIZE_BITS>
void GeneralizedZipTrie<CHAR_T, RANK_T, MEMORY_EFFICIENT, CHAR_SIZE_BITS>::insert(const BitString<CHAR_T, CHAR_SIZE_BITS>* key) noexcept
{
	_buckets.push_back({ key, get_random_rank() });
	unsigned x_index = size() - 1;
	Bucket& x = _buckets.back();
	Memory memory;
	unsigned spine_index = NULLPTR;
	_root_index = insert_recursive(&x, x_index, { _root_index }, {}, memory, spine_index);

	if (_root_index == x_index)
	{
		if (x.left != NULLPTR)
		{
			_buckets[x.left].parent_lcp = memory.left_lcp_shared_with_x;
			_buckets[x.left].spine_lcp = {};
		}

		if (x.right != NULLPTR)
		{
			_buckets[x.right].parent_lcp = memory.right_lcp_shared_with_x;
			_buckets[x.right].spine_lcp = {};
		}
	}
}

template<typename CHAR_T, typename RANK_T, bool MEMORY_EFFICIENT, unsigned CHAR_SIZE_BITS>
unsigned GeneralizedZipTrie<CHAR_T, RANK_T, MEMORY_EFFICIENT, CHAR_SIZE_BITS>::insert_recursive(Bucket* x, unsigned x_index, const TraversalKnowledge& tk, const Ancestors& ancestors, Memory& memory, unsigned& spine_index) noexcept
{
	if (tk.index == NULLPTR)
	{
		return x_index;
	}

	auto next_tk = compare(x->key, tk);

	if (next_tk.comparison == std::strong_ordering::equal)
	{
		return tk.index;
	}

	Bucket* root = &_buckets[tk.index];

	unsigned subroot_index = insert_recursive(x, x_index, next_tk, { tk.index, next_tk.comparison == tk.comparison ? ancestors.spine : ancestors.parent }, memory, spine_index);
	Bucket* subroot = &_buckets[subroot_index];

	if (spine_index == tk.index)
	{
		x->spine_lcp = next_tk.lcp_shared;
	}

	if (next_tk.comparison == std::strong_ordering::less)
	{
		if (subroot_index == x_index && x->rank >= root->rank) // x is new root, we are in P
		{
			unsigned x_right_index = x->right;
			if (x_right_index != NULLPTR)
			{
				auto& x_right_bucket = _buckets[x_right_index];
				const auto& x_right_ancestors = memory.right_ancestors;

				if (tk.index != x_right_ancestors.parent)
				{
					x_right_bucket.parent_lcp = x_right_bucket.spine_lcp;
				}

				x_right_bucket.spine_lcp = memory.right_lcp_shared_with_x;
			}

			root->left = x->right;
			x->right = tk.index;

			memory.right_ancestors = ancestors;
			memory.right_lcp_shared_with_x = next_tk.lcp_shared;

			return x_index;
		}
		else
		{
			root->left = subroot_index;
		}
	}
	else
	{
		if (subroot_index == x_index && x->rank > root->rank) // x is new root, we are in P
		{
			unsigned x_left_index = x->left;
			if (x_left_index != NULLPTR)
			{
				auto& x_left_bucket = _buckets[x_left_index];
				const auto& x_left_ancestors = memory.left_ancestors;

				if (tk.index != x_left_ancestors.parent)
				{
					x_left_bucket.parent_lcp = x_left_bucket.spine_lcp;
				}

				x_left_bucket.spine_lcp = memory.left_lcp_shared_with_x;
			}

			root->right = x->left;
			x->left = tk.index;

			memory.left_ancestors = ancestors;
			memory.left_lcp_shared_with_x = next_tk.lcp_shared;

			return x_index;
		}
		else
		{
			root->right = subroot_index;
		}
	}

	if (subroot_index == x_index)
	{
		x->parent_lcp = next_tk.lcp_shared;

		spine_index = tk.comparison == next_tk.comparison ? ancestors.spine : ancestors.parent;

		unsigned left_index = x->left;
		if (left_index != NULLPTR)
		{
			auto& left_bucket = _buckets[left_index];
			const auto& left_ancestors = memory.left_ancestors;

			if (next_tk.comparison == std::strong_ordering::greater && tk.index != left_ancestors.spine)
			{
				left_bucket.spine_lcp = left_bucket.parent_lcp;
			}

			left_bucket.parent_lcp = memory.left_lcp_shared_with_x;
		}

		unsigned right_index = x->right;
		if (right_index != NULLPTR)
		{
			auto& right_bucket = _buckets[right_index];
			const auto& right_ancestors = memory.right_ancestors;

			if (next_tk.comparison == std::strong_ordering::less && tk.index != right_ancestors.spine)
			{
				right_bucket.spine_lcp = right_bucket.parent_lcp;
			}

			right_bucket.parent_lcp = memory.right_lcp_shared_with_x;
		}
	}

	return tk.index;
}

template<typename CHAR_T, typename RANK_T, bool MEMORY_EFFICIENT, unsigned CHAR_SIZE_BITS>
unsigned GeneralizedZipTrie<CHAR_T, RANK_T, MEMORY_EFFICIENT, CHAR_SIZE_BITS>::size() const noexcept
{
	return _buckets.size();
}

template <typename CHAR_T, typename RANK_T, bool MEMORY_EFFICIENT, unsigned CHAR_SIZE_BITS>
int GeneralizedZipTrie<CHAR_T, RANK_T, MEMORY_EFFICIENT, CHAR_SIZE_BITS>::height() const noexcept
{
	return height(_root_index);
}

template <typename CHAR_T, typename RANK_T, bool MEMORY_EFFICIENT, unsigned CHAR_SIZE_BITS>
int GeneralizedZipTrie<CHAR_T, RANK_T, MEMORY_EFFICIENT, CHAR_SIZE_BITS>::height(unsigned node_index) const noexcept
{
	if (node_index == NULLPTR)
	{
		return -1;
	}

	return std::max(height(_buckets[node_index].left), height(_buckets[node_index].right)) + 1;
}

template <typename CHAR_T, typename RANK_T, bool MEMORY_EFFICIENT, unsigned CHAR_SIZE_BITS>
int GeneralizedZipTrie<CHAR_T, RANK_T, MEMORY_EFFICIENT, CHAR_SIZE_BITS>::get_depth(const BitString<CHAR_T, CHAR_SIZE_BITS>* key) const noexcept
{
	unsigned curr_index = _root_index;
	int depth = 0;

	while (curr_index != NULLPTR)
	{
		if (*key < *_buckets[curr_index].key)
		{
			curr_index = _buckets[curr_index].left;
		}
		else if (*_buckets[curr_index].key < *key)
		{
			curr_index = _buckets[curr_index].right;
		}
		else
		{
			return depth;
		}

		++depth;
	}

	return -1;
}

template <typename CHAR_T, typename RANK_T, bool MEMORY_EFFICIENT, unsigned CHAR_SIZE_BITS>
double GeneralizedZipTrie<CHAR_T, RANK_T, MEMORY_EFFICIENT, CHAR_SIZE_BITS>::get_average_height() const noexcept
{
	return static_cast<double>(get_total_depth(_root_index, 0)) / size();
}

template <typename CHAR_T, typename RANK_T, bool MEMORY_EFFICIENT, unsigned CHAR_SIZE_BITS>
uint64_t GeneralizedZipTrie<CHAR_T, RANK_T, MEMORY_EFFICIENT, CHAR_SIZE_BITS>::get_total_depth(unsigned node_index, uint64_t depth) const noexcept
{
	if (node_index == NULLPTR)
	{
		return 0;
	}

	return get_total_depth(_buckets[node_index].left, depth + 1) + get_total_depth(_buckets[node_index].right, depth + 1) + depth;
}

template <typename CHAR_T, typename RANK_T, bool MEMORY_EFFICIENT, unsigned CHAR_SIZE_BITS>
MemoryEfficientLCP GeneralizedZipTrie<CHAR_T, RANK_T, MEMORY_EFFICIENT, CHAR_SIZE_BITS>::get_memory_efficient_lcp(size_t num) const noexcept
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

template <typename CHAR_T, typename RANK_T, bool MEMORY_EFFICIENT, unsigned CHAR_SIZE_BITS>
void GeneralizedZipTrie<CHAR_T, RANK_T, MEMORY_EFFICIENT, CHAR_SIZE_BITS>::to_dot(const std::string& file_path) const noexcept
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

template <typename CHAR_T, typename RANK_T, bool MEMORY_EFFICIENT, unsigned CHAR_SIZE_BITS>
void GeneralizedZipTrie<CHAR_T, RANK_T, MEMORY_EFFICIENT, CHAR_SIZE_BITS>::to_dot_recursive(std::ofstream& os, unsigned node_index) const noexcept
{
	if (node_index == NULLPTR)
	{
		return;
	}

	const Bucket& node = _buckets[node_index];

	os << *node.key << " [label=\"" << *node.key << "\\n(" << node.spine_lcp << ")\"];\n";

	if (node.left != NULLPTR)
	{
		const Bucket& left = _buckets[node.left];

		os << *node.key << " -> " << *left.key << " [label=\"" << left.parent_lcp << "\"];\n";
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

		os << *node.key << " -> " << *right.key << " [label=\"" << right.parent_lcp << "\"];\n";
		to_dot_recursive(os, node.right);
	}
	else
	{
		os << *node.key << "_right [label=\"Ø\"];\n";
		os << *node.key << " -> " << *node.key << "_right;\n";
	}
}

