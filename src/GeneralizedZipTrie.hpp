#pragma once

#include "BitString.hpp"

#include <limits>
#include <vector>

template <typename CHAR_T, typename RANK_T, typename LCP_T, unsigned CHAR_SIZE_BITS = sizeof(CHAR_T) * 8>
class GeneralizedZipTrie
{
public:
	GeneralizedZipTrie(unsigned max_size, LCP_T max_lcp_length);

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

protected:
	unsigned _root_index;
	LCP_T _max_lcp_length;

	static constexpr unsigned NULLPTR = std::numeric_limits<unsigned>::max();

	struct Bucket
	{
		const BitString<CHAR_T, CHAR_SIZE_BITS>* key;
		RANK_T rank;
		unsigned left = NULLPTR, right = NULLPTR;
		LCP_T parent_lcp = 0, spine_lcp = 0;
	};

	std::vector<Bucket> _buckets;

	virtual RANK_T get_random_rank() const noexcept = 0;

private:
	int height(unsigned nodeIndex) const noexcept;
	uint64_t get_total_depth(unsigned nodeIndex, uint64_t depth) const noexcept;
};

template<typename CHAR_T, typename RANK_T, typename LCP_T, unsigned CHAR_SIZE_BITS>
GeneralizedZipTrie<CHAR_T, RANK_T, LCP_T, CHAR_SIZE_BITS>::GeneralizedZipTrie(unsigned max_size, LCP_T max_lcp_length): _root_index(NULLPTR), _max_lcp_length(max_lcp_length)
{
	_buckets.reserve(max_size);
}

template<typename CHAR_T, typename RANK_T, typename LCP_T, unsigned CHAR_SIZE_BITS>
bool GeneralizedZipTrie<CHAR_T, RANK_T, LCP_T, CHAR_SIZE_BITS>::contains(const BitString<CHAR_T, CHAR_SIZE_BITS>* key) const noexcept
{
	return contains_lcp(key).contains;
}

template<typename CHAR_T, typename RANK_T, typename LCP_T, unsigned CHAR_SIZE_BITS>
LCP_T GeneralizedZipTrie<CHAR_T, RANK_T, LCP_T, CHAR_SIZE_BITS>::lcp(const BitString<CHAR_T, CHAR_SIZE_BITS>* key) const noexcept
{
	return contains_lcp(key).max_lcp;
}

template<typename CHAR_T, typename RANK_T, typename LCP_T, unsigned CHAR_SIZE_BITS>
typename GeneralizedZipTrie<CHAR_T, RANK_T, LCP_T, CHAR_SIZE_BITS>::ContainsLCP GeneralizedZipTrie<CHAR_T, RANK_T, LCP_T, CHAR_SIZE_BITS>::contains_lcp(const BitString<CHAR_T, CHAR_SIZE_BITS>* key) const noexcept
{
	if (_buckets.empty())
	{
		return { false, 0 };
	}

	LCP_T delta = 0;

	// the below two booleans don't really matter for initialization, since delta is 0
	bool lcp_shared_with_parent = true;
	bool smaller_than_parent = false;

	unsigned curr_index = _root_index;

	while (curr_index != NULLPTR)
	{
		const auto& curr = _buckets[curr_index];
		LCP_T relevant_lcp = lcp_shared_with_parent ? curr.parent_lcp : curr.spine_lcp;

		if (delta == relevant_lcp) // we must actually compare when equal
		{
			auto [comparison, lcp_length] = key->compare(*curr.key, delta);

			if (comparison == std::strong_ordering::equal)
			{
				return { true, delta };
			}

			delta = lcp_length;
			smaller_than_parent = comparison == std::strong_ordering::less;
			lcp_shared_with_parent = true;
		}
		else
		{
			smaller_than_parent = smaller_than_parent == (lcp_shared_with_parent == (delta < relevant_lcp));
			lcp_shared_with_parent = delta < relevant_lcp;
		}

		curr_index = smaller_than_parent ? curr.left : curr.right;
	}

	return { false, delta };
}


// template <typename KeyType, typename RANK_T>
// void GeneralizedZipTrie<KeyType, RANK_T>::insert(const KeyType& key) noexcept
// {
// 	Bucket x = { key, get_random_rank(&_total_comparisons, &_firstTies, &_bothTies) };
// 	unsigned xIndex = _buckets.size();

// 	if (xIndex == 0)
// 	{
// 		_root_index = xIndex;
// 		_buckets.emplace_back(x);
// 		return;
// 	}

// 	auto& rank = x.rank;

// 	unsigned curr_index = _root_index;
// 	unsigned prevIndex = NULLPTR;

// 	while (curr_index != NULLPTR && (rank < _buckets[curr_index].rank || (rank == _buckets[curr_index].rank && key > _buckets[curr_index].key)))
// 	{
// 		prevIndex = curr_index;
// 		curr_index = key < _buckets[curr_index].key ? _buckets[curr_index].left : _buckets[curr_index].right;
// 	}

// 	_buckets.emplace_back(x);

// 	if (curr_index == _root_index)
// 	{
// 		_root_index = xIndex;
// 	}
// 	else if (key < _buckets[prevIndex].key)
// 	{
// 		_buckets[prevIndex].left = xIndex;
// 	}
// 	else
// 	{
// 		_buckets[prevIndex].right = xIndex;
// 	}

// 	if (curr_index == NULLPTR)
// 	{
// 		return;
// 	}

// 	if (key < _buckets[curr_index].key)
// 	{
// 		_buckets[xIndex].right = curr_index;
// 	}
// 	else
// 	{
// 		_buckets[xIndex].left = curr_index;
// 	}

// 	prevIndex = xIndex;

// 	while (curr_index != NULLPTR)
// 	{
// 		unsigned fixIndex = prevIndex;

// 		if (_buckets[curr_index].key < key)
// 		{
// 			do
// 			{
// 				prevIndex = curr_index;
// 				curr_index = _buckets[curr_index].right;
// 			}
// 			while (curr_index != NULLPTR && _buckets[curr_index].key < key);
// 		}
// 		else
// 		{
// 			do
// 			{
// 				prevIndex = curr_index;
// 				curr_index = _buckets[curr_index].left;
// 			}
// 			while (curr_index != NULLPTR && _buckets[curr_index].key > key);
// 		}

// 		if (_buckets[fixIndex].key > key || (fixIndex == xIndex && _buckets[prevIndex].key > key))
// 		{
// 			_buckets[fixIndex].left = curr_index;
// 		}
// 		else
// 		{
// 			_buckets[fixIndex].right = curr_index;
// 		}
// 	}
// }


// template <typename KeyType, typename RANK_T>
// unsigned GeneralizedZipTrie<KeyType, RANK_T>::size() const noexcept
// {
// 	return _buckets.size();
// }

// template <typename KeyType, typename RANK_T>
// int GeneralizedZipTrie<KeyType, RANK_T>::height() const noexcept
// {
// 	return height(_root_index);
// }

// template <typename KeyType, typename RANK_T>
// int GeneralizedZipTrie<KeyType, RANK_T>::height(unsigned nodeIndex) const noexcept
// {
// 	if (nodeIndex == NULLPTR)
// 	{
// 		return -1;
// 	}

// 	return std::max(height(_buckets[nodeIndex].left), height(_buckets[nodeIndex].right)) + 1;
// }

// template <typename KeyType, typename RANK_T>
// int GeneralizedZipTrie<KeyType, RANK_T>::get_depth(const KeyType& key) const noexcept
// {
// 	unsigned curr_index = _root_index;
// 	int depth = 0;

// 	while (curr_index != NULLPTR)
// 	{
// 		if (key < _buckets[curr_index].key)
// 		{
// 			curr_index = _buckets[curr_index].left;
// 		}
// 		else if (_buckets[curr_index].key < key)
// 		{
// 			curr_index = _buckets[curr_index].right;
// 		}
// 		else
// 		{
// 			return depth;
// 		}

// 		++depth;
// 	}

// 	return -1;
// }

// template <typename KeyType, typename RANK_T>
// double GeneralizedZipTrie<KeyType, RANK_T>::get_average_height() const noexcept
// {
// 	return static_cast<double>(get_total_depth(_root_index, 0)) / size();
// }

// template <typename KeyType, typename RANK_T>
// uint64_t GeneralizedZipTrie<KeyType, RANK_T>::get_total_depth(unsigned nodeIndex, uint64_t depth) const noexcept
// {
// 	if (nodeIndex == NULLPTR)
// 	{
// 		return 0;
// 	}

// 	return get_total_depth(_buckets[nodeIndex].left, depth + 1) + get_total_depth(_buckets[nodeIndex].right, depth + 1) + depth;
// }
