/**
 * @file ZipTrie.hpp
 * @brief Defines the ZipTrie class, a trie structure combined with Zip Tree properties for efficient string storage and retrieval.
 */
#pragma once

#include "BitString.cuh" // Assumes BitString class definition is here

#include <compare>   // For std::strong_ordering
#include <limits>    // For std::numeric_limits
#include <vector>    // For std::vector
#include <cmath>     // For std::log2
#include <random>    // For random rank generation
#include <fstream>   // For file output (to_dot)
#include <string>    // For std::string
#include <type_traits> // For std::conditional_t
#include <ostream>   // For std::ostream
#include <algorithm> // For std::max, std::min

/**
 * @struct MemoryEfficientLCP
 * @brief Represents an approximate Longest Common Prefix (LCP) in a memory-efficient format.
 * Stores the LCP value as multiple * 2^exp_of_2.
 */
struct MemoryEfficientLCP
{
	/** @brief The exponent part of the LCP representation (2^exp_of_2). */
	uint8_t exp_of_2 = 0;
	/** @brief The multiple part of the LCP representation. */
	uint8_t multiple = 0;

	/** @brief Default comparison operator. */
	auto operator<=>(const MemoryEfficientLCP&) const = default;

	/**
	 * @brief Calculates the actual LCP value.
	 * @return unsigned The calculated LCP value.
	 */
	unsigned value() const noexcept
	{
		return (1u << exp_of_2) * multiple;
	}

	/**
	 * @brief Prints the LCP representation to an output stream.
	 * @param os The output stream.
	 * @return std::ostream& Reference to the output stream.
	 */
	std::ostream& print(std::ostream& os) const noexcept;
};

/**
 * @struct GeometricRank
 * @brief Represents a rank used in Zip Trees, typically generated randomly.
 * Combines a geometric and a uniform random component for tie-breaking.
 */
struct GeometricRank
{
	/** @brief Rank component typically following a geometric distribution. */
	uint8_t geometric_rank;
	/** @brief Rank component typically following a uniform distribution (for tie-breaking). */
	uint8_t uniform_rank;

	/** @brief Default comparison operator. Compares geometric_rank first, then uniform_rank. */
	auto operator<=>(const GeometricRank&) const = default;

	/**
	 * @brief Generates a random rank using static random number generators.
	 * @return GeometricRank A randomly generated rank.
	 * @note Uses static generators, making it potentially not thread-safe if called concurrently.
	 */
	static GeometricRank get_random()
	{
		// Static generators ensure the sequence is consistent across calls within a run,
		// but also mean it's not thread-safe without external locking.
		static std::random_device rd;
		static std::default_random_engine generator(rd());
		// static std::default_random_engine generator(1); // For deterministic testing
		static std::geometric_distribution<uint8_t> g_dist(0.5); // p=0.5 means ~50% chance of 0, ~25% chance of 1, etc.
		static std::uniform_int_distribution<uint8_t> u_dist(0, std::numeric_limits<uint8_t>::max());

		return { g_dist(generator), u_dist(generator) };
	}
};

/**
 * @class ZipTrie
 * @brief Implements a Zip Trie, combining properties of Tries and Zip Trees.
 *
 * A Zip Trie stores keys (represented by BitString) and associated ranks.
 * It maintains tree structure based on both key prefixes (like a Trie) and
 * random ranks (like a Zip Tree/Treap), aiming for balanced expected depth.
 *
 * @tparam CHAR_T The character type used in the keys (BitString).
 * @tparam MEMORY_EFFICIENT If true, uses `MemoryEfficientLCP` for LCP storage, otherwise uses `unsigned`.
 * @tparam RANK_T The type used for node ranks (default: `GeometricRank`). Must provide `get_random()` and comparison.
 * @tparam CHAR_SIZE_BITS The number of significant bits per character in `CHAR_T`.
 */
template <typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T = GeometricRank, unsigned CHAR_SIZE_BITS = sizeof(CHAR_T) * 8>
class ZipTrie
{
public:
	/** @brief The type used to store Longest Common Prefix values. Conditional based on `MEMORY_EFFICIENT`. */
	using LCP_T = std::conditional_t<MEMORY_EFFICIENT, MemoryEfficientLCP, unsigned>;
	/** @brief The type used for keys, based on the BitString template. */
	using KEY_T = BitString<CHAR_T, CHAR_SIZE_BITS>;

	/**
	 * @brief Constructs a ZipTrie.
	 * @param max_size Hint for reserving memory for the expected maximum number of nodes.
	 * @param max_lcp_length The maximum possible LCP length, used for memory-efficient LCP calculation.
	 */
	ZipTrie(unsigned max_size, unsigned max_lcp_length);

	/**
	 * @brief Calculates the depth of a given key in the trie.
	 * @param key Pointer to the key to search for.
	 * @return int The depth of the key (root is at depth 0), or -1 if the key is not found.
	 */
	int get_depth(const KEY_T* key) const noexcept;

	/**
	 * @brief Calculates the height of the trie.
	 * @return int The height of the trie (longest path from root to a leaf), or -1 for an empty trie.
	 */
	int height() const noexcept;

	/**
	 * @brief Calculates the average depth of all nodes in the trie.
	 * @return double The average depth. Returns NaN or 0 for an empty trie.
	 */
	double get_average_height() const noexcept;

	/**
	 * @brief Returns the number of nodes currently stored in the trie.
	 * @return unsigned The number of nodes.
	 */
	unsigned size() const noexcept;

	/**
	 * @brief Checks if a key exists in the trie.
	 * @param key Pointer to the key to search for.
	 * @return bool True if the key is found, false otherwise.
	 */
	bool contains(const KEY_T* key) const noexcept;

	/**
	 * @brief Finds the Longest Common Prefix (LCP) between the search key and the path taken in the trie.
	 * If the key is found, returns the LCP with the found key. If not found, returns the LCP with the
	 * path leading to the insertion point.
	 * @param key Pointer to the key to search for.
	 * @return LCP_T The calculated LCP value.
	 */
	LCP_T lcp(const KEY_T* key) const noexcept;

	/**
	 * @brief Inserts a key into the zip trie.
	 *
	 * Assigns a random rank to the new node and inserts it based on zip tree and trie properties.
	 * Assumes keys are unique; inserting duplicate keys leads to undefined behavior.
	 *
	 * @param key Pointer to the new key to insert. Must be non-null and point to a valid key.
	 * The ZipTrie stores this pointer directly; the key's lifetime must exceed the trie's.
	 */
	void insert(const KEY_T* key) noexcept;

	/**
	 * @brief Removes a node with a given key from the zip tree. (Currently commented out)
	 *
	 * @param  key key of node to remove
	 * @return     true if a node was removed, false otherwise
	 */
	// bool remove(const KEY_T& key) noexcept; // @note Functionality not implemented.

	/**
	 * @brief Gets the rank of the root node.
	 * @return const RANK_T& A constant reference to the root node's rank. Behavior is undefined if the trie is empty.
	 */
	const RANK_T& get_root_rank() const noexcept
	{
		// Undefined behavior if _root_index == NULLPTR
		return _buckets[_root_index].rank;
	}

	/**
	 * @brief Directly adds a bucket to the internal storage.
	 * @note Intended for testing purposes only. Bypasses normal insertion logic.
	 * @param key Pointer to the key.
	 * @param rank The rank to assign.
	 * @param left Index of the left child.
	 * @param right Index of the right child.
	 * @param predecessor_lcp LCP with the predecessor ancestor path.
	 * @param successor_lcp LCP with the successor ancestor path.
	 */
	void set(const KEY_T* key, RANK_T rank, unsigned left, unsigned right, LCP_T predecessor_lcp, LCP_T successor_lcp) noexcept
	{
		_buckets.push_back({ key, rank, left, right, {predecessor_lcp, successor_lcp} });
	}

	/**
	 * @brief Directly sets the root index of the trie.
	 * @note Intended for testing purposes only.
	 * @param root_index The index to set as the root.
	 */
	void set_root_index(unsigned root_index) noexcept
	{
		_root_index = root_index;
	}

	/**
	 * @struct SearchResults
	 * @brief Holds the results of a search operation in the ZipTrie.
	 */
	struct SearchResults
	{
		/** @brief True if the key was found in the trie, false otherwise. */
		bool contains;
		/** @brief The Longest Common Prefix (LCP) calculated during the search. See `lcp()` method documentation. */
		LCP_T max_lcp;
		/** @brief The depth at which the key was found, or -1 if not found. */
		int depth;
	};

	/**
	 * @brief Performs a search for a key in the trie.
	 * @param key Pointer to the key to search for.
	 * @return SearchResults A struct containing the search outcome (found status, LCP, depth).
	 */
	SearchResults search(const KEY_T* key) const noexcept;

	/**
	 * @brief Generates a DOT representation of the trie for visualization.
	 * @param file_path The path to the output DOT file.
	 */
	void to_dot(const std::string& file_path) const noexcept;

protected:
	/** @brief Index of the root node in the `_buckets` vector. `NULLPTR` if the trie is empty. */
	unsigned _root_index;
	/** @brief Stores the maximum possible LCP length, used for memory-efficient LCP calculation. */
	unsigned _max_lcp_length; // TODO: Currently unused, potentially remove or integrate?
	/** @brief Stores the approximate log base 2 of the maximum expected size, used for memory-efficient LCP calculation. */
	uint8_t _log_max_size;

	/** @brief Represents a null child pointer index. */
	static constexpr unsigned NULLPTR = std::numeric_limits<unsigned>::max();

	/**
	 * @struct AncestorLCPs
	 * @brief Stores the LCP values associated with the path from the root to a node's parent.
	 * These represent the LCPs with the inorder predecessor and successor paths discovered so far.
	 */
	struct AncestorLCPs
	{
		/** @brief LCP with the path corresponding to the inorder predecessor. */
		LCP_T predecessor = {};
		/** @brief LCP with the path corresponding to the inorder successor. */
		LCP_T successor = {};
	};

	/**
	 * @struct Bucket
	 * @brief Represents a node in the ZipTrie. Stores key, rank, children indices, and ancestor LCPs.
	 */
	struct Bucket
	{
		/** @brief Pointer to the key associated with this node. Lifetime managed externally. */
		const KEY_T* key;
		/** @brief The rank associated with this node (used for zip property). */
		RANK_T rank;
		/** @brief Index of the left child in `_buckets`, or `NULLPTR`. */
		unsigned left = NULLPTR;
		/** @brief Index of the right child in `_buckets`, or `NULLPTR`. */
		unsigned right = NULLPTR;
		/** @brief LCPs inherited from the parent node during insertion/search. */
		AncestorLCPs ancestor_lcps = {};
	};

	/** @brief Vector storing all the nodes (buckets) of the trie. */
	std::vector<Bucket> _buckets;

private:
	/**
	 * @struct ComparisonResult
	 * @brief Holds the result of comparing a search key `x` with a node `v`.
	 */
	struct ComparisonResult
	{
		/** @brief The comparison ordering (less, greater, equal). */
		std::strong_ordering comparison;
		/** @brief The Longest Common Prefix calculated during the comparison. */
		LCP_T lcp;
	};

	/**
	 * @brief Compares a key `x_key` with the key in a node `v`, considering ancestor LCPs.
	 *
	 * This is a key function for navigating the Zip Trie. It determines the relative order
	 * and LCP based on previously computed LCPs along the ancestor path, potentially avoiding
	 * full string comparisons.
	 *
	 * @param x_key Pointer to the key being searched for or inserted.
	 * @param v Constant reference to the current node (bucket) being compared against.
	 * @param ancestor_lcps LCPs computed along the path from the root to `v`'s parent.
	 * @return ComparisonResult The result of the comparison (ordering and LCP).
	 */
	inline ComparisonResult k_compare(const KEY_T* x_key, const Bucket& v, const AncestorLCPs& ancestor_lcps) const noexcept;

	/**
	 * @brief Recursive helper function to calculate the height of a subtree.
	 * @param node_index Index of the root of the subtree.
	 * @return int The height of the subtree, or -1 if `node_index` is `NULLPTR`.
	 */
	int height(unsigned node_index) const noexcept;

	/**
	 * @brief Recursive helper function to calculate the sum of depths of all nodes in a subtree.
	 * @param node_index Index of the root of the subtree.
	 * @param depth The depth of the current `node_index`.
	 * @return uint64_t The sum of depths of all nodes in the subtree rooted at `node_index`.
	 */
	uint64_t get_total_depth(unsigned node_index, uint64_t depth) const noexcept;


	/**
	 * @brief Recursive helper function for inserting a new node `x`.
	 *
	 * Traverses the trie based on key comparisons and LCPs. When the correct position
	 * is found, it inserts the node `x`. During backtracking, it performs rotations (zips)
	 * if the rank of `x` is higher than the current node `v`, maintaining the zip property.
	 *
	 * @param x Pointer to the new bucket being inserted.
	 * @param x_index Index of the new bucket `x` in the `_buckets` vector.
	 * @param v_index Index of the current node being visited in the recursion.
	 * @param ancestor_lcps LCPs inherited from the path above `v_index`.
	 * @return unsigned The index of the root of the modified subtree (could be `x_index` or `v_index`).
	 */
	unsigned insert_recursive(Bucket* x, unsigned x_index, unsigned v_index, AncestorLCPs ancestor_lcps = {}) noexcept;

	/**
	 * @brief Recursive helper function to generate DOT output for a subtree.
	 * @param os The output file stream.
	 * @param node_index Index of the root of the subtree to output.
	 */
	void to_dot_recursive(std::ofstream& os, unsigned node_index) const noexcept;

	/**
	 * @brief Converts a standard LCP length (`size_t`) into the memory-efficient format.
	 * @param num The LCP length to convert.
	 * @return MemoryEfficientLCP The LCP value in the memory-efficient representation.
	 * @note This function is only relevant when `MEMORY_EFFICIENT` is true.
	 */
	inline MemoryEfficientLCP get_memory_efficient_lcp(size_t num) const noexcept;
};

//-----------------------------------------------------------------------------
// ZipTrie Method Implementations
//-----------------------------------------------------------------------------

/**
 * @brief Constructor implementation. Initializes members and reserves bucket capacity.
 */
template<typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::ZipTrie(unsigned max_size, unsigned max_lcp_length)
	: _root_index(NULLPTR), _max_lcp_length(max_lcp_length), _log_max_size(max_size > 0 ? static_cast<uint8_t>(std::log2(max_size)) : 0) // Avoid log2(0)
{
	_buckets.reserve(max_size);
}

/**
 * @brief Implementation of `contains`. Calls `search` and returns the boolean result.
 */
template<typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
bool ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::contains(const BitString<CHAR_T, CHAR_SIZE_BITS>* key) const noexcept
{
	return search(key).contains;
}

/**
 * @brief Implementation of `lcp`. Calls `search` and returns the LCP result.
 */
template<typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::LCP_T ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::lcp(const BitString<CHAR_T, CHAR_SIZE_BITS>* key) const noexcept
{
	return search(key).max_lcp;
}

/**
 * @brief Implementation of the key comparison logic using ancestor LCPs.
 */
template<typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
typename ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::ComparisonResult ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::k_compare(const KEY_T* x, const Bucket& v, const AncestorLCPs& ancestor_lcps) const noexcept
{
	// Determine the relevant LCPs from the ancestor path and the current node v
	auto predecessor_lcp = ancestor_lcps.predecessor;
	auto successor_lcp = ancestor_lcps.successor;

	// Find the maximum LCP from the path taken so far to reach v's parent
	auto x_max_lcp = std::max(predecessor_lcp, successor_lcp);
	// Find the corresponding LCP stored *at* node v that relates to the path taken
	// If the predecessor path had the max LCP, use v's predecessor LCP, otherwise use v's successor LCP.
	auto corr_v_lcp = predecessor_lcp > successor_lcp ? v.ancestor_lcps.predecessor : v.ancestor_lcps.successor;

	// If the path LCP (x_max_lcp) doesn't match the LCP stored at v (corr_v_lcp),
	// we can determine the order without further string comparison.
	if (x_max_lcp != corr_v_lcp)
	{
		// Determine order based on which LCP was larger and which path (predecessor/successor) it belonged to.
		// If x_max_lcp > corr_v_lcp, and the path taken was the predecessor path (pred > succ), then x < v.key.
		// If x_max_lcp > corr_v_lcp, and the path taken was the successor path (pred <= succ), then x > v.key.
		// And vice versa if x_max_lcp < corr_v_lcp.
		return {
			(x_max_lcp > corr_v_lcp) == (predecessor_lcp > successor_lcp) ? std::strong_ordering::less : std::strong_ordering::greater,
			std::min(x_max_lcp, corr_v_lcp) // The actual LCP is the minimum of the two differing LCPs.
		};
	}

	// If the path LCP matches the node's stored LCP, we need to perform actual string comparison
	// starting from that known common prefix length (x_max_lcp).
	if constexpr (MEMORY_EFFICIENT)
	{
		// Use the .value() method if LCPs are memory-efficient objects
		auto [comparison, lcp_val] = x->seq_k_compare(*v.key, x_max_lcp.value());
		// Convert the resulting LCP length back to the memory-efficient format
		return { comparison, get_memory_efficient_lcp(lcp_val) };
	}
	else
	{
		// Use the LCP value directly if it's just an unsigned integer
		auto [comparison, lcp_val] = x->seq_k_compare(*v.key, x_max_lcp);
		return { comparison, static_cast<LCP_T>(lcp_val) };
	}
}

/**
 * @brief Implementation of the search logic. Traverses the trie based on `k_compare`.
 */
template<typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
typename ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::SearchResults ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::search(const BitString<CHAR_T, CHAR_SIZE_BITS>* key) const noexcept
{
	if (_root_index == NULLPTR) // Check if trie is empty
	{
		return { false, {}, -1 }; // Return default LCP_T
	}

	unsigned v_index = _root_index;
	AncestorLCPs ancestor_lcps = {}; // Start with zero ancestor LCPs at root
	int depth = 0;

	while (v_index != NULLPTR)
	{
		const Bucket& v_node = _buckets[v_index];
		auto [ comparison, lcp ] = k_compare(key, v_node, ancestor_lcps);

		// Key found exactly
		if (comparison == std::strong_ordering::equal)
		{
			// The LCP with the found key is the one calculated in the final k_compare
			return { true, lcp, depth };
		}

		// Navigate left or right based on comparison, updating the relevant ancestor LCP
		if (comparison == std::strong_ordering::less)
		{
			// Going left: the current node v becomes a potential successor boundary.
			// The LCP with this boundary is `lcp`. Update the successor path LCP.
			ancestor_lcps.successor = std::max(ancestor_lcps.successor, lcp);
			v_index = v_node.left;
		}
		else // comparison == std::strong_ordering::greater
		{
			// Going right: the current node v becomes a potential predecessor boundary.
			// The LCP with this boundary is `lcp`. Update the predecessor path LCP.
			ancestor_lcps.predecessor = std::max(ancestor_lcps.predecessor, lcp);
			v_index = v_node.right;
		}

		++depth;
	}

	// Key not found. The max LCP is the maximum of the predecessor and successor LCPs
	// accumulated along the search path.
	return { false, std::max(ancestor_lcps.predecessor, ancestor_lcps.successor), -1 };
}

/**
 * @brief Implementation of insert. Creates a node and calls the recursive helper.
 */
template <typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
void ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::insert(const BitString<CHAR_T, CHAR_SIZE_BITS>* key) noexcept
{
	// Add the new bucket to the vector. It gets a random rank.
	_buckets.push_back({ key, RANK_T::get_random() });
	unsigned new_node_index = size() - 1;
	// Call the recursive insertion function starting from the root.
	_root_index = insert_recursive(&_buckets[new_node_index], new_node_index, _root_index);
}

/**
 * @brief Implementation of the recursive insertion logic with zipping.
 */
template<typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
unsigned ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::insert_recursive(Bucket* x, unsigned x_index, unsigned v_index, AncestorLCPs ancestor_lcps) noexcept
{
	// Base case: Found the insertion point (null link)
	if (v_index == NULLPTR)
	{
		// Store the final ancestor LCPs calculated along the path in the new node.
		x->ancestor_lcps = ancestor_lcps;
		return x_index; // Return the index of the newly inserted node.
	}

	Bucket& v_node = _buckets[v_index]; // Get mutable reference to current node v

	// Compare the new key x->key with the current node v_node's key
	auto [ comparison, lcp ] = k_compare(x->key, v_node, ancestor_lcps);

	// Key already exists (or is equivalent based on comparison) - stop insertion.
	if (comparison == std::strong_ordering::equal)
	{
		// Note: Behavior for duplicate keys is undefined. This simply stops recursion.
		// Could potentially update value or handle differently if needed.
		// Need to potentially deallocate or mark the pushed-back bucket x as unused.
		// Current implementation might leave an unused bucket if duplicates are inserted.
		return v_index; // Return index of existing node
	}

	// Navigate left
	if (comparison == std::strong_ordering::less)
	{
		// Recursively insert into the left subtree.
		// Update the successor ancestor LCP for the recursive call.
		unsigned subroot_index = insert_recursive(x, x_index, v_node.left, { ancestor_lcps.predecessor, std::max(ancestor_lcps.successor, lcp) });

		// Check if the new node x became the root of the subtree returned by the recursive call,
		// AND if x's rank requires it to be "zipped" up past the current node v.
		if (subroot_index == x_index && x->rank >= v_node.rank)
		{
			// Perform rotation (Right rotation around v)
			v_node.left = x->right; // v adopts x's right child
			x->right = v_index;     // x becomes parent of v

			// Update ancestor LCPs stored *at the nodes* based on the rotation
			v_node.ancestor_lcps.predecessor = lcp; // v's predecessor LCP is now the LCP calculated during comparison
			x->ancestor_lcps = ancestor_lcps;       // x inherits the original ancestor LCPs from above v

			return x_index; // x is the new root of this subtree
		}
		else
		{
			// No rotation needed, just update v's left child pointer
			v_node.left = subroot_index;
		}
	}
	// Navigate right
	else // comparison == std::strong_ordering::greater
	{
		// Recursively insert into the right subtree.
		// Update the predecessor ancestor LCP for the recursive call.
		unsigned subroot_index = insert_recursive(x, x_index, v_node.right, { std::max(ancestor_lcps.predecessor, lcp), ancestor_lcps.successor });

		// Check if the new node x became the root of the subtree returned by the recursive call,
		// AND if x's rank requires it to be "zipped" up past the current node v.
		if (subroot_index == x_index && x->rank >= v_node.rank)
		{
			// Perform rotation (Left rotation around v)
			v_node.right = x->left; // v adopts x's left child
			x->left = v_index;      // x becomes parent of v

			// Update ancestor LCPs stored *at the nodes* based on the rotation
			v_node.ancestor_lcps.successor = lcp; // v's successor LCP is now the LCP calculated during comparison
			x->ancestor_lcps = ancestor_lcps;     // x inherits the original ancestor LCPs from above v

			return x_index; // x is the new root of this subtree
		}
		else
		{
			// No rotation needed, just update v's right child pointer
			v_node.right = subroot_index;
		}
	}

	// If no rotation occurred, v remains the root of this subtree
	return v_index;
}

/**
 * @brief Returns the number of nodes.
 */
template<typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
unsigned ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::size() const noexcept
{
	return _buckets.size();
}

/**
 * @brief Calculates the height of the trie by calling the recursive helper.
 */
template <typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
int ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::height() const noexcept
{
	return height(_root_index);
}

/**
 * @brief Recursive height calculation.
 */
template <typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
int ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::height(unsigned node_index) const noexcept
{
	if (node_index == NULLPTR)
	{
		return -1; // Height of an empty tree is -1
	}

	// Height is 1 + max height of subtrees
	return std::max(height(_buckets[node_index].left), height(_buckets[node_index].right)) + 1;
}

/**
 * @brief Gets the depth of a key by calling `search`.
 */
template <typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
int ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::get_depth(const BitString<CHAR_T, CHAR_SIZE_BITS>* key) const noexcept
{
	return search(key).depth;
}

/**
 * @brief Calculates the average height (average depth) of nodes.
 */
template <typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
double ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::get_average_height() const noexcept
{
	unsigned current_size = size();
	if (current_size == 0)
	{
		return 0.0; // Or NaN, depending on desired behavior for empty trie
	}
	// Calculate sum of depths using recursive helper, then divide by size.
	return static_cast<double>(get_total_depth(_root_index, 0)) / current_size;
}

/**
 * @brief Recursive calculation of the sum of depths of all nodes in a subtree.
 */
template <typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
uint64_t ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::get_total_depth(unsigned node_index, uint64_t depth) const noexcept
{
	if (node_index == NULLPTR)
	{
		return 0;
	}

	// Sum of depths = current node's depth + sum in left subtree + sum in right subtree
	return get_total_depth(_buckets[node_index].left, depth + 1) + get_total_depth(_buckets[node_index].right, depth + 1) + depth;
}

/**
 * @brief Converts a size_t LCP value to an approximate MemoryEfficientLCP.
 */
template <typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
MemoryEfficientLCP ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::get_memory_efficient_lcp(size_t num) const noexcept
{
	if (num < _log_max_size)
	{
		return { 0, static_cast<uint8_t>(num) };
	}

	// Otherwise, calculate exponent and multiple.
	MemoryEfficientLCP lcp;
	lcp.exp_of_2 = std::log2(num / _log_max_size);
	lcp.multiple = num / (1u << lcp.exp_of_2);

	return lcp;
}

/**
 * @brief Generates the DOT file header and calls the recursive helper.
 */
template <typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
void ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::to_dot(const std::string& file_path) const noexcept
{
	std::ofstream os(file_path);
	if (!os) {
		// Handle error: Failed to open file
		std::cerr << "Error: Could not open file " << file_path << " for writing." << std::endl;
		return;
	}
	os << "digraph G {\n";
	os << "node [shape=record];\n"; // Use record shape for better labeling potentially
	to_dot_recursive(os, _root_index);
	os << "}\n";
}

//-----------------------------------------------------------------------------
// Standalone Function Implementations
//-----------------------------------------------------------------------------

/**
 * @brief Implementation of the print method for MemoryEfficientLCP.
 */
inline std::ostream& MemoryEfficientLCP::print(std::ostream& os) const noexcept
{
	return os << "2^" << static_cast<unsigned>(exp_of_2) << " * " << static_cast<unsigned>(multiple);
}

/**
 * @brief Overloaded stream insertion operator for MemoryEfficientLCP.
 * @param os The output stream.
 * @param lcp The MemoryEfficientLCP object to print.
 * @return std::ostream& Reference to the output stream.
 */
inline std::ostream& operator<<(std::ostream& os, const MemoryEfficientLCP& lcp)
{
	return lcp.print(os);
}

/**
 * @brief Implementation of the recursive DOT file generation.
 */
template <typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
void ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::to_dot_recursive(std::ofstream& os, unsigned node_index) const noexcept
{
	if (node_index == NULLPTR)
	{
		return;
	}

	const Bucket& node = _buckets[node_index];

	// Define the node appearance in DOT format
	// Using pointer value as a unique ID in the DOT file to handle potential duplicate keys visually
	os << "  node_" << node_index << " [label=\"{Key: " << *node.key
	   << "|Rank: (" << static_cast<unsigned>(node.rank.geometric_rank) << "," << static_cast<unsigned>(node.rank.uniform_rank) << ")"
	   << "|LCPs: (" << node.ancestor_lcps.predecessor << ", " << node.ancestor_lcps.successor << ")}\"];\n";

	// Define edges to children
	if (node.left != NULLPTR)
	{
		const Bucket& left = _buckets[node.left];
		os << "  node_" << node_index << " -> node_" << node.left << " [label=\" L\"];\n";
		to_dot_recursive(os, node.left);
	}
	else
	{
		// Optional: Draw null pointers explicitly
		os << "  null_left_" << node_index << " [shape=point];\n";
		// os << "  node_" << node_index << " -> null_left_" << node_index << " [label=\" L\"];\n";
	}

	if (node.right != NULLPTR)
	{
		const Bucket& right = _buckets[node.right];
		os << "  node_" << node_index << " -> node_" << node.right << " [label=\" R\"];\n";
		to_dot_recursive(os, node.right);
	}
	else
	{
		// Optional: Draw null pointers explicitly
		// os << "  null_right_" << node_index << " [shape=point];\n";
		// os << "  node_" << node_index << " -> null_right_" << node_index << " [label=\" R\"];\n";
	}
}
