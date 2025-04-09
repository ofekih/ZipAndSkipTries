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
#include <optional>  // For std::optional
#include <string>    // For std::string
#include <type_traits> // For std::conditional_t
#include <ostream>   // For std::ostream
#include <algorithm> // For std::max, std::min
#include <iostream>  // For std::cerr (used in to_dot error handling)

/**
 * @struct MemoryEfficientLCP
 * @brief Represents an approximate Longest Common Prefix (LCP) in a memory-efficient format.
 * Stores the LCP value as multiple * 2^exp_of_2. This allows representing larger LCP values
 * with less precision using fewer bits.
 */
struct MemoryEfficientLCP
{
	/** @brief The exponent part of the LCP representation (base 2). */
	uint8_t exp_of_2 = 0;
	/** @brief The multiple (mantissa) part of the LCP representation. */
	uint8_t multiple = 0;

	/** @brief Default comparison operator. Compares `exp_of_2` then `multiple`. */
	auto operator<=>(const MemoryEfficientLCP&) const = default;

	/**
	 * @brief Calculates the approximate actual LCP value represented by this struct.
	 * @return unsigned The calculated LCP value (multiple * 2^exp_of_2).
	 */
	unsigned value() const noexcept
	{
		// Calculate the value by left-shifting 1 by exp_of_2 and multiplying by multiple.
		return (1u << exp_of_2) * multiple;
	}

	/**
	 * @brief Prints the LCP representation (e.g., "2^3 * 5") to an output stream.
	 * @param os The output stream to write to.
	 * @return std::ostream& Reference to the output stream after writing.
	 */
	std::ostream& print(std::ostream& os) const noexcept;
};

/**
 * @struct GeometricRank
 * @brief Represents a rank used in Zip Trees, typically generated randomly.
 * Combines a geometric and a uniform random component for tie-breaking, crucial for the Zip Tree property.
 */
struct GeometricRank
{
	/** @brief Rank component typically following a geometric distribution (primary rank). */
	uint8_t geometric_rank;
	/** @brief Rank component typically following a uniform distribution (used for tie-breaking). */
	uint8_t uniform_rank;

	/**
	 * @brief Default comparison operator. Compares `geometric_rank` first, then `uniform_rank` for tie-breaking.
	 * Higher ranks are considered "greater".
	 */
	auto operator<=>(const GeometricRank&) const = default;

	/**
	 * @brief Generates a random rank using static random number generators.
	 * The geometric component helps achieve logarithmic expected depth, while the uniform component
	 * ensures uniqueness with high probability.
	 * @return GeometricRank A randomly generated rank.
	 * @note Uses static generators internally. This means successive calls produce different random numbers,
	 * but it is **not thread-safe** if called from multiple threads concurrently without external locking.
	 * Consider providing thread-local or externally managed generators for concurrent use.
	 */
	static GeometricRank get_random()
	{
		// Static generators ensure the sequence is consistent across calls within a single thread,
		// but also mean it's not thread-safe without external locking.
		static std::random_device rd; // Obtain a random seed from the hardware device
		static std::default_random_engine generator(rd()); // Seed the default random engine
		// static std::default_random_engine generator(1); // Use for deterministic testing (fixed seed)
		static std::geometric_distribution<uint8_t> g_dist(0.5); // p=0.5: ~50% chance of 0, ~25% of 1, etc.
		static std::uniform_int_distribution<uint8_t> u_dist(0, std::numeric_limits<uint8_t>::max()); // Uniform distribution over all possible uint8_t values

		// Generate and return a rank with both components
		return { g_dist(generator), u_dist(generator) };
	}
};

/**
 * @class ZipTrie
 * @brief Implements a Zip Trie, a data structure combining properties of Tries and Zip Trees (Treaps).
 *
 * A Zip Trie stores keys (represented by `BitString`) and associated ranks (`RANK_T`).
 * It maintains a tree structure based on both key prefixes (like a standard Trie) and
 * random ranks (like a Zip Tree/Treap). This combination aims to achieve balanced expected depth
 * (logarithmic in the number of elements) while efficiently handling shared prefixes.
 *
 * The structure relies on comparing keys based on their Longest Common Prefix (LCP) with ancestor paths
 * to navigate the trie efficiently. Rotations (zips) based on ranks ensure the heap property (Zip property)
 * is maintained, leading to the balanced structure on average.
 *
 * @tparam CHAR_T The underlying character type used in the keys (`BitString`). Typically `char` or `uint8_t`.
 * @tparam MEMORY_EFFICIENT If true, uses `MemoryEfficientLCP` for storing LCP values, trading precision for memory.
 * If false, uses `unsigned` for exact LCP storage.
 * @tparam RANK_T The type used for node ranks (default: `GeometricRank`). Must provide a static `get_random()`
 * method and support comparison operators (`<=>`).
 * @tparam CHAR_SIZE_BITS The number of significant bits per character in `CHAR_T`. Defaults to `sizeof(CHAR_T) * 8`.
 * Allows handling types where not all bits are used (e.g., 7-bit ASCII in an 8-bit char).
 */
template <typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T = GeometricRank, unsigned CHAR_SIZE_BITS = sizeof(CHAR_T) * 8>
class ZipTrie
{
public:
	/** @brief The type used to store Longest Common Prefix values. Conditional based on `MEMORY_EFFICIENT`. */
	using LCP_T = std::conditional_t<MEMORY_EFFICIENT, MemoryEfficientLCP, unsigned>;
	/** @brief The type used for keys, an instantiation of the `BitString` template. */
	using KEY_T = BitString<CHAR_T, CHAR_SIZE_BITS>;

	/**
	 * @brief Constructs an empty ZipTrie.
	 * @param max_size Hint for reserving memory. Pre-allocates internal storage (`_buckets`)
	 * to potentially avoid reallocations during insertions.
	 * @param max_lcp_length The maximum possible LCP length between any two keys expected to be stored.
	 * This is primarily used for the `MemoryEfficientLCP` calculation if enabled.
	 * Setting this accurately can improve the precision of approximate LCPs.
	 */
	ZipTrie(unsigned max_size, unsigned max_lcp_length);

	/**
	 * @brief Calculates the depth of a given key in the trie.
	 * The depth is the number of edges from the root to the node containing the key.
	 * @param key Pointer to the key to search for. The key must exist for a valid depth.
	 * @return int The depth of the key (root is at depth 0), or -1 if the key is not found in the trie.
	 */
	int get_depth(const KEY_T* key) const noexcept;

	/**
	 * @brief Calculates the height of the trie.
	 * The height is the number of edges on the longest path from the root to any leaf node.
	 * @return int The height of the trie, or -1 if the trie is empty.
	 */
	int height() const noexcept;

	/**
	 * @brief Calculates the average depth of all nodes currently in the trie.
	 * This provides a measure of the average search cost.
	 * @return double The average depth. Returns 0.0 for an empty trie.
	 */
	double get_average_height() const noexcept;

	/**
	 * @brief Returns the number of nodes (keys) currently stored in the trie.
	 * @return unsigned The number of nodes.
	 */
	unsigned size() const noexcept;

	/**
	 * @brief Checks if a key exists within the trie.
	 * @param key Pointer to the key to search for.
	 * @return bool True if a node with the exact key exists, false otherwise.
	 */
	bool contains(const KEY_T* key) const noexcept;

	/**
	 * @brief Finds the Longest Common Prefix (LCP) between the search key and the path taken in the trie.
	 * If the key is found, this returns the LCP between the search key and the found key (which is the key's full length).
	 * If the key is not found, it returns the LCP between the search key and the path leading to the
	 * point where the search terminated (i.e., the LCP with the closest existing keys).
	 * @param key Pointer to the key to search for.
	 * @return LCP_T The calculated LCP value, potentially approximate if `MEMORY_EFFICIENT` is true.
	 */
	LCP_T lcp(const KEY_T* key) const noexcept;

	/**
	 * @brief Inserts a key into the Zip Trie.
	 *
	 * A new node is created for the key, assigned a random rank using `RANK_T::get_random()`,
	 * and inserted into the trie. The insertion maintains both the Trie property (ordering based on keys)
	 * and the Zip property (heap ordering based on ranks) using rotations (zips).
	 *
	 * **Assumes keys are unique.** Inserting a key that already exists (or compares as equal)
	 * results in the insertion being skipped, and the existing tree structure remains unchanged.
	 * The provided `key` pointer is stored directly in the node.
	 * **The lifetime of the object pointed to by `key` must exceed the lifetime of the ZipTrie
	 * or at least until the key is removed (if removal is implemented).**
	 *
	 * @param key Pointer to the new key to insert. Must be non-null and point to a valid `KEY_T` object.
	 */
	void insert(const KEY_T* key) noexcept;

	/**
	 * @brief Removes a node with a given key from the zip tree. (Currently commented out / Not implemented)
	 * @param key The key of the node to remove.
	 * @return bool True if a node was found and removed, false otherwise.
	 * @note This functionality is declared but not implemented in the provided code.
	 */
	// bool remove(const KEY_T& key) noexcept;

	/**
	 * @brief Gets the rank of the root node.
	 * @return const RANK_T& A constant reference to the root node's rank.
	 * @warning Behavior is undefined if the trie is empty (`_root_index == NULLPTR`).
	 */
	const RANK_T& get_root_rank() const noexcept
	{
		// Assert or check for empty trie might be useful here in a debug build.
		return _buckets[_root_index].rank;
	}

	/**
	 * @brief Directly adds a pre-constructed bucket to the internal storage vector.
	 * @warning Intended for testing or specific initialization scenarios **only**.
	 * This method bypasses the standard insertion logic (`insert`), including rank generation
	 * and maintaining the Zip/Trie properties. Using it improperly can corrupt the trie structure.
	 * @param key Pointer to the key for the new bucket.
	 * @param rank The rank to assign to the new bucket.
	 * @param left Index of the left child (use `NULLPTR` for none).
	 * @param right Index of the right child (use `NULLPTR` for none).
	 * @param predecessor_lcp The ancestor LCP with the predecessor path for this bucket.
	 * @param successor_lcp The ancestor LCP with the successor path for this bucket.
	 */
	void set(const KEY_T* key, RANK_T rank, unsigned left, unsigned right, LCP_T predecessor_lcp, LCP_T successor_lcp) noexcept
	{
		_buckets.push_back({ key, rank, left, right, {predecessor_lcp, successor_lcp} });
	}

	/**
	 * @brief Directly sets the root index of the trie.
	 * @warning Intended for testing or specific initialization scenarios **only**.
	 * This method bypasses the standard insertion logic. Using it improperly can lead to
	 * an invalid trie state (e.g., pointing to an invalid index or violating trie properties).
	 * @param root_index The index in `_buckets` to set as the new root.
	 */
	void set_root_index(unsigned root_index) noexcept
	{
		_root_index = root_index;
	}

	/**
	 * @struct SearchResults
	 * @brief Holds the results returned by a search operation (`search` method).
	 */
	struct SearchResults
	{
		/** @brief True if a node with the exact key was found in the trie, false otherwise. */
		bool contains;
		/** @brief The Longest Common Prefix (LCP) calculated during the search. See `lcp()` method documentation for details on what this represents when `contains` is true or false. */
		LCP_T max_lcp;
		/** @brief The depth at which the key was found (number of edges from the root). Set to -1 if the key was not found (`contains` is false). */
		int depth;
	};

	/**
	 * @brief Performs a search for a key within the ZipTrie.
	 * Traverses the trie based on key comparisons (`k_compare`), updating ancestor LCPs along the path.
	 * @param key Pointer to the key to search for.
	 * @return SearchResults A struct containing the outcome: whether the key was found (`contains`),
	 * the calculated LCP (`max_lcp`), and the depth if found (`depth`).
	 */
	virtual SearchResults search(const KEY_T* key) const noexcept;

	/**
	 * @brief Generates a DOT language representation of the trie structure for visualization.
	 * Writes the DOT output to the specified file. This file can then be processed by tools
	 * like Graphviz (e.g., `dot -Tpng input.dot -o output.png`) to create a visual diagram of the trie.
	 * @param file_path The path (including filename) where the output DOT file should be written.
	 * If the file cannot be opened, an error message is printed to `std::cerr`.
	 */
	void to_dot(const std::string& file_path) const noexcept;

protected:
	/** @brief Index of the root node within the `_buckets` vector. Initialized to `NULLPTR` for an empty trie. */
	unsigned _root_index;
	/** @brief Stores the maximum possible LCP length provided during construction. Used for `MemoryEfficientLCP` scaling. */
	unsigned _max_lcp_length; // TODO: Review if this is actively used or if _log_max_size suffices.
	/** @brief Stores the approximate log base 2 of the `max_size` provided during construction. Used for `MemoryEfficientLCP` scaling. */
	uint8_t _log_max_size;

	/** @brief Represents a null child pointer index. Uses the maximum value of `unsigned`. */
	static constexpr unsigned NULLPTR = std::numeric_limits<unsigned>::max();

	/**
	 * @struct AncestorLCPs
	 * @brief Stores the Longest Common Prefix (LCP) values associated with the path from the root
	 * down to a node's *parent*. These represent the LCPs calculated so far along the paths
	 * corresponding to the inorder predecessor and successor branches relative to the current search/insertion path.
	 * These are crucial for the `k_compare` optimization.
	 */
	struct AncestorLCPs
	{
		/** @brief LCP accumulated along the path corresponding to the inorder predecessor branch. */
		LCP_T predecessor = {}; // Default initialize (0 for unsigned, {0,0} for MemoryEfficientLCP)
		/** @brief LCP accumulated along the path corresponding to the inorder successor branch. */
		LCP_T successor = {};   // Default initialize
	};

	/**
	 * @struct Bucket
	 * @brief Represents a single node within the ZipTrie.
	 * Contains the key, rank, child indices, and the ancestor LCPs *at the time of its insertion*.
	 */
	struct Bucket
	{
		/** @brief Pointer to the key associated with this node. The ZipTrie does not own the key data. */
		const KEY_T* key;
		/** @brief The rank associated with this node, used for maintaining the Zip Tree (heap) property. */
		RANK_T rank;
		/** @brief Index of the left child in the `_buckets` vector, or `NULLPTR` if no left child. */
		unsigned left = NULLPTR;
		/** @brief Index of the right child in the `_buckets` vector, or `NULLPTR` if no right child. */
		unsigned right = NULLPTR;
		/** @brief Stores the `AncestorLCPs` passed down from the parent during this node's insertion.
		 * Used by `k_compare` when comparing against this node `v`. */
		AncestorLCPs ancestor_lcps = {};
	};

	/** @brief A contiguous vector storing all the nodes (`Bucket`s) of the trie. Node relationships (parent/child) are managed via indices (`_root_index`, `left`, `right`). */
	std::vector<Bucket> _buckets;

	/**
	 * @brief Converts a standard LCP length (`size_t`) into the `MemoryEfficientLCP` format.
	 * This is an internal helper function used only when `MEMORY_EFFICIENT` is true.
	 * It approximates the LCP value to fit into the smaller `MemoryEfficientLCP` struct.
	 * @param num The exact LCP length (number of matching bits) to convert.
	 * @return MemoryEfficientLCP The approximate LCP value in the memory-efficient representation.
	 * @note The conversion logic depends on `_log_max_size` for scaling.
	 */
	inline MemoryEfficientLCP get_memory_efficient_lcp(size_t num) const noexcept;

	/**
	 * @brief Recursive helper function for the public `search` method.
	 * Implements the core traversal logic using a provided comparison function.
	 * @tparam CompareFunction A callable type (e.g., lambda, function pointer) that takes
	 * `(const KEY_T* key, const Bucket& v_node, const AncestorLCPs& ancestor_lcps)`
	 * and returns a `ComparisonResult`.
	 * @param key The key being searched for.
	 * @param compare_func The comparison function to use for navigating the trie.
	 * @return SearchResults The result of the search operation.
	 */
	template <typename CompareFunction>
	SearchResults search_recursive(const KEY_T* key, CompareFunction compare_func) const noexcept;

	/**
	 * @brief Recursive helper function for inserting a new node `x`.
	 *
	 * Traverses the trie based on key comparisons (`compare_func`) and LCPs, updating `ancestor_lcps`
	 * along the path. When the correct insertion position (a `NULLPTR` link) is found, it inserts the node `x`.
	 * During the backtracking phase (as the recursion unwinds), it performs rotations (zips)
	 * if the rank of the newly inserted node `x` (or the root of the subtree returned by the recursive call)
	 * is higher than or equal to the rank of the current node `v`. This maintains the Zip property (max-heap on ranks).
	 *
	 * @tparam CompareFunction A callable type for comparing keys, same signature as in `search_recursive`.
	 * @param x Pointer to the new bucket being inserted (already added to `_buckets`).
	 * @param x_index Index of the new bucket `x` in the `_buckets` vector.
	 * @param v_index Index of the current node (`v`) being visited in the recursion. Starts with `_root_index`.
	 * @param ancestor_lcps LCPs inherited from the path above `v_index`. These are updated based on the traversal direction.
	 * @param compare_func The comparison function to use.
	 * @return unsigned The index of the root of the modified subtree after potential insertion and rotations.
	 * This will be either `v_index` (if no rotation occurred at this level or insertion happened deeper)
	 * or `x_index` (if `x` was rotated up to become the new root of this subtree).
	 */
	template<typename CompareFunction>
	unsigned insert_recursive(Bucket* x, unsigned x_index, unsigned v_index, AncestorLCPs ancestor_lcps, CompareFunction compare_func) noexcept;

// Protected Member Variables and Helper Structures/Functions used internally

protected:
	/**
	 * @struct ComparisonResult
	 * @brief Holds the result of comparing a search/insert key `x` with a node `v`'s key,
	 * considering the context of ancestor LCPs.
	 */
	struct ComparisonResult
	{
		/** @brief The relative ordering of `x` compared to `v` (`std::strong_ordering::less`, `greater`, or `equal`). */
		std::strong_ordering comparison;
		/** @brief The Longest Common Prefix (LCP) calculated between `x` and `v` during this specific comparison,
		 * potentially starting the comparison from an offset determined by `k_compare_prefix_check`. */
		LCP_T lcp;
	};

	/**
	 * @brief Performs an initial check based on ancestor LCPs before full key comparison.
	 * This is a core optimization in Zip Trie navigation (`k_compare`). It compares the LCPs
	 * accumulated along the search path (`ancestor_lcps`) with the LCPs stored in the current node `v`
	 * (`v.ancestor_lcps`, which were the path LCPs when `v` was inserted).
	 * If these path LCPs differ significantly, the relative order of the search key `x` and node `v`
	 * can often be determined *without* comparing the actual key data (`x->key` vs `v.key`).
	 *
	 * @param v The current node (`Bucket`) being considered.
	 * @param ancestor_lcps The LCPs accumulated along the current search/insert path down to `v`'s parent.
	 * @param[out] out_start_lcp_val If a full key comparison is needed (returns `std::nullopt`), this parameter
	 * is populated with the length (in bits) from which the actual key comparison
	 * should start. This avoids re-comparing the prefix already known to match.
	 * @return std::optional<ComparisonResult> If the comparison can be decided based on LCPs alone, returns
	 * a `ComparisonResult` containing the determined ordering and the
	 * minimum of the differing LCPs. If a full key comparison is required,
	 * returns `std::nullopt`.
	 */
	std::optional<ComparisonResult> k_compare_prefix_check(
		const Bucket& v,
		const AncestorLCPs& ancestor_lcps,
		size_t& out_start_lcp_val) const noexcept;

	/**
	 * @brief Converts an exact LCP length (number of bits) to the `LCP_T` type.
	 * If `MEMORY_EFFICIENT` is true, this calls `get_memory_efficient_lcp` for approximation.
	 * Otherwise, it performs a simple cast to `unsigned`.
	 * @param val The exact LCP length (number of matching bits).
	 * @return LCP_T The LCP value in the appropriate storage format (`MemoryEfficientLCP` or `unsigned`).
	 */
	inline LCP_T convert_lcp(size_t val) const noexcept;

private:
	// Private helper methods used internally by public or protected methods.

	/**
	 * @brief Compares a key `x` with the key in a node `v`, utilizing ancestor LCPs for optimization.
	 *
	 * This is the primary comparison function used for navigating the Zip Trie during search and insertion.
	 * It first calls `k_compare_prefix_check` to see if the comparison result can be determined solely
	 * based on the LCP information stored along the ancestor paths. If not (`k_compare_prefix_check`
	 * returns `std::nullopt`), it proceeds to perform a direct comparison of the key data (`x->key` vs `v.key`),
	 * starting from the bit offset (`out_start_lcp_val`) determined by the matching prefix check.
	 *
	 * @param x Pointer to the key being searched for or inserted.
	 * @param v Constant reference to the current node (`Bucket`) being compared against.
	 * @param ancestor_lcps The LCPs computed along the path from the root down to `v`'s parent.
	 * @return ComparisonResult A struct containing the relative ordering (`comparison`) of `x` vs `v`
	 * and the calculated Longest Common Prefix (`lcp`) between them in the context of this comparison.
	 */
	inline ComparisonResult k_compare(const KEY_T* x, const Bucket& v, const AncestorLCPs& ancestor_lcps) const noexcept;

	/**
	 * @brief Recursive helper function to calculate the height of the subtree rooted at `node_index`.
	 * @param node_index Index of the root of the subtree in the `_buckets` vector.
	 * @return int The height of the subtree (longest path from `node_index` to a leaf in its subtree),
	 * or -1 if `node_index` is `NULLPTR` (representing an empty subtree).
	 */
	int height(unsigned node_index) const noexcept;

	/**
	 * @brief Recursive helper function to calculate the sum of the depths of all nodes within a subtree.
	 * Used by `get_average_height`.
	 * @param node_index Index of the root of the subtree.
	 * @param depth The depth of the `node_index` node itself (relative to the overall trie root).
	 * @return uint64_t The sum of the depths of all nodes in the subtree rooted at `node_index`.
	 * Each node's depth contributes to the sum. Returns 0 for an empty subtree (`NULLPTR`).
	 */
	uint64_t get_total_depth(unsigned node_index, uint64_t depth) const noexcept;

	/**
	 * @brief Recursive helper function to generate the DOT language representation for a subtree.
	 * Called by the public `to_dot` method. It defines the current node and its edges to its children
	 * in the DOT format and then recursively calls itself for the children.
	 * @param os The output file stream (`std::ofstream`) to write the DOT commands to.
	 * @param node_index Index of the root of the subtree to output.
	 */
	void to_dot_recursive(std::ofstream& os, unsigned node_index) const noexcept;
};

//-----------------------------------------------------------------------------
// ZipTrie Method Implementations
// (Definitions are usually placed in the header for templates)
//-----------------------------------------------------------------------------

/**
 * @brief Constructor implementation. Initializes members and reserves bucket capacity.
 */
template<typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::ZipTrie(unsigned max_size, unsigned max_lcp_length)
	: _root_index(NULLPTR),
	  _max_lcp_length(max_lcp_length),
	  // Calculate log2(max_size) safely, handling max_size = 0 case.
	  _log_max_size(max_size > 0 ? static_cast<uint8_t>(std::log2(static_cast<double>(max_size))) : 0)
{
	// Reserve space in the vector to potentially avoid reallocations if max_size is a good estimate.
	_buckets.reserve(max_size);
}

/**
 * @brief Implementation of `contains`. Calls `search` and returns only the boolean `contains` field from the result.
 */
template<typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
bool ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::contains(const KEY_T* key) const noexcept
{
	// Delegate the actual search logic to the search method.
	return search(key).contains;
}

/**
 * @brief Implementation of `lcp`. Calls `search` and returns only the `max_lcp` field from the result.
 */
template<typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
typename ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::LCP_T ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::lcp(const KEY_T* key) const noexcept
{
	// Delegate the actual search logic to the search method.
	return search(key).max_lcp;
}

/**
 * @brief Implementation of the LCP-based prefix check optimization.
 * Compares ancestor LCPs to potentially determine comparison outcome without full key check.
 */
template<typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
std::optional<typename ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::ComparisonResult> ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::k_compare_prefix_check(
	const Bucket& v,
	const AncestorLCPs& ancestor_lcps,
	size_t& out_start_lcp_val) const noexcept
{
	// Get the LCPs from the current search path.
	auto predecessor_lcp = ancestor_lcps.predecessor;
	auto successor_lcp = ancestor_lcps.successor;

	// Determine the maximum LCP encountered on the current search path.
	auto x_max_lcp = std::max(predecessor_lcp, successor_lcp);
	// Determine the corresponding LCP stored in node 'v' (from the path when 'v' was inserted).
	// If the search path's max LCP came from the predecessor side, use v's predecessor LCP, otherwise use v's successor LCP.
	auto corr_v_lcp = (predecessor_lcp > successor_lcp) ? v.ancestor_lcps.predecessor : v.ancestor_lcps.successor;

	// If the maximum LCP on the current path (x_max_lcp) differs from the
	// corresponding LCP recorded when node 'v' was inserted (corr_v_lcp)...
	if (x_max_lcp != corr_v_lcp)
	{
		// ...we can determine the order based on which LCP is larger and which path it came from.
		// The logic determines if x should be less or greater than v.
		// `(x_max_lcp > corr_v_lcp)` is true if the current path's LCP is longer.
		// `(predecessor_lcp > successor_lcp)` is true if the current path's max LCP came from the predecessor side.
		// If these booleans match, x is less; otherwise, x is greater.
		std::strong_ordering comparison_result =
			((x_max_lcp > corr_v_lcp) == (predecessor_lcp > successor_lcp)) ? std::strong_ordering::less : std::strong_ordering::greater;

		// The resulting LCP for this comparison is the minimum of the two differing LCPs.
		LCP_T result_lcp = std::min(x_max_lcp, corr_v_lcp);

		// Return the determined comparison result and LCP. No need for full key comparison.
		return ComparisonResult{ comparison_result, result_lcp };
	}

	// If the LCPs matched (x_max_lcp == corr_v_lcp), we cannot decide the order yet.
	// We need to perform a full key comparison, but we know the keys match up to `x_max_lcp`.
	// Set the starting offset for the key comparison.
	if constexpr (MEMORY_EFFICIENT)
	{
		// If using approximate LCPs, get the actual value.
		out_start_lcp_val = x_max_lcp.value();
	}
	else
	{
		// If using exact LCPs (unsigned), the value is directly usable.
		out_start_lcp_val = x_max_lcp;
	}

	// Indicate that a full key comparison should proceed by returning nullopt.
	return std::nullopt;
}

/**
 * @brief Implementation of the LCP value conversion helper.
 * Handles conditional conversion based on the `MEMORY_EFFICIENT` template parameter.
 */
template<typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
typename ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::LCP_T ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::convert_lcp(size_t val) const noexcept
{
	if constexpr (MEMORY_EFFICIENT) {
		// If memory efficiency is enabled, call the approximation function.
		return get_memory_efficient_lcp(val);
	} else {
		// Otherwise, just cast the exact value to the LCP type (unsigned).
		// Ensure the value fits; consider potential truncation if LCP_T were smaller than size_t.
		return static_cast<LCP_T>(val);
	}
}

/**
 * @brief Implementation of the core key comparison logic (`k_compare`).
 * Uses `k_compare_prefix_check` first, then performs `BitString::seq_k_compare` if needed.
 */
template<typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
typename ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::ComparisonResult ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::k_compare(const KEY_T* x, const Bucket& v, const AncestorLCPs& ancestor_lcps) const noexcept
{
	size_t start_lcp_val = 0; // Initialize starting bit offset for comparison.

	// Attempt to determine the result based on ancestor LCPs first.
	if (auto prefix_result = k_compare_prefix_check(v, ancestor_lcps, start_lcp_val))
	{
		// If k_compare_prefix_check returned a result, use it directly.
		return *prefix_result;
	}

	// If prefix check didn't yield a result, perform sequential key comparison
	// starting from the calculated bit offset 'start_lcp_val'.
	// Assumes KEY_T has a method `seq_k_compare(const KEY_T& other, size_t start_bit)`
	// that returns a pair {std::strong_ordering, size_t actual_lcp}.
	auto [comparison, actual_lcp_val] = x->seq_k_compare(*v.key, start_lcp_val);

	// Convert the resulting exact LCP length to the appropriate LCP_T format.
	return { comparison, convert_lcp(actual_lcp_val) };
}

/**
 * @brief Implementation of the public `search` method. Sets up and calls the recursive search helper.
 */
template<typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
typename ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::SearchResults ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::search(const KEY_T* key) const noexcept
{
	// Define a lambda function that encapsulates the call to the standard k_compare method.
	// This lambda matches the signature required by search_recursive.
	auto sequential_compare = [&](const KEY_T* k, const Bucket& v, const AncestorLCPs& ancestor_lcps)
		{
			// Call the member function k_compare for the actual comparison logic.
			return this->k_compare(k, v, ancestor_lcps);
		};

	// Call the recursive search helper, passing the key and the comparison lambda.
	return search_recursive(key, sequential_compare);
}

/**
 * @brief Implementation of the recursive search logic.
 * Traverses the trie iteratively based on the comparison function's results.
 */
template<typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
template<typename CompareFunction>
typename ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::SearchResults ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::search_recursive(const KEY_T* key, CompareFunction compare_func) const noexcept
{
	// Handle empty trie case.
	if (_root_index == NULLPTR)
	{
		// Return not found, default LCP, depth -1.
		return { false, {}, -1 };
	}

	unsigned v_index = _root_index; // Start traversal at the root.
	AncestorLCPs ancestor_lcps = {}; // Initialize ancestor LCPs to zero/default at the root.
	int depth = 0; // Initialize depth counter.

	// Iterate down the trie until a null pointer is reached or the key is found.
	while (v_index != NULLPTR)
	{
		// Get a reference to the current node.
		const Bucket& v_node = _buckets[v_index];

		// Perform the comparison using the provided comparison function.
		auto [ comparison, lcp ] = compare_func(key, v_node, ancestor_lcps);

		// Check if the key is an exact match.
		if (comparison == std::strong_ordering::equal)
		{
			// Key found. Return true, the final LCP, and the current depth.
			return { true, lcp, depth };
		}

		// If not equal, decide whether to go left or right.
		if (comparison == std::strong_ordering::less)
		{
			// Go left. Update the successor LCP for the path: it's the maximum of the
			// current successor LCP and the LCP computed in this step.
			ancestor_lcps.successor = std::max(ancestor_lcps.successor, lcp);
			// Move to the left child.
			v_index = v_node.left;
		}
		else // comparison == std::strong_ordering::greater
		{
			// Go right. Update the predecessor LCP similarly.
			ancestor_lcps.predecessor = std::max(ancestor_lcps.predecessor, lcp);
			// Move to the right child.
			v_index = v_node.right;
		}

		// Increment depth for the next level.
		++depth;
	}

	// Key not found after traversing the path.
	// Return false, the maximum LCP encountered along the path, and depth -1.
	return { false, std::max(ancestor_lcps.predecessor, ancestor_lcps.successor), -1 };
}

/**
 * @brief Implementation of the public `insert` method.
 * Creates a new node, assigns a random rank, adds it to the bucket list,
 * and calls the recursive insertion helper to place it correctly and perform rotations.
 */
template <typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
void ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::insert(const KEY_T* key) noexcept
{
	// Add a new bucket to the vector. Initialize with key and a new random rank.
	// Children and ancestor LCPs will be set during the recursive insertion.
	_buckets.push_back({ key, RANK_T::get_random() });
	// Get the index of the newly added bucket.
	unsigned new_node_index = _buckets.size() - 1;

	// Define the comparison function lambda, similar to the search method.
	auto sequential_compare = [&](const KEY_T* k, const Bucket& v, const AncestorLCPs& ancestor_lcps)
		{
			return this->k_compare(k, v, ancestor_lcps); // Use the standard k_compare
		};

	// Call the recursive insert helper to place the new node and update the tree structure.
	// Start recursion from the root (_root_index), passing the new node's pointer and index,
	// initial empty ancestor LCPs, and the comparison function.
	// Update the _root_index with the result, as the root might change due to rotations.
	_root_index = insert_recursive(&_buckets[new_node_index], new_node_index, _root_index, {}, sequential_compare);
}

/**
 * @brief Implementation of the recursive insertion logic with rotations (zips).
 * Traverses, inserts, and performs rotations on the way back up to maintain Zip property.
 */
template <typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
template<typename CompareFunction>
unsigned ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::insert_recursive(Bucket* x, unsigned x_index, unsigned v_index, AncestorLCPs ancestor_lcps, CompareFunction compare_func) noexcept
{
	// Base case: Reached a null pointer, meaning this is the insertion spot.
	if (v_index == NULLPTR)
	{
		// Set the ancestor LCPs for the new node 'x' based on the path taken.
		x->ancestor_lcps = ancestor_lcps;
		// Return the index of the new node 'x', indicating it's the root of this (now non-empty) subtree.
		return x_index;
	}

	// Get a reference to the current node 'v'.
	Bucket& v_node = _buckets[v_index];

	// Compare the new key 'x->key' with the current node's key 'v_node.key'.
	auto [ comparison, lcp ] = compare_func(x->key, v_node, ancestor_lcps);

	// If keys are equal, the key already exists. Do not insert duplicates.
	// Return the current node's index, leaving the subtree unchanged.
	if (comparison == std::strong_ordering::equal)
	{
		// Note: Could potentially update value associated with key here if storing values.
		// Since we only store keys, we just return.
		// Also, need to handle the bucket pushed back in `insert`. A robust implementation
		// might check for existence *before* adding to _buckets or remove the unused bucket here.
		// Current implementation leaves the unused bucket if key exists.
		return v_index;
	}

	// Navigate left or right based on the comparison result.
	if (comparison == std::strong_ordering::less)
	{
		// Go left. Recursively insert into the left subtree.
		// Update the ancestor LCPs for the recursive call: predecessor stays the same,
		// successor becomes the max of the current path's successor and the LCP just computed.
		unsigned subroot_index = insert_recursive(x, x_index, v_node.left, { ancestor_lcps.predecessor, std::max(ancestor_lcps.successor, lcp) }, compare_func);

		// After the recursive call returns, check if a rotation (zip) is needed.
		// Rotation occurs if the returned subtree root is the new node 'x' (meaning 'x' was inserted
		// somewhere below and might need to move up) AND 'x' has a rank greater than or equal to 'v'.
		if (subroot_index == x_index && x->rank >= v_node.rank)
		{
			// Perform a right rotation (promoting 'x' over 'v').
			v_node.left = x->right; // v adopts x's right child as its left child.
			x->right = v_index;     // x takes v as its right child.
			// Update LCPs stored *in the nodes* due to the structural change.
			// v's predecessor LCP is now the LCP computed between x and v.
			v_node.ancestor_lcps.predecessor = lcp;
			// x's ancestor LCPs become those passed into this level (ancestor_lcps).
			x->ancestor_lcps = ancestor_lcps;
			// Return x_index as the new root of this subtree.
			return x_index;
		}
		else
		{
			// No rotation needed, or rotation happened deeper in the recursion.
			// Simply update v's left child pointer to the returned subtree root.
			v_node.left = subroot_index;
		}
	}
	else // comparison == std::strong_ordering::greater
	{
		// Go right. Recursively insert into the right subtree.
		// Update ancestor LCPs: successor stays same, predecessor updates.
		unsigned subroot_index = insert_recursive(x, x_index, v_node.right, { std::max(ancestor_lcps.predecessor, lcp), ancestor_lcps.successor }, compare_func);

		// Check if a left rotation (promoting 'x' over 'v') is needed.
		if (subroot_index == x_index && x->rank >= v_node.rank)
		{
			// Perform a left rotation.
			v_node.right = x->left; // v adopts x's left child as its right child.
			x->left = v_index;      // x takes v as its left child.
			// Update LCPs stored in the nodes.
			v_node.ancestor_lcps.successor = lcp; // v's successor LCP is now the computed LCP.
			x->ancestor_lcps = ancestor_lcps;     // x inherits parent's LCPs.
			// Return x_index as the new root.
			return x_index;
		}
		else
		{
			// No rotation needed here. Update v's right child pointer.
			v_node.right = subroot_index;
		}
	}
	// If no rotation occurred at this level, 'v' remains the root of this subtree.
	return v_index;
}

/**
 * @brief Returns the number of nodes currently stored in the trie.
 */
template<typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
unsigned ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::size() const noexcept
{
	// Simply return the current size of the underlying vector storing the nodes.
	return _buckets.size();
	// Note: This assumes no gaps or marked-deleted nodes if removal were implemented differently.
	// If duplicate keys caused buckets to be added but not used in `insert_recursive`,
	// this size might be slightly inflated compared to unique keys.
}

/**
 * @brief Calculates the height of the trie by calling the recursive helper starting from the root.
 */
template <typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
int ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::height() const noexcept
{
	// Delegate to the recursive helper function starting at the root index.
	return height(_root_index);
}

/**
 * @brief Recursive implementation for calculating the height of a subtree.
 */
template <typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
int ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::height(unsigned node_index) const noexcept
{
	// Base case: An empty subtree (represented by NULLPTR) has a height of -1.
	if (node_index == NULLPTR)
	{
		return -1;
	}

	// Recursive step: The height of a non-empty subtree is 1 (for the current node)
	// plus the maximum of the heights of its left and right subtrees.
	return std::max(height(_buckets[node_index].left), height(_buckets[node_index].right)) + 1;
}

/**
 * @brief Gets the depth of a specific key by performing a search.
 */
template <typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
int ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::get_depth(const KEY_T* key) const noexcept
{
	// Perform a standard search for the key and return the depth field from the results.
	// If the key is not found, search returns depth -1, which is correctly propagated here.
	return search(key).depth;
}

/**
 * @brief Calculates the average depth (average height) of all nodes in the trie.
 */
template <typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
double ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::get_average_height() const noexcept
{
	// Get the total number of nodes.
	unsigned current_size = size();
	// Handle the case of an empty trie to avoid division by zero.
	if (current_size == 0)
	{
		return 0.0; // Average depth of an empty trie is 0.
	}
	// Calculate the sum of depths of all nodes using the recursive helper, starting at the root (depth 0).
	uint64_t total_depth_sum = get_total_depth(_root_index, 0);
	// Compute the average by dividing the total depth sum by the number of nodes.
	return static_cast<double>(total_depth_sum) / current_size;
}

/**
 * @brief Recursive implementation for calculating the sum of depths of all nodes in a subtree.
 */
template <typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
uint64_t ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::get_total_depth(unsigned node_index, uint64_t depth) const noexcept
{
	// Base case: An empty subtree contributes 0 to the total depth sum.
	if (node_index == NULLPTR)
	{
		return 0;
	}

	// Recursive step: The total depth sum for the subtree rooted at node_index is:
	// - The depth of the current node (`depth`)
	// - Plus the total depth sum of the left subtree (nodes in the left subtree are one level deeper)
	// - Plus the total depth sum of the right subtree (nodes in the right subtree are one level deeper)
	return depth + get_total_depth(_buckets[node_index].left, depth + 1) + get_total_depth(_buckets[node_index].right, depth + 1);
}

/**
 * @brief Converts a size_t LCP value to an approximate MemoryEfficientLCP.
 * Implements the specific approximation logic based on `_log_max_size`.
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
 * @brief Generates the DOT file header and initiates the recursive DOT generation process.
 * Handles file opening and basic DOT structure.
 */
template <typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
void ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::to_dot(const std::string& file_path) const noexcept
{
	// Attempt to open the specified file path for writing.
	std::ofstream os(file_path);
	if (!os) {
		// If the file cannot be opened, print an error message to standard error.
		std::cerr << "Error: Could not open file " << file_path << " for writing." << std::endl;
		return; // Exit the function if file opening failed.
	}
	// Write the standard DOT graph header.
	os << "digraph G {\n";
	// Define default node attributes (shape=record allows complex labels).
	os << "  node [shape=record, fontname=\"Arial\"];\n";
	// Define default edge attributes.
	os << "  edge [fontname=\"Arial\"];\n";
	// Call the recursive helper function to write the nodes and edges, starting from the root.
	to_dot_recursive(os, _root_index);
	// Write the closing brace for the DOT graph definition.
	os << "}\n";
	// File stream `os` is automatically closed when it goes out of scope.
}

/**
 * @brief Recursive implementation for generating DOT output for each node and its children.
 */
template <typename CHAR_T, bool MEMORY_EFFICIENT, typename RANK_T, unsigned CHAR_SIZE_BITS>
void ZipTrie<CHAR_T, MEMORY_EFFICIENT, RANK_T, CHAR_SIZE_BITS>::to_dot_recursive(std::ofstream& os, unsigned node_index) const noexcept
{
	// Base case: If the node index represents a null pointer, do nothing.
	if (node_index == NULLPTR)
	{
		return;
	}

	// Get a constant reference to the current node's data.
	const Bucket& node = _buckets[node_index];

	// Define the DOT representation for the current node.
	// Use a unique ID "node_<index>" for each node.
	// The label uses DOT's record shape syntax: { field1 | field2 | ... }
	// Display the Key, Rank (geometric, uniform), and stored Ancestor LCPs (predecessor, successor).
	// Note: Assumes KEY_T has an overloaded operator<< for printing.
	// Note: Explicitly cast uint8_t ranks to unsigned for printing numeric values.
	os << "  node_" << node_index << " [label=\"{Idx: " << node_index
	   << "|Key: " << *node.key // Assumes KEY_T is printable
	   << "|Rank: (" << static_cast<unsigned>(node.rank.geometric_rank) << "," << static_cast<unsigned>(node.rank.uniform_rank) << ")"
	   << "|LCPs: (" << node.ancestor_lcps.predecessor << ", " << node.ancestor_lcps.successor << ")}\"];\n";

	// Process the left child.
	if (node.left != NULLPTR)
	{
		// Draw an edge from the current node to the left child node.
		// Label the edge " L" for clarity.
		os << "  node_" << node_index << " -> node_" << node.left << " [label=\" L\"];\n";
		// Recursively call the function for the left child.
		to_dot_recursive(os, node.left);
	}
	else
	{
		// Optional: Explicitly draw null pointers for visualization.
		// Create a small, invisible point node for the null left child.
		// os << "  null_left_" << node_index << " [shape=point, style=invis];\n";
		// Draw an edge to the null point.
		// os << "  node_" << node_index << " -> null_left_" << node_index << " [label=\" L\", style=dashed];\n";
	}

	// Process the right child.
	if (node.right != NULLPTR)
	{
		// Draw an edge from the current node to the right child node.
		// Label the edge " R".
		os << "  node_" << node_index << " -> node_" << node.right << " [label=\" R\"];\n";
		// Recursively call the function for the right child.
		to_dot_recursive(os, node.right);
	}
	else
	{
		// Optional: Explicitly draw null pointers for visualization.
		// os << "  null_right_" << node_index << " [shape=point, style=invis];\n";
		// os << "  node_" << node_index << " -> null_right_" << node_index << " [label=\" R\", style=dashed];\n";
	}
}


//-----------------------------------------------------------------------------
// Standalone Function Implementations (related to structs defined earlier)
//-----------------------------------------------------------------------------

/**
 * @brief Implementation of the print method for `MemoryEfficientLCP`.
 * Formats the output as "2^exp * mult".
 */
inline std::ostream& MemoryEfficientLCP::print(std::ostream& os) const noexcept
{
	// Cast uint8_t to unsigned for correct stream output as numbers.
	return os << "2^" << static_cast<unsigned>(exp_of_2) << "*" << static_cast<unsigned>(multiple);
}

/**
 * @brief Overloaded stream insertion operator (`<<`) for `MemoryEfficientLCP`.
 * Allows printing `MemoryEfficientLCP` objects directly to output streams (e.g., `std::cout << lcp;`).
 * @param os The output stream (e.g., `std::cout`, `std::ofstream`).
 * @param lcp The `MemoryEfficientLCP` object to print.
 * @return std::ostream& Reference to the output stream after printing.
 */
inline std::ostream& operator<<(std::ostream& os, const MemoryEfficientLCP& lcp)
{
	// Delegate the actual printing logic to the struct's print method.
	return lcp.print(os);
}

/**
 * @brief Overloaded stream insertion operator (`<<`) for `GeometricRank`.
 * Allows printing `GeometricRank` objects directly to output streams (e.g., `std::cout << rank;`).
 * @param os The output stream (e.g., `std::cout`, `std::ofstream`).
 * @param rank The `GeometricRank` object to print.
 * @return std::ostream& Reference to the output stream after printing.
 */
inline std::ostream& operator<<(std::ostream& os, const GeometricRank& rank)
{
	return os << "(" << static_cast<unsigned>(rank.geometric_rank) << ", "
		   << static_cast<unsigned>(rank.uniform_rank) << ")";
}
