/**
 * @file SkipTrie.hpp
 * @brief Defines the SkipTrie class, a skip list based data structure optimized for string storage and retrieval.
 * @details This file contains the definition of the SkipTrie template class, which uses a skip list
 * structure combined with bit-level string representation (`BitString`) for efficient operations like
 * insertion, search, and removal. It aims for logarithmic average time complexity by utilizing
 * skip links and efficient LCP (Longest Common Prefix) calculations.
 * @see BitString
 */
#pragma once

#include <cmath> // for log2
#include <cstddef> // for size_t
#include <memory> // for std::shared_ptr
#include <random> // for random height generation
#include <vector> // for std::vector
#include <unordered_map> // for std::unordered_map
#include <fstream> // for std::ostream
#include <string> // for std::string
#include <sstream> // for std::ostringstream in to_string
#include <compare> // for std::strong_ordering

#include "BitString.cuh"

/**
 * @enum Direction
 * @brief Enumerates the possible traversal directions within the skip list structure.
 * @details Used to indicate the direction of movement (forward or backward) or state (inplace)
 * during operations like insertion or search.
 */
enum class Direction
{
	FORWARD, ///< Represents traversal towards the end of the list (successor nodes).
	BACKWARD, ///< Represents traversal towards the beginning of the list (predecessor nodes).
	INPLACE ///< Indicates no traversal needed, typically when the exact key is found.
};

/**
 * @class SkipTrie
 * @brief Implements a skip list data structure optimized for storing and retrieving `BitString` keys.
 * @details This class uses a multi-level linked list (skip list) structure to store keys.
 * Each node can have connections to the next/previous node on the same level and a connection
 * to a corresponding node on the level below. Random heights are assigned to nodes upon insertion
 * to maintain probabilistic balance, aiming for logarithmic average time complexity for search,
 * insertion, and deletion. The structure leverages Longest Common Prefix (LCP) information
 * stored between adjacent nodes to speed up comparisons during traversal.
 *
 * @tparam CHAR_T The underlying character type used in the keys (`BitString`).
 * @tparam CHAR_SIZE_BITS The number of significant bits per character in `CHAR_T`. Defaults to `sizeof(CHAR_T) * 8`.
 * @see BitString
 */
template<typename CHAR_T, unsigned CHAR_SIZE_BITS = sizeof(CHAR_T) * 8>
class SkipTrie
{
public:
	/** @brief Alias for the key type used in the SkipTrie, based on `BitString`. */
	using KEY_T = BitString<CHAR_T, CHAR_SIZE_BITS>;

	/**
	 * @brief Constructs an empty SkipTrie.
	 * @details Initializes the size and height to zero. Head and tail pointers are initially null.
	 * Layers are added dynamically as needed during insertion.
	 */
	SkipTrie();

	/**
	 * @brief Generates a random height for a new node based on a geometric distribution.
	 * @details The height determines how many levels the new node will participate in.
	 * Uses a static random number generator (geometric distribution with p=0.5).
	 * @return size_t The randomly generated height (number of levels above the base level).
	 * @warning Not thread-safe if called concurrently from multiple threads due to static generator usage.
	 */
	static size_t get_random_height();

	/**
	 * @brief Inserts a new key into the skip list with a randomly generated height.
	 * @details Calls `get_random_height()` to determine the height and then calls the height-specific insert method.
	 * @param key A pointer to the `BitString` key to insert. The SkipTrie stores this pointer directly.
	 * @return bool True if the key was inserted successfully, false if the key already exists.
	 * @warning Key Lifetime: The lifetime of the object pointed to by `key` must exceed the lifetime of the SkipTrie.
	 * @see insert(const KEY_T*, size_t)
	 * @see get_random_height()
	 */
	bool insert(const KEY_T* key) noexcept;

	/**
	 * @brief Inserts a new key into the skip list with a specified height.
	 * @details Adds necessary layers if the specified height exceeds the current maximum height.
	 * Then, calls the recursive insertion helper function to place the node at the correct position
	 * across all levels up to the specified height. Increments the size if insertion is successful.
	 * @param key A pointer to the `BitString` key to insert. The SkipTrie stores this pointer directly.
	 * @param height The height (number of levels above the base) for the new node.
	 * @return bool True if the key was inserted successfully, false if the key already exists.
	 * @warning Key Lifetime: The lifetime of the object pointed to by `key` must exceed the lifetime of the SkipTrie.
	 * @see insert_recursive
	 * @see add_layer
	 */
	bool insert(const KEY_T* key, size_t height) noexcept;

	/**
	 * @brief Checks if a specific key exists within the skip list.
	 * @details Uses `find_first` to locate the key.
	 * @param key A pointer to the `BitString` key to search for.
	 * @return bool True if the key is found in the skip list, false otherwise.
	 * @see find_first
	 */
	bool contains(const KEY_T* key) const noexcept;

	/**
	 * @brief Removes a key from the skip list.
	 * @details Finds the node containing the key using `find_first`. If found, removes the node
	 * and its corresponding nodes on all levels below it by adjusting pointers. Decrements the size.
	 * Removes empty top layers if necessary.
	 * @param key A pointer to the `BitString` key to remove.
	 * @return bool True if the key was found and removed, false otherwise.
	 * @see find_first
	 * @see remove_layer
	 */
	bool remove(const KEY_T* key) noexcept;

	/**
	 * @brief Prints a textual representation of the skip list structure, layer by layer, to an output stream.
	 * @details Traverses each layer from head to tail, printing node keys and the LCP values between adjacent nodes.
	 * Useful for debugging and visualizing the structure.
	 * @param os The output stream (e.g., `std::cout`) to write the representation to.
	 * @return std::ostream& A reference to the output stream `os` for chaining.
	 */
	std::ostream& print(std::ostream& os) const noexcept;

	/**
	 * @brief Generates a string representation of the skip list structure.
	 * @details Similar to `print`, but captures the output into a `std::string`.
	 * @return std::string A string containing the textual representation of the skip list.
	 * @see print
	 */
	std::string to_string() const noexcept;

	/**
	 * @brief Returns the number of keys currently stored in the skip list.
	 * @return size_t The total number of keys (nodes at the base level).
	 */
	size_t size() const noexcept { return m_size; }

	/**
	 * @brief Returns the current maximum height of the skip list.
	 * @details The height is the number of levels in the skip list (0-indexed).
	 * An empty list has height 0 (after first insertion, height becomes >= 1).
	 * @return size_t The maximum height.
	 */
	size_t height() const noexcept { return m_height; }

	/**
	 * @brief Calculates the Longest Common Prefix (LCP) between a given key and its position in the skip list search path.
	 * @details Performs a search (`find_first`) for the key and returns the LCP value computed during that search.
	 * If the key is found, this LCP is typically the full length of the key. If not found, it's the LCP
	 * with the path leading to the insertion point.
	 * @param key A pointer to the `BitString` key.
	 * @return size_t The calculated LCP length (in characters).
	 * @see find_first
	 */
	size_t lcp(const KEY_T* key) const noexcept;

	/**
	 * @brief Calculates the maximum Longest Common Prefix (LCP) between a given key and its immediate neighbors (predecessor and successor) at the base level.
	 * @details Finds the node for the key (or its successor if not present) using `find_first` at level 0.
	 * Returns the maximum of the LCP with the next node and the LCP with the previous node at the base level.
	 * If the key is not found, it returns the LCP computed during the search for its potential position.
	 * @param key A pointer to the `BitString` key.
	 * @return size_t The maximum LCP length with its neighbors (or the search path LCP if not found).
	 * @see find_first
	 * @see Node::lcp
	 */
	size_t lcp_with_others(const KEY_T* key) const noexcept;

	/**
	 * @brief Finds all keys in the skip list that have the given key as a prefix.
	 * @details Performs a search (`find_equal_or_successor`) to locate the first potential match.
	 * If the LCP at that point equals the length of the input `key`, it iterates forward
	 * through the base level as long as the LCP between adjacent nodes is greater than or equal
	 * to the input key's length, collecting all matching keys.
	 * @param key A pointer to the `BitString` representing the prefix to search for.
	 * @return std::vector<const KEY_T*> A vector containing pointers to all keys that have `key` as a prefix.
	 * Returns an empty vector if no such keys are found.
	 * @see find_equal_or_successor
	 */
	std::vector<const KEY_T*> suffix_search(const KEY_T* key) const noexcept;

	/**
	 * @brief Finds all keys within a specified lexicographical range [key1, key2).
	 * @details Locates the starting node (equal to or successor of `key1`) using `find_equal_or_successor`.
	 * Locates the ending node (equal to or successor of `key2`). Iterates forward from the start node
	 * up to (but not including) the end node at the base level, collecting pointers to the keys.
	 * @param key1 A pointer to the `BitString` representing the lower bound of the range (inclusive).
	 * @param key2 A pointer to the `BitString` representing the upper bound of the range (exclusive).
	 * @return std::vector<const KEY_T*> A vector containing pointers to all keys within the specified range.
	 * @see find_equal_or_successor
	 */
	std::vector<const KEY_T*> range_search(const KEY_T* key1, const KEY_T* key2) const noexcept;

	/**
	 * @brief Calculates the sizes of groups of consecutive nodes sharing the same LCP value at the base level.
	 * @details Iterates through the base level of the skip list, tracking changes in LCP values between adjacent nodes.
	 * Uses a stack-like approach to determine the extent (size) of consecutive nodes that share a common LCP value
	 * greater than or equal to the LCP values of surrounding nodes.
	 * @return std::unordered_map<size_t, size_t> A map where the key is the LCP length and the value is the maximum
	 * size of a consecutive group found sharing at least that LCP length with their neighbors.
	 * @note This can be useful for analyzing the structure and prefix sharing within the dataset.
	 */
	std::unordered_map<size_t, size_t> get_lcp_group_sizes() const noexcept;

protected:
	/** @brief Stores the total number of keys (nodes) in the base level of the skip list. */
	size_t m_size;
	/** @brief Stores the current maximum height (number of levels) of the skip list. */
	size_t m_height;

	/**
	 * @struct Node
	 * @brief Represents a single node within the SkipTrie structure.
	 * @details Each node stores a pointer to a key, pointers to its neighbors (next, previous, down),
	 * and the Longest Common Prefix (LCP) length with its immediate successor on the same level.
	 */
	struct Node
	{
		/** @brief Pointer to the key associated with this node. The SkipTrie does not own the key data. */
		const KEY_T* key = nullptr;
		/** @brief Shared pointer to the next node on the same level. */
		std::shared_ptr<Node> next{nullptr};
		/** @brief Raw pointer to the previous node on the same level. */
		Node* prev = nullptr;
		/** @brief Shared pointer to the corresponding node on the level below. Null for base level nodes. */
		std::shared_ptr<Node> down{nullptr};
		/** @brief The Longest Common Prefix (LCP) length between this node's key and the next node's key on the same level. */
		size_t lcp_next{0};

		/**
		 * @brief Calculates the longest common prefix (LCP) with the neighbor node in the specified direction.
		 * @param direction The direction (`FORWARD` or `BACKWARD`) relative to this node.
		 * @return size_t The LCP length with the neighbor in the given direction. Returns 0 if there is no neighbor in that direction (e.g., `lcp(BACKWARD)` for the first node).
		 */
		size_t lcp(Direction direction) const noexcept;

		/**
		 * @brief Retrieves the neighboring node in the specified direction on the same level.
		 * @param direction The direction (`FORWARD` or `BACKWARD`) to retrieve the neighbor from.
		 * @return Node* A pointer to the neighboring node, or `nullptr` if no neighbor exists in that direction (e.g., `next_node(FORWARD)` for the tail sentinel).
		 */
		Node* next_node(Direction direction) const noexcept;
	};

	/** @brief Shared pointer to the head sentinel node of the top-most layer. */
	std::shared_ptr<Node> m_head = nullptr;
	/** @brief Shared pointer to the head sentinel node of the base layer (level 0). */
	std::shared_ptr<Node> m_lower_head = nullptr;
	/** @brief Shared pointer to the tail sentinel node of the base layer (level 0). */
	std::shared_ptr<Node> m_lower_tail = nullptr;

	/**
	 * @brief Adds a new, empty layer on top of the current skip list structure.
	 * @details Creates new head and tail sentinel nodes for the new top layer. Links the new head's `down` pointer
	 * to the old head and the new tail's `down` pointer to the old tail. Updates `m_head` and increments `m_height`.
	 * Initializes `m_lower_head` and `m_lower_tail` if this is the first layer being added.
	 * @note Called by `insert` when a node's random height exceeds the current maximum height.
	 */
	void add_layer() noexcept;

	/**
	 * @brief Removes the top-most layer of the skip list if it is empty (contains only sentinels).
	 * @details Updates `m_head` to point to the head of the layer below and decrements `m_height`.
	 * Resets `m_lower_head` and `m_lower_tail` if the height becomes zero.
	 * @note Called by `remove` to potentially shrink the skip list height after deletions.
	 */
	void remove_layer() noexcept;

	/**
	 * @brief Recursively finds the correct position and inserts links for a new node at a specific level and below.
	 * @details Traverses the current level (`current_height`) using `iter_layer` to find the predecessor node (`curr`)
	 * for the `key`. Recursively calls itself for the level below (`current_height - 1`). If the insertion height
	 * matches the `current_height`, it creates the new node at this level, links it between `curr` and `curr->next`,
	 * sets its `down` pointer to the node returned by the recursive call (or null at base level), and updates LCP values.
	 * @param key A pointer to the `BitString` key being inserted.
	 * @param curr The starting node for traversal at the current level.
	 * @param height The target insertion height for the new node (determines up to which level it's inserted).
	 * @param direction The initial traversal direction inherited from the level above.
	 * @param lcp The LCP length calculated so far along the search path from the level above.
	 * @param current_height The current level being processed in the recursion (starts from `m_height - 1`).
	 * @return std::shared_ptr<Node> A shared pointer to the newly created node at this level, or `nullptr` if the key already existed
	 * or if `current_height` is above the target `height`. This pointer is used to link the `down` pointer from the level above.
	 * @see iter_layer
	 */
	std::shared_ptr<Node> insert_recursive(const KEY_T* key, Node* curr, size_t height, Direction direction, size_t lcp, size_t current_height) noexcept;

	/**
	 * @brief Iterates horizontally along a single layer to find the node preceding the insertion/search point for a given key.
	 * @details Starts from `curr` and moves in the specified `direction`. Compares the `key` with the keys of subsequent nodes,
	 * using LCP values for optimization. Updates `curr` to point to the immediate predecessor of the `key`'s position on this layer.
	 * Updates the `lcp` parameter with the LCP value computed during the final relevant comparison.
	 * @param key A pointer to the `BitString` key being searched for or inserted.
	 * @param[in,out] curr A reference to the pointer to the current node. Updated by the function to the predecessor node.
	 * @param direction The initial traversal direction (usually `FORWARD` or `BACKWARD` inherited from the level above).
	 * @param[in,out] lcp A reference to the LCP length calculated so far. Updated by the function.
	 * @return Direction The direction the search should proceed on the level below (`FORWARD` or `BACKWARD`), or `INPLACE` if the exact key was found.
	 * @see compare
	 */
	Direction iter_layer(const KEY_T* key, Node*& curr, Direction direction, size_t& lcp) const noexcept;

	/**
	 * @struct NodeLCP
	 * @brief Helper struct to return a Node pointer and an associated LCP value.
	 * @details Used primarily by `find_first`.
	 */
	struct NodeLCP
	{
		/** @brief Pointer to the found node (or null if not found). */
		Node* node;
		/** @brief The LCP value computed during the search leading to this node. */
		size_t lcp;
	};

	/**
	 * @brief Finds the node containing the exact key, returning the node and the search path LCP.
	 * @details Uses `find_equal_or_successor` to perform the search. If the key is found (`is_equal` is true),
	 * returns the node and LCP. Otherwise, returns `nullptr` for the node.
	 * @param key A pointer to the `BitString` key to find.
	 * @param require_level0 If true, ensures the returned node (if found) is from the base level (level 0). Defaults to false.
	 * @return NodeLCP A struct containing the pointer to the found node (or `nullptr`) and the LCP computed during the search.
	 * @see find_equal_or_successor
	 */
	virtual NodeLCP find_first(const KEY_T* key, const bool require_level0 = false) const noexcept;

	/**
	 * @struct EqualOrSuccessor
	 * @brief Helper struct to return a Node pointer, a flag indicating equality, and an LCP value.
	 * @details Used by `find_equal_or_successor`.
	 */
	struct EqualOrSuccessor
	{
		/** @brief Pointer to the node found (either exact match or the successor). */
		Node* node;
		/** @brief True if `node` points to an exact match for the search key, false if it points to the successor. */
		bool is_equal;
		/** @brief The LCP value computed during the search leading to this node. */
		size_t lcp;
	};

	/**
	 * @brief Finds the node containing the exact key or its immediate successor.
	 * @details Traverses the skip list from the top layer down, using `iter_layer` at each level to narrow down the position.
	 * Returns the node found (either the exact match or the first node lexicographically greater than the key)
	 * and indicates whether an exact match was found. Also returns the LCP computed along the search path.
	 * @param key A pointer to the `BitString` key to search for.
	 * @param require_level0 If true, ensures the search descends fully to level 0 before returning the node. Defaults to false.
	 * @return EqualOrSuccessor A struct containing the found node (match or successor), an equality flag, and the search path LCP.
	 * @note This is a core search function used by `contains`, `insert`, `remove`, etc.
	 * @see iter_layer
	 */
	virtual EqualOrSuccessor find_equal_or_successor(const KEY_T* key, const bool require_level0) const noexcept;

	/**
	 * @struct ResultLCP
	 * @brief Helper struct to return a comparison result (ordering) and an LCP value.
	 * @details Used by `compare`.
	 */
	struct ResultLCP
	{
		/** @brief The result of the comparison (`std::strong_ordering::less`, `equal`, or `greater`). */
		std::strong_ordering result;
		/** @brief The Longest Common Prefix length calculated during the comparison. */
		size_t lcp;
	};

	/**
	 * @brief Compares two keys starting from a known common prefix length.
	 * @details Performs a sequential comparison of `key1` and `key2`, assuming the first `lcp` characters are identical.
	 * Uses `BitString::seq_k_compare` for the actual comparison.
	 * @param key1 Pointer to the first `BitString` key.
	 * @param key2 Pointer to the second `BitString` key.
	 * @param lcp The length of the known common prefix (in characters) from which to start the comparison.
	 * @return ResultLCP A struct containing the comparison result (`std::strong_ordering`) and the total LCP length found.
	 * @note This is the default comparison implementation. It can be overridden by derived classes (like `ParallelSkipTrie`).
	 * @see BitString::seq_k_compare
	 */
	virtual ResultLCP compare(const KEY_T* key1, const KEY_T* key2, size_t lcp) const noexcept;

public:
	/**
	 * @class Iterator
	 * @brief Provides forward/backward iteration over the keys in the base level of the SkipTrie.
	 * @details Allows traversal through the nodes at level 0 using standard iterator operations
	 * (`++`, `*`, `->`, `==`, `!=`). Can be constructed to move `FORWARD` or `BACKWARD`.
	 */
	class Iterator
	{
	public:
		/**
		 * @brief Constructs an iterator.
		 * @param node A pointer to the starting `Node` for the iterator (usually a node at level 0).
		 * @param direction The direction (`FORWARD` or `BACKWARD`) the iterator will move upon incrementing.
		 */
		Iterator(Node* node, Direction direction) noexcept;

		/**
		 * @brief Pre-increments the iterator to the next node in its designated direction.
		 * @return Iterator& A reference to the incremented iterator.
		 */
		Iterator& operator++() noexcept;

		/**
		 * @brief Post-increments the iterator.
		 * @return Iterator A copy of the iterator *before* it was incremented.
		 * @note Prefer pre-increment (`++it`) for potentially better performance.
		 */
		Iterator operator++(int) noexcept;

		/**
		 * @brief Dereferences the iterator to access the key of the current node.
		 * @return const KEY_T& A constant reference to the `BitString` key at the iterator's current position.
		 */
		const KEY_T& operator*() const noexcept;

		/**
		 * @brief Dereferences the iterator to access the key pointer of the current node.
		 * @return const KEY_T* A constant pointer to the `BitString` key at the iterator's current position.
		 */
		const KEY_T* operator->() const noexcept;

		/**
		 * @brief Compares this iterator with another for equality.
		 * @param other The iterator to compare against.
		 * @return bool True if both iterators point to the same `Node`, false otherwise.
		 */
		bool operator==(const Iterator& other) const noexcept;

		/**
		 * @brief Compares this iterator with another for inequality.
		 * @param other The iterator to compare against.
		 * @return bool True if the iterators point to different `Node`s, false otherwise.
		 */
		bool operator!=(const Iterator& other) const noexcept;

	private:
		/** @brief Pointer to the current node the iterator is positioned at. */
		Node* m_node;
		/** @brief The direction the iterator moves upon incrementing. */
		Direction m_direction;
	};

	/**
	 * @brief Returns an iterator pointing to the first key in the skip list (at the base level).
	 * @return Iterator An iterator positioned at the first element after the head sentinel.
	 */
	Iterator begin() const noexcept;

	/**
	 * @brief Returns an iterator pointing past the last key in the skip list (the tail sentinel at the base level).
	 * @return Iterator An iterator positioned at the tail sentinel node.
	 */
	Iterator end() const noexcept;
};

/**
 * @brief Overload of the stream insertion operator (`<<`) for `SkipTrie`.
 * @details Allows printing a `SkipTrie` object directly to an output stream using `os << skipTrie;`.
 * Delegates the actual printing logic to the `SkipTrie::print` method.
 * @tparam CHAR_T The character type of the keys.
 * @tparam CHAR_SIZE_BITS The number of significant bits per character.
 * @param os The output stream (e.g., `std::cout`).
 * @param ssl The `SkipTrie` object to print.
 * @return std::ostream& A reference to the output stream `os`.
 * @relates SkipTrie
 */
template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
std::ostream& operator<<(std::ostream& os, const SkipTrie<CHAR_T, CHAR_SIZE_BITS>& ssl);


// ========================================================================== //
//                           Implementation Details                           //
// ========================================================================== //

// --- Iterator Implementation ---

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
SkipTrie<CHAR_T, CHAR_SIZE_BITS>::Iterator::Iterator(Node* node, Direction direction) noexcept
	: m_node(node), m_direction(direction)
{
	// Constructor body intentionally empty
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
typename SkipTrie<CHAR_T, CHAR_SIZE_BITS>::Iterator& SkipTrie<CHAR_T, CHAR_SIZE_BITS>::Iterator::operator++() noexcept
{
	m_node = m_node->next_node(m_direction);
	return *this;
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
typename SkipTrie<CHAR_T, CHAR_SIZE_BITS>::Iterator SkipTrie<CHAR_T, CHAR_SIZE_BITS>::Iterator::operator++(int) noexcept
{
	Iterator temp = *this;
	++(*this); // Use pre-increment
	return temp;
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
const BitString<CHAR_T, CHAR_SIZE_BITS>& SkipTrie<CHAR_T, CHAR_SIZE_BITS>::Iterator::operator*() const noexcept
{
	return *m_node->key;
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
const BitString<CHAR_T, CHAR_SIZE_BITS>* SkipTrie<CHAR_T, CHAR_SIZE_BITS>::Iterator::operator->() const noexcept
{
	return m_node->key;
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
bool SkipTrie<CHAR_T, CHAR_SIZE_BITS>::Iterator::operator==(const Iterator& other) const noexcept
{
	return m_node == other.m_node;
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
bool SkipTrie<CHAR_T, CHAR_SIZE_BITS>::Iterator::operator!=(const Iterator& other) const noexcept
{
	return !(*this == other); // Delegate to operator==
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
typename SkipTrie<CHAR_T, CHAR_SIZE_BITS>::Iterator SkipTrie<CHAR_T, CHAR_SIZE_BITS>::begin() const noexcept
{
	// Returns iterator starting from the first actual node after the head sentinel at the base level.
	return Iterator(m_lower_head->next.get(), Direction::FORWARD);
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
typename SkipTrie<CHAR_T, CHAR_SIZE_BITS>::Iterator SkipTrie<CHAR_T, CHAR_SIZE_BITS>::end() const noexcept
{
	// Returns iterator pointing to the tail sentinel at the base level.
	return Iterator(m_lower_tail.get(), Direction::FORWARD);
}

// --- Node Implementation ---

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
size_t SkipTrie<CHAR_T, CHAR_SIZE_BITS>::Node::lcp(Direction direction) const noexcept
{
	if (direction == Direction::FORWARD)
	{
		// LCP with the next node is stored directly.
		return lcp_next;
	}

	// LCP with the previous node is stored in the previous node's lcp_next.
	return prev ? prev->lcp_next : 0;
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
typename SkipTrie<CHAR_T, CHAR_SIZE_BITS>::Node* SkipTrie<CHAR_T, CHAR_SIZE_BITS>::Node::next_node(Direction direction) const noexcept
{
	if (direction == Direction::FORWARD)
	{
		return next.get();
	}

	return prev;
}

// --- SkipTrie Implementation ---

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
SkipTrie<CHAR_T, CHAR_SIZE_BITS>::SkipTrie()
	: m_size(0), m_height(0)
{
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
size_t SkipTrie<CHAR_T, CHAR_SIZE_BITS>::get_random_height()
{
	// Use static variables for random number generation.
	// Note: Not thread-safe if multiple threads call this concurrently.
	static std::random_device rd; // Obtain a random seed from the hardware device
	static std::mt19937 gen(rd()); // Seed the default random engine
	// static std::mt19937 gen(1); // Use for deterministic testing (fixed seed)
	static std::geometric_distribution<size_t> dist(0.5); // p=0.5 gives 50% chance of 0, 25% of 1, etc.

	return dist(gen);
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
void SkipTrie<CHAR_T, CHAR_SIZE_BITS>::add_layer() noexcept
{
	// Create new sentinel nodes for the layer.
	std::shared_ptr<Node> new_head = std::make_shared<Node>();
	std::shared_ptr<Node> new_tail = std::make_shared<Node>();

	// Link the new head down to the old head (or null if no layers yet).
	new_head->down = m_head;
	// Link the new tail down to the old tail (or null if no layers yet).
	// Access old tail via old head's next pointer.
	if (m_head)
	{
		new_tail->down = m_head->next;
	}
	// Link new head and tail horizontally.
	new_head->next = new_tail;
	new_tail->prev = new_head.get();

	// Update the main head pointer to the new head.
	m_head = new_head;
	// Increment the height counter.
	++m_height;

	// If this is the very first layer added, set the lower head/tail pointers.
	if (!m_lower_head)
	{
		m_lower_head = m_head;
		m_lower_tail = new_tail;
	}
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
void SkipTrie<CHAR_T, CHAR_SIZE_BITS>::remove_layer() noexcept
{
	// Move the head pointer down to the next layer.
	m_head = m_head->down;
	// Decrement the height.
	--m_height;

	// If height becomes 0, reset lower head/tail pointers.
	if (m_height == 0)
	{
		m_lower_head = nullptr;
		m_lower_tail = nullptr;
	}
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
auto SkipTrie<CHAR_T, CHAR_SIZE_BITS>::compare(const KEY_T* key1, const KEY_T* key2, size_t lcp) const noexcept -> ResultLCP
{
	// Delegate comparison to the BitString's sequential compare method.
	auto [comparison, next_lcp] = key1->seq_k_compare(*key2, lcp);
	return { comparison, next_lcp };
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
Direction SkipTrie<CHAR_T, CHAR_SIZE_BITS>::iter_layer(const KEY_T* key, Node*& curr, Direction direction, size_t& lcp) const noexcept
{
	// If already found inplace, no need to iterate.
	if (direction == Direction::INPLACE)
	{
		return Direction::INPLACE;
	}

	// Get the next node in the current direction.
	Node* next = curr->next_node(direction);

	// Iterate while the next node exists and the LCP allows potential match/crossing.
	// Check `next->key` to ensure we don't compare against the tail sentinel.
	// `curr->lcp(direction)` gives LCP between `curr` and `next`.
	while (next->key && curr->lcp(direction) >= lcp)
	{
		// Optimization: If LCP between curr and next is strictly greater than current path LCP,
		// we know the target key cannot be between curr and next based on prefixes. Skip over next.
		if (curr->lcp(direction) > lcp)
		{
			curr = next;
			next = curr->next_node(direction);
			continue;
		}

		// If LCPs are equal, we must perform a full comparison starting from `lcp`.
		auto [comparison, next_lcp] = compare(key, next->key, lcp);
		// Update the path LCP with the result of the comparison.
		lcp = next_lcp;

		// Move to the next node.
		curr = next;

		// Check comparison result:
		if (comparison == std::strong_ordering::equal)
		{
			// Exact match found.
			return Direction::INPLACE;
		}
		else if (comparison == std::strong_ordering::less)
		{
			// Key is less than curr->key.
			if (direction == Direction::FORWARD)
			{
				// If moving forward, we've gone past the insertion point. Need to go backward on level below.
				return Direction::BACKWARD;
			}
			// If already moving backward, continue moving backward.
		}
		else // comparison == std::strong_ordering::greater
		{
			// Key is greater than curr->key.
			if (direction == Direction::BACKWARD)
			{
				// If moving backward, we've gone past the insertion point. Need to go forward on level below.
				return Direction::FORWARD;
			}
			// If already moving forward, continue moving forward.
		}

		// Get the next node for the next iteration.
		next = curr->next_node(direction);
	}

	// Loop finished: `next` is either the tail sentinel or `curr->lcp(direction) < lcp`.
	// The correct direction to proceed on the level below is the current `direction`.
	return direction;
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
typename SkipTrie<CHAR_T, CHAR_SIZE_BITS>::NodeLCP SkipTrie<CHAR_T, CHAR_SIZE_BITS>::find_first(const KEY_T* key, const bool require_level0) const noexcept
{
	// Use find_equal_or_successor to perform the search.
	auto [node, is_equal, lcp] = find_equal_or_successor(key, require_level0);

	// Return the node only if an exact match was found.
	return { is_equal ? node : nullptr, lcp };
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
bool SkipTrie<CHAR_T, CHAR_SIZE_BITS>::contains(const KEY_T* key) const noexcept
{
	// Check if find_first returns a non-null node.
	return find_first(key).node != nullptr;
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
bool SkipTrie<CHAR_T, CHAR_SIZE_BITS>::insert(const KEY_T* key) noexcept
{
	// Insert with a randomly generated height.
	return insert(key, get_random_height());
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
bool SkipTrie<CHAR_T, CHAR_SIZE_BITS>::insert(const KEY_T* key, size_t height) noexcept
{
	// Ensure the skip list has enough layers for the specified height.
	// Need height + 1 levels (0 to height). m_height is number of levels.
	while (height + 1 >= m_height)
	{
		add_layer();
	}

	// Call the recursive insertion helper. Starts from top head node (m_head).
	// Initial direction is FORWARD, initial LCP is 0, current height starts at m_height - 1.
	// insert_recursive returns nullptr if key already exists.
	bool already_exists = insert_recursive(key, m_head.get(), height, Direction::FORWARD, 0, m_height - 1) == nullptr;

	if (already_exists)
	{
		return false; // Key was not inserted because it already existed.
	}

	// Increment size only if insertion was successful.
	++m_size;
	return true;
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
std::shared_ptr<typename SkipTrie<CHAR_T, CHAR_SIZE_BITS>::Node> SkipTrie<CHAR_T, CHAR_SIZE_BITS>::insert_recursive(const KEY_T* key, Node* curr, size_t height, Direction direction, size_t lcp, size_t current_height) noexcept
{
	// Find the predecessor node (`curr`) on the current level.
	Direction next_direction = iter_layer(key, curr, direction, lcp);

	// If iter_layer found an exact match, return nullptr to signal duplication.
	if (next_direction == Direction::INPLACE)
	{
		return nullptr;
	}

	// Recursively insert on the level below, unless we are at the base level.
	std::shared_ptr<Node> child = nullptr;
	if (current_height > 0)
	{
		child = insert_recursive(key, curr->down.get(), height, next_direction, lcp, current_height - 1);

		// If recursive call returned nullptr, propagate nullptr up.
		if (!child)
		{
			return nullptr;
		}
	}


	// If the current level is above the target insertion height, just return the child pointer
	// (or nullptr if base level) without inserting at this level.
	if (current_height > height)
	{
		return child;
	}

	// --- Insert node at the current level ---

	// Create the new node.
	std::shared_ptr<Node> new_node = std::make_shared<Node>();
	new_node->key = key;
	new_node->down = child; // Link down to the node created/returned from the level below.

	// Link the new node horizontally based on the direction determined by iter_layer.
	// `curr` is the predecessor node found by iter_layer.
	switch (next_direction)
	{
	case Direction::FORWARD: // Insert new_node after curr
		new_node->next = curr->next;
		new_node->prev = curr;
		curr->next->prev = new_node.get(); // Update next node's prev pointer
		curr->next = new_node;

		// Update LCP values.
		new_node->lcp_next = curr->lcp_next; // New node inherits LCP curr had with its old next.
		curr->lcp_next = lcp; // LCP between curr and new_node is the final LCP from iter_layer.
		break;
	case Direction::BACKWARD: // Insert new_node before curr (i.e., after curr->prev)
		new_node->next = curr->prev->next; // Should be `curr` itself
		new_node->prev = curr->prev;
		curr->prev->next = new_node; // Update previous node's next pointer
		curr->prev = new_node.get();

		// Update LCP values.
		new_node->lcp_next = lcp; // LCP between new_node and curr is the final LCP from iter_layer.
		// new_node->prev->lcp_next remains unchanged (LCP between prev and new_node).
		break;
	case Direction::INPLACE:
		// This case should not be reached here due to the check after iter_layer.
		// If it were, it would indicate an error.
		break;
	}

	// Return the newly created node so the level above can link its `down` pointer.
	return new_node;
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
bool SkipTrie<CHAR_T, CHAR_SIZE_BITS>::remove(const KEY_T* key) noexcept
{
	Node* curr = find_first(key).node;

	// If key not found, return false.
	if (!curr)
	{
		return false;
	}

	// remove the node and all its children
	while (curr)
	{
		auto prev = curr->prev;
		auto next = curr->next;
		prev->lcp_next = std::min(prev->lcp_next, curr->lcp_next);
		curr = curr->down.get();
		
		prev->next = next;
		next->prev = prev;
	}

	// remove empty layers
	while (!m_head->down || !m_head->down->next->key)
	{
		remove_layer();

		if (m_height == 0)
		{
			break;
		}
	}

	--m_size;

	return true;
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
size_t SkipTrie<CHAR_T, CHAR_SIZE_BITS>::lcp(const KEY_T* key) const noexcept
{
	// Find the node (or insertion point) and return the LCP computed during search.
	return find_first(key).lcp;
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
size_t SkipTrie<CHAR_T, CHAR_SIZE_BITS>::lcp_with_others(const KEY_T* key) const noexcept
{
	auto [node, lcp] = find_first(key, true);

	if (!node)
	{
		return lcp;
	}

	return std::max(node->lcp(Direction::FORWARD), node->lcp(Direction::BACKWARD));
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
std::vector<const BitString<CHAR_T, CHAR_SIZE_BITS>*> SkipTrie<CHAR_T, CHAR_SIZE_BITS>::suffix_search(const KEY_T* key) const noexcept
{
	std::vector<const KEY_T*> keys;

	auto [curr, _, lcp] = find_equal_or_successor(key, true).node;

	if (lcp != key->size())
	{
		return keys; // No keys with this prefix found.
	}

	keys.push_back(curr->key);

	// iterate through the list until the edge lcp is less than the key length
	while (curr->next->key && curr->lcp_next >= lcp)
	{
		curr = curr->next.get();
		keys.push_back(curr->key);
	}

	return keys;
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
typename SkipTrie<CHAR_T, CHAR_SIZE_BITS>::EqualOrSuccessor SkipTrie<CHAR_T, CHAR_SIZE_BITS>::find_equal_or_successor(const KEY_T* key, const bool require_level0) const noexcept
{
	// Handle empty list case.
	if (!m_head)
	{
		return { nullptr, false, 0 };
	}

	Node* curr = m_head.get();
	Node* prev = nullptr;
	Direction direction = Direction::FORWARD;
	size_t lcp = 0;

	// Traverse down the levels.
	while (curr)
	{
		// Find predecessor on current level.
		Direction next_direction = iter_layer(key, curr, direction, lcp);

		if (next_direction == Direction::INPLACE)
		{
			while (require_level0 && curr->down)
			{
				curr = curr->down.get();
			}

			return { curr, true, lcp };
		}

		prev = curr;
		curr = curr->down.get();
		direction = next_direction;
	}

	return { direction == Direction::FORWARD ? prev->next.get() : prev, false, lcp };
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
std::vector<const BitString<CHAR_T, CHAR_SIZE_BITS>*> SkipTrie<CHAR_T, CHAR_SIZE_BITS>::range_search(const KEY_T* key1, const KEY_T* key2) const noexcept
{
	std::vector<const KEY_T*> keys;

	// Find the starting node (>= key1) at the base level.
	Node* curr = find_equal_or_successor(key1, true).node;

	if (!curr)
	{
		return keys;
	}

	auto [end, is_equal] = find_equal_or_successor(key2, true);

	if (is_equal)
	{
		end = end->next.get();
	}

	while (curr != end)
	{
		keys.push_back(curr->key);
		curr = curr->next.get();
	}

	return keys;
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
std::unordered_map<size_t, size_t> SkipTrie<CHAR_T, CHAR_SIZE_BITS>::get_lcp_group_sizes() const noexcept
{
	std::unordered_map<size_t, size_t> lcp_group_sizes;

	std::vector<size_t> lcps_we_are_greater_than;
	std::vector<size_t> index_first_encountered;

	Node* curr = m_lower_head->next.get();

	size_t i = 0;

	while (curr->key)
	{
		size_t lcp = curr->lcp_next;

		while (!lcps_we_are_greater_than.empty() && lcp < lcps_we_are_greater_than.back())
		{
			size_t popped_lcp = lcps_we_are_greater_than.back();
			lcps_we_are_greater_than.pop_back();

			size_t group_size = i - index_first_encountered.back();
			index_first_encountered.pop_back();

			lcp_group_sizes[popped_lcp] = std::max(lcp_group_sizes[popped_lcp], group_size);
		}

		if (lcps_we_are_greater_than.empty() || lcp > lcps_we_are_greater_than.back())
		{
			lcps_we_are_greater_than.push_back(lcp);
			index_first_encountered.push_back(i);
		}
		// If current_lcp is equal to top, do nothing (extend the current group).

		// Move to the next node and increment the gap index.
		curr = curr->next.get();
	}

	return lcp_group_sizes;
}

// --- Stream Operator Implementation ---

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
std::ostream& operator<<(std::ostream& os, const SkipTrie<CHAR_T, CHAR_SIZE_BITS>& ssl) {
	// Delegate to the print method.
	return ssl.print(os);
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
std::ostream& SkipTrie<CHAR_T, CHAR_SIZE_BITS>::print(std::ostream& os) const noexcept
{
	// print layer by layer
	auto left = m_head;
	size_t height = m_height - 1;
	while (left)
	{
		os << "Height: " << height << "|\t";
		os << "HEAD -0-> ";
		auto curr = left->next;
		while (curr->next)
		{
			// Print the key and the LCP to the next node.
			// Assumes KEY_T has an overloaded operator<<.
			os << *curr->key << " -" << curr->lcp_next << "-> ";
			curr = curr->next;
		}

		os << "TAIL\n";
		left = left->down;
		--height;
	}

	return os;
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
std::string SkipTrie<CHAR_T, CHAR_SIZE_BITS>::to_string() const noexcept
{
	// Use a string stream to capture the output of the print method.
	std::ostringstream oss;
	print(oss);
	return oss.str();
}
