/**
 * @file SkipTrie.hpp
 * @brief Definition of the SkipTrie class for efficient string storage and retrieval.
 *
 * The SkipTrie is a templated class designed to store strings in a manner that allows for fast insertion, deletion,
 * and lookup operations. It is particularly optimized for strings, utilizing a bit-level representation to minimize
 * comparison operations. The implementation uses a skip list data structure, where each node has links that skip over
 * multiple elements, allowing for logarithmic average time complexities for its primary operations.
 */

#pragma once

#include <cmath> // for log2
#include <cstddef> // for size_t
#include <memory> // for smart pointers
#include <random> // for random height generation
#include <vector> // for dynamic arrays

#include <unordered_map>

#include <fstream> // for file operations

#include "BitString.cuh" // for string representation

/**
 * @enum Direction
 * @brief Enumerates the possible traversal directions within the skip list.
 */
enum class Direction
{
	FORWARD, ///< Represents forward traversal direction.
	BACKWARD, ///< Represents backward traversal direction.
	INPLACE ///< Indicates no traversal, used for operations that do not require movement.
};

/**
 * @class SkipTrie
 * @brief Implements a skip list specifically optimized for string storage.
 *
 * @tparam CHAR_T The character type of the strings to be stored. Defaults to char.
 * @tparam CHAR_SIZE_BITS The size in bits of the character type. Defaults to the bit-size of CHAR_T.
 */
template<typename CHAR_T, unsigned CHAR_SIZE_BITS = sizeof(CHAR_T) * 8>
class SkipTrie
{
public:
	using KEY_T = BitString<CHAR_T, CHAR_SIZE_BITS>;

	/**
	 * @brief Constructs a new String Skip List object.
	 */
	SkipTrie();

	/**
	 * @brief Generates a random height for a new node based on a geometric distribution.
	 * @return size_t The randomly generated height.
	 */
	static size_t get_random_height();

	/**
	 * @brief Inserts a new key into the skip list.
	 * @param key A pointer to the BitString representing the key to insert.
	 */
	bool insert(const KEY_T* key) noexcept;

	/**
	 * @brief Inserts a new key into the skip list with a specified height.
	 * @param key A pointer to the BitString representing the key to insert.
	 * @param height The height at which to insert the key.
	 */
	bool insert(const KEY_T* key, size_t height) noexcept;

	bool contains(const KEY_T* key) const noexcept;

	bool remove(const KEY_T* key) noexcept;

	/**
	 * @brief Prints the contents of the skip list to the specified output stream.
	 * @param os The output stream to which the skip list should be printed.
	 * @return std::ostream& The same output stream for chaining.
	 */
	std::ostream& print(std::ostream& os) const noexcept;

	std::string to_string() const noexcept;

	size_t size() const noexcept { return m_size; }
	size_t height() const noexcept { return m_height; }

	size_t lcp(const KEY_T* key) const noexcept;
	size_t lcp_with_others(const KEY_T* key) const noexcept;

	// return all keys with the same prefix
	std::vector<const KEY_T*> suffix_search(const KEY_T* key) const noexcept;

	// return all keys within a certain range
	std::vector<const KEY_T*> range_search(const KEY_T* key1, const KEY_T* key2) const noexcept;

	std::unordered_map<size_t, size_t> get_lcp_group_sizes() const noexcept;

protected:
	size_t m_size; ///< Stores the number of elements in the skip list.
	size_t m_height; ///< Stores the current height of the skip list.

	/**
	 * @struct Node
	 * @brief Represents a node within the skip list.
	 */
	struct Node
	{
		const KEY_T* key; ///< The key stored in this node.
		std::shared_ptr<Node> next{nullptr}; ///< Pointer to the next node in the list.
		Node* prev = nullptr; ///< Pointer to the previous node in the list.
		std::shared_ptr<Node> down{nullptr}; ///< Pointer to the node below in the list.
		size_t lcp_next{0}; ///< The longest common prefix length with the next node.

		/**
		 * @brief Calculates the longest common prefix (LCP) with the specified direction's node.
		 * @param direction The direction in which to calculate the LCP.
		 * @return size_t The length of the LCP.
		 */
		size_t lcp(Direction direction) const noexcept;

		/**
		 * @brief Retrieves the next node in the specified direction.
		 * @param direction The direction in which to retrieve the next node.
		 * @return Node* A pointer to the next node.
		 */
		Node* next_node(Direction direction) const noexcept;
	};

	std::shared_ptr<Node> m_head = nullptr; ///< Pointer to the head node of the skip list.
	std::shared_ptr<Node> m_lower_head = nullptr; ///< Pointer to the head of the lowest level.
	std::shared_ptr<Node> m_lower_tail = nullptr; ///< Pointer to the tail of the lowest level.

	/**
	 * @brief Adds a new layer on top of the skip list to increase its height.
	 */
	void add_layer() noexcept;

	/**
	 * @brief Removes the top layer of the skip list to decrease its height.
	 */
	void remove_layer() noexcept;

	/**
	 * @brief Recursively inserts a key into the skip list at a specified height and direction.
	 * @param key A pointer to the BitString representing the key to insert.
	 * @param curr The current node being examined.
	 * @param height The height at which to insert the key.
	 * @param direction The current direction of insertion.
	 * @param current_height The current height in the recursive stack.
	 * @return std::shared_ptr<Node> A shared pointer to the newly inserted node, or nullptr if insertion was not needed.
	 */
	std::shared_ptr<Node> insert_recursive(const KEY_T* key, Node* curr, size_t height, Direction direction, size_t lcp, size_t current_height) noexcept;
	
	/**
	 * @brief Iterates through a layer of the skip list to find the appropriate insertion point for a key.
	 *
	 * This method iterates over the nodes in the specified direction until it finds the position where a new key should be
	 * inserted based on its lexicographical order. It handles the traversal logic, including decision-making when a key
	 * is equal to, less than, or greater than the current node's key during traversal.
	 *
	 * @param key A pointer to the BitString representing the key to insert.
	 * @param curr Reference to a pointer to the current node being examined. Updated as the method iterates.
	 * @param direction The initial direction for the iteration.
	 * @return Direction The direction in which the next iteration should proceed.
	 */
	Direction iter_layer(const KEY_T* key, Node*& curr, Direction direction, size_t& lcp) const noexcept;

	struct NodeLCP
	{
		Node* node;
		size_t lcp;
	};

	virtual NodeLCP find_first(const KEY_T* key, const bool require_level0 = false) const noexcept;
	
	struct EqualOrSuccessor
	{
		Node* node;
		bool is_equal;
		size_t lcp;
	};

	// NOTE: RESETS LCP
	virtual EqualOrSuccessor find_equal_or_successor(const KEY_T* key, const bool require_level0) const noexcept;

	struct ResultLCP
	{
		std::strong_ordering result;
		size_t lcp;
	};

	virtual ResultLCP compare(const KEY_T* key1, const KEY_T* key2, size_t lcp) const noexcept;

public:
	/**
	 * @class Iterator
	 * @brief An iterator for the SkipTrie class, allowing forward and backward traversal of the list.
	 *
	 * The Iterator class provides a means to traverse through the elements of the SkipTrie, either in forward
	 * or backward direction, depending on the instantiation. It encapsulates the logic for moving between nodes and
	 * accessing their keys, offering both pre-increment and post-increment operations.
	 */
	class Iterator
	{
	public:
		/**
		 * @brief Constructs a new Iterator object.
		 * @param node A pointer to the starting node for this iterator.
		 * @param direction The direction in which this iterator will move.
		 */
		Iterator(Node* node, Direction direction) noexcept;

		/**
		 * @brief Increments the iterator to the next node in its direction.
		 * @return Iterator& A reference to this iterator after moving to the next node.
		 */
		Iterator& operator++() noexcept;

		/**
		 * @brief Post-increment operator. Moves the iterator to the next node and returns an iterator to the previous node.
		 * @return Iterator An iterator to the node before this operation was performed.
		 */
		Iterator operator++(int) noexcept;

		/**
		 * @brief Dereference operator. Provides access to the key of the node currently pointed to by the iterator.
		 * @return const KEY_T& A reference to the key of the current node.
		 */
		const KEY_T& operator*() const noexcept;

		/**
		 * @brief Arrow operator. Provides pointer-like access to the key of the node currently pointed to by the iterator.
		 * @return const KEY_T* A pointer to the key of the current node.
		 */
		const KEY_T* operator->() const noexcept;

		/**
		 * @brief Equality comparison operator. Checks if two iterators are the same by comparing their current nodes.
		 * @param other The other iterator to compare with.
		 * @return bool True if both iterators point to the same node, false otherwise.
		 */
		bool operator==(const Iterator& other) const noexcept;

		/**
		 * @brief Inequality comparison operator. Checks if two iterators are not the same.
		 * @param other The other iterator to compare with.
		 * @return bool True if the iterators point to different nodes, false if they are the same.
		 */
		bool operator!=(const Iterator& other) const noexcept;

	private:
		Node* m_node; ///< Pointer to the current node in the traversal.
		Direction m_direction; ///< The direction of traversal.
	};

	/**
	 * @brief Retrieves an iterator to the beginning of the list.
	 * @return Iterator An iterator pointing to the first element of the list.
	 */
	Iterator begin() const noexcept;

	/**
	 * @brief Retrieves an iterator to the end of the list.
	 * @return Iterator An iterator pointing to the past-the-end element of the list.
	 */
	Iterator end() const noexcept;
};

/**
 * @brief Stream insertion operator for printing the SkipTrie.
 * @tparam CHAR_T The character type of the strings stored in the list.
 * @tparam CHAR_SIZE_BITS The size in bits of the character type.
 * @param os The output stream to which the skip list is to be printed.
 * @param ssl The SkipTrie object to print.
 * @return std::ostream& The output stream, for chaining.
 */
template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
std::ostream& operator<<(std::ostream& os, const SkipTrie<CHAR_T, CHAR_SIZE_BITS>& ssl);


template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
SkipTrie<CHAR_T, CHAR_SIZE_BITS>::Iterator::Iterator(Node* node, Direction direction) noexcept
	: m_node(node), m_direction(direction)
{
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
	++(*this);
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
	return !(*this == other);
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
typename SkipTrie<CHAR_T, CHAR_SIZE_BITS>::Iterator SkipTrie<CHAR_T, CHAR_SIZE_BITS>::begin() const noexcept
{
	return Iterator(m_lower_head->next.get(), Direction::FORWARD);
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
typename SkipTrie<CHAR_T, CHAR_SIZE_BITS>::Iterator SkipTrie<CHAR_T, CHAR_SIZE_BITS>::end() const noexcept
{
	return Iterator(m_lower_tail.get(), Direction::FORWARD);
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
size_t SkipTrie<CHAR_T, CHAR_SIZE_BITS>::Node::lcp(Direction direction) const noexcept
{
	if (direction == Direction::FORWARD)
	{
		return lcp_next;
	}

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

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
SkipTrie<CHAR_T, CHAR_SIZE_BITS>::SkipTrie()
	: m_size(0), m_height(0)
{
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
size_t SkipTrie<CHAR_T, CHAR_SIZE_BITS>::get_random_height()
{
	static std::random_device rd;
	static std::mt19937 gen(rd());
	// static std::mt19937 gen(1);
	static std::geometric_distribution<size_t> dist(0.5);
	
	return dist(gen);
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
void SkipTrie<CHAR_T, CHAR_SIZE_BITS>::add_layer() noexcept
{
	std::shared_ptr<Node> new_head = std::make_shared<Node>();
	std::shared_ptr<Node> new_tail = std::make_shared<Node>();

	new_head->down = m_head;
	if (m_head)
	{
		new_tail->down = m_head->next;
	}
	new_head->next = new_tail;
	new_tail->prev = new_head.get();

	m_head = new_head;
	++m_height;

	if (!m_lower_head)
	{
		m_lower_head = m_head;
		m_lower_tail = new_tail;
	}
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
void SkipTrie<CHAR_T, CHAR_SIZE_BITS>::remove_layer() noexcept
{
	m_head = m_head->down;
	--m_height;

	if (m_height == 0)
	{
		m_lower_head = nullptr;
		m_lower_tail = nullptr;
	}
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
auto SkipTrie<CHAR_T, CHAR_SIZE_BITS>::compare(const KEY_T* key1, const KEY_T* key2, size_t lcp) const noexcept -> ResultLCP
{
	auto [comparison, next_lcp] = key1->seq_k_compare(*key2, lcp);

	return { comparison, next_lcp };
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
Direction SkipTrie<CHAR_T, CHAR_SIZE_BITS>::iter_layer(const KEY_T* key, Node*& curr, Direction direction, size_t& lcp) const noexcept
{
	if (direction == Direction::INPLACE)
	{
		return Direction::INPLACE;
	}

	Node* next = curr->next_node(direction);

	while (next->key && curr->lcp(direction) >= lcp)
	{
		// if strictly greater, skip over
		if (curr->lcp(direction) > lcp)
		{
			curr = next;
			next = curr->next_node(direction);
			continue;
		}
		
		// if equal, must actually compare
		auto [comparison, next_lcp] = compare(key, next->key, lcp);
		lcp = next_lcp;

		curr = next;

		if (comparison == std::strong_ordering::equal)
		{
			return Direction::INPLACE;
		}
		else if (comparison == std::strong_ordering::less)
		{
			if (direction == Direction::FORWARD) // must return
			{
				return Direction::BACKWARD;
			}
		}
		else
		{
			if (direction == Direction::BACKWARD) // must return
			{
				return Direction::FORWARD;
			}
		}

		next = curr->next_node(direction);
	}

	return direction;
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
typename SkipTrie<CHAR_T, CHAR_SIZE_BITS>::NodeLCP SkipTrie<CHAR_T, CHAR_SIZE_BITS>::find_first(const KEY_T* key, const bool require_level0) const noexcept
{
	auto [node, is_equal, lcp] = find_equal_or_successor(key, require_level0);

	return { is_equal ? node : nullptr, lcp };
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
bool SkipTrie<CHAR_T, CHAR_SIZE_BITS>::contains(const KEY_T* key) const noexcept
{
	return find_first(key).node;
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
bool SkipTrie<CHAR_T, CHAR_SIZE_BITS>::insert(const KEY_T* key) noexcept
{
	return insert(key, get_random_height());
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
bool SkipTrie<CHAR_T, CHAR_SIZE_BITS>::insert(const KEY_T* key, size_t height) noexcept
{
	while (height + 1 >= m_height)
	{
		add_layer();
	}

	bool already_exists = insert_recursive(key, m_head.get(), height, Direction::FORWARD, 0, m_height - 1) == nullptr;

	if (already_exists)
	{
		return false;
	}

	++m_size;
	return true;
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
std::shared_ptr<typename SkipTrie<CHAR_T, CHAR_SIZE_BITS>::Node> SkipTrie<CHAR_T, CHAR_SIZE_BITS>::insert_recursive(const KEY_T* key, Node* curr, size_t height, Direction direction, size_t lcp, size_t current_height) noexcept
{
	Direction next_direction = iter_layer(key, curr, direction, lcp);

	if (next_direction == Direction::INPLACE) // key already exists
	{
		return nullptr;
	}

	// current node is the one to drop down from

	auto child = current_height == 0 ? nullptr : insert_recursive(key, curr->down.get(), height, next_direction, lcp, current_height - 1);

	if (current_height > height) // no need to insert at this level
	{
		return child;
	}

	if (!child && current_height > 0)
	{
		return nullptr;
	}

	std::shared_ptr<Node> new_node = std::make_shared<Node>();
	new_node->key = key;
	new_node->down = child;

	switch (next_direction)
	{
	case Direction::FORWARD:
		new_node->next = curr->next;
		new_node->prev = curr;
		curr->next = new_node;
		new_node->next->prev = new_node.get();

		new_node->lcp_next = curr->lcp_next;
		curr->lcp_next = lcp;
		break;
	case Direction::BACKWARD:
		new_node->next = curr->prev->next;
		new_node->prev = curr->prev;
		curr->prev = new_node.get();
		new_node->prev->next = new_node;

		new_node->lcp_next = lcp;
		break;
	case Direction::INPLACE:
		break; // should never happen
	}

	return new_node;
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
bool SkipTrie<CHAR_T, CHAR_SIZE_BITS>::remove(const KEY_T* key) noexcept
{
	Node* curr = find_first(key).node;

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
		return keys;
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
	if (!m_head)
	{
		return { nullptr, false, 0 };
	}

	Node* curr = m_head.get();
	Node* prev = nullptr;
	Direction direction = Direction::FORWARD;
	size_t lcp = 0;

	while (curr)
	{
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

		++i;
		curr = curr->next.get();
	}

	return lcp_group_sizes;
}

template<typename CHAR_T, unsigned CHAR_SIZE_BITS>
std::ostream& operator<<(std::ostream& os, const SkipTrie<CHAR_T, CHAR_SIZE_BITS>& ssl) {
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
	return "";
}