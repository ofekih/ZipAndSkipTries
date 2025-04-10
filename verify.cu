#include "src/BitString.cuh"
#include "src/SkipTrie.hpp"
#include "src/ParallelSkipTrie.cuh"
#include "src/synthetic.hpp"
#include "src/ZipTrie.hpp"
#include "src/ParallelZipTrie.cuh"
#include "src/data.hpp"

#include <stdio.h>
#include <cassert>
#include <vector>
#include <string>
#include <algorithm>

int main()
{
	size_t word_length = 10000;
	size_t num_words = 10000;
	size_t avg_lcp_length = 500;
	auto words = get_random_words(word_length, num_words, avg_lcp_length);

	std::vector<BitString<char>> bs(num_words);
	std::transform(words.begin(), words.end(), bs.begin(), [](const std::string& word) {
		return BitString<char>(word);
	});

	CPUTimer timer;

	// verify par_find_mismatch
	{
		timer.start("\nVerifying par_find_mismatch");

		size_t max_size_words = (word_length + BitString<char>::ALPHA - 1) / BitString<char>::ALPHA;
		uintmax_t *d_a = copy_to_device(bs[0].data(), bs[0].size());
		uintmax_t *d_largeblock = alloc_large_block_to_device(max_size_words);
		size_t max_copied = 0;

		for (size_t i = 20; i < bs.size(); i++)
		{
			size_t expected = bs[0].seq_k_compare(bs[i], 0, bs[0].size()).lcp;
			size_t result = bs[0].par_k_compare(bs[i], 0, bs[0].size(), d_a, d_largeblock, max_copied).lcp;

			if (expected != result)
			{
				printf("Verification failed for words %zu and %zu\n", 0, i);
				exit(EXIT_FAILURE);
			}
		}

		device_free(d_a);
		device_free(d_largeblock);

		timer.print();
	}

	// Verify SkipTrie
	{
		printf("\nVerifying SkipTrie\n");

		SkipTrie<char> trie;

		timer.start("\tInserting words");

		for (const auto& word : bs)
		{
			trie.insert(&word);
		}

		timer.print();

		timer.start("\tVerifying words");

		for (const auto& word : bs)
		{
			if (!trie.contains(&word)) {
				printf("Verification failed for word: %s\n", word.to_string().c_str());
				exit(EXIT_FAILURE);
			}
		}

		timer.print();
	}

	// Verify ParallelSkipTrie
	{
		printf("\nVerifying ParallelSkipTrie\n");

		ParallelSkipTrie<char> trie(word_length);

		timer.start("\tInserting words");

		for (const auto& word : bs)
		{
			trie.insert(&word);
		}

		timer.print();

		timer.start("\tVerifying words");

		for (const auto& word : bs)
		{
			if (!trie.contains(&word)) {
				printf("Verification failed for word: %s\n", word.to_string().c_str());
				exit(EXIT_FAILURE);
			}
		}

		timer.print();
	}

	// Verify ZipTrie
	{
		printf("\nVerifying ZipTrie\n");

		ZipTrie<char, false> trie(num_words, word_length);

		timer.start("\tInserting words");

		for (const auto& word : bs)
		{
			trie.insert(&word);
		}

		timer.print();

		timer.start("\tVerifying words");

		for (const auto& word : bs)
		{
			if (!trie.contains(&word)) {
				printf("Verification failed for word: %s\n", word.to_string().c_str());
				exit(EXIT_FAILURE);
			}
		}

		timer.print();
	}

	// Verify ParallelZipTrie
	{
		printf("\nVerifying ParallelZipTrie\n");

		ParallelZipTrie<char, false> trie(num_words, word_length);

		timer.start("\tInserting words");

		for (const auto& word : bs)
		{
			trie.insert(&word);
		}

		timer.print();

		timer.start("\tVerifying words");

		for (const auto& word : bs)
		{
			if (!trie.contains(&word)) {
				printf("Verification failed for word: %s\n", word.to_string().c_str());
				exit(EXIT_FAILURE);
			}
		}

		timer.print();
	}

	// Verify Memory-Efficient ZipTrie
	{
		printf("\nVerifying Memory-Efficient ZipTrie\n");

		ZipTrie<char, true> trie(num_words, word_length);

		timer.start("\tInserting words");

		for (const auto& word : bs)
		{
			trie.insert(&word);
		}

		timer.print();

		timer.start("\tVerifying words");

		for (const auto& word : bs)
		{
			if (!trie.contains(&word)) {
				printf("Verification failed for word: %s\n", word.to_string().c_str());
				exit(EXIT_FAILURE);
			}
		}

		timer.print();
	}

	// Verify Memory-Efficient ParallelZipTrie
	{
		printf("\nVerifying Memory-Efficient ParallelZipTrie\n");

		ParallelZipTrie<char, true> trie(num_words, word_length);

		timer.start("\tInserting words");

		for (const auto& word : bs)
		{
			trie.insert(&word);
		}

		timer.print();

		timer.start("\tVerifying words");

		for (const auto& word : bs)
		{
			if (!trie.contains(&word)) {
				printf("Verification failed for word: %s\n", word.to_string().c_str());
				exit(EXIT_FAILURE);
			}
		}

		timer.print();
	}

	printf("\nAll tests passed!\n");

	return 0;
}
