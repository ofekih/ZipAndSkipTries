/**
 * @file genetic_benchmark.cu
 * @brief Benchmarking program for trie data structures using genetic data.
 *
 * @details This program benchmarks various trie implementations (SkipTrie, ZipTrie,
 * ParallelSkipTrie, ParallelZipTrie, and C-Trie++) on genetic data from the ABC-HuMi dataset.
 * It measures and compares construction and search performance across different
 * implementations, with a focus on genetic sequence data represented as nucleotides.
 * Results are saved to CSV files for further analysis and visualization.
 *
 * @see synthetic_benchmark.cu
 * @see src/genetics.cuh
 * @see src/BitString.cuh
 * @see src/SkipTrie.hpp
 * @see src/ZipTrie.hpp
 * @see src/ParallelSkipTrie.cuh
 * @see src/ParallelZipTrie.cuh
 */

#include "src/BitString.cuh"
#include "src/SkipTrie.hpp"
#include "src/ParallelSkipTrie.cuh"
#include "src/ZipTrie.hpp"
#include "src/ParallelZipTrie.cuh"
#include "src/genetics.cuh"
#include "src/data.hpp"

#include "ctriepp/ctriepp/CTriePP.hpp"
#include "ctriepp/ctriepp/LongString.hpp"

#include <chrono>
#include <algorithm>
#include <random>
#include <string>
#include <vector>
#include <locale.h>

#include <iostream>

using namespace genetics;
using namespace ctriepp;

/**
 * @struct DataPair
 * @brief A structure holding both BitString and LongString representations of genetic data.
 *
 * @details This structure is used to store the same genetic sequence in two different formats:
 * as a BitString (for our trie implementations) and as a LongString (for C-Trie++),
 * allowing for fair comparison between different data structures.
 */
struct DataPair
{
	Gene gene_bs; ///< Gene represented as BitStrings (for SkipTrie/ZipTrie variants)
	LongString gene_ls; ///< Gene represented as a LongString (for C-Trie++)
};

/**
 * @brief Loads genetic data from the ABC-HuMi dataset.
 *
 * @details Reads all genes from the ABC-HuMi dataset, converts them to both BitString and
 * LongString formats, and returns them as a vector of DataPair objects. The function
 * measures and reports the time taken to load and process the data.
 *
 * @return std::vector<DataPair> A vector containing all genes from the ABC-HuMi dataset
 * in both BitString and LongString formats.
 *
 * @see genetics::GeneManager
 */
std::vector<DataPair> load_abchumi_data()
{
	GeneManager gm(ABC_HUMI_DIRECTORY + "ABC-HuMi");

	CPUTimer timer;

	timer.start("Appending ABC-HuMi data");

	auto genes = gm.all_genes();

	static std::vector<std::string> words;
	std::vector<DataPair> data;
	words.reserve(genes.size());
	data.reserve(genes.size());
	for (const auto& gene : genes)
	{
		std::string gene_str;
		for (auto n : gene)
		{
			gene_str += nucleotide_to_char[static_cast<unsigned>(n)];
		}

		words.push_back(gene_str);
		data.push_back({ gene, LongString(&words.back()) });
	}

	timer.print();

	return data;
}

/**
 * @brief Shuffles the first n elements of the data vector.
 *
 * @details Uses the Fisher-Yates algorithm to randomize the order of the first n elements
 * in the data vector. This is used to create different insertion orders for benchmarking.
 *
 * @param[in,out] data The vector of DataPair objects to shuffle.
 * @param n The number of elements to shuffle.
 */
void shuffle_data(std::vector<DataPair>& data, size_t n)
{
	static std::mt19937 gen(std::random_device{}());

	// use fisher-yates to randomize the first n elements of the data
	for (size_t i = 0; i < n; ++i)
	{
		std::uniform_int_distribution<size_t> dist(i, data.size() - 1);
		size_t j = dist(gen);
		std::swap(data[i], data[j]);
	}
}

/**
 * @brief Generates a vector of indices from 0 to n-1.
 *
 * @details Creates a vector containing integers from 0 to n-1 in ascending order.
 * This is used for creating random access patterns for search benchmarks.
 *
 * @param n The number of indices to generate.
 * @return std::vector<size_t> A vector containing integers from 0 to n-1.
 */
std::vector<size_t> generate_indices(size_t n)
{
	std::vector<size_t> indices(n);
	std::iota(indices.begin(), indices.end(), 0);
	return indices;
}

/**
 * @brief Runs construction benchmarks on various trie implementations using genetic data.
 *
 * @details Measures the time taken to construct different trie data structures (C-Trie++,
 * ZipTrie, memory-intensive ZipTrie, ParallelZipTrie, memory-intensive ParallelZipTrie,
 * SkipTrie, and ParallelSkipTrie) with the same genetic data. For each implementation,
 * the benchmark is repeated multiple times and the results are saved to CSV files.
 *
 * @param[in,out] data The vector of DataPair objects containing the genetic data.
 * @param n The number of genes to use in the benchmark.
 * @param num_trials The number of trials to run (each with a different shuffle).
 * @param num_repetitions The number of repetitions for each trial (default: 100).
 *
 * @see save_construction_data
 */
void run_construction_benchmark(std::vector<DataPair>& data, size_t n, size_t num_trials, size_t num_repetitions = 100)
{
	static std::random_device rd;
	static std::mt19937 gen(rd());

	CPUTimer timer;

	auto indices = generate_indices(n);

	for (size_t i = 0; i < num_trials; ++i)
	{
		shuffle_data(data, n);
		std::shuffle(indices.begin(), indices.end(), gen);

		size_t max_m = 0;
		size_t N = 0;
		size_t L = 0;

		{
			SkipTrie<Nucleotide, 2> st;

			for (size_t j = 0; j < n; ++j)
			{
				st.insert(&data[j].gene_bs);
			}

			for (size_t j = 0; j < n; ++j)
			{
				const auto& gene = data[j].gene_bs;
				max_m = std::max(max_m, gene.size());
				N += gene.size();
				L += st.lcp_with_others(&gene);
			}
		}

		timer.start();

		for (size_t _ = 0; _ < num_repetitions; ++_)
		{
			CTriePP<bool, false> ctriepp;

			for (size_t j = 0; j < n; ++j)
			{
				ctriepp.insert(data[j].gene_ls, true);
			}
		}

		save_construction_data("c-trie++", n, N, L, timer.elapsed_nanoseconds(), 1, 0, true);

		timer.start();

		// next, test ZipTrie (sequential)
		for (size_t _ = 0; _ < num_repetitions; ++_)
		{
			ZipTrie<Nucleotide, true, GeometricRank, 2> zt(max_m, max_m);

			for (size_t j = 0; j < n; ++j)
			{
				zt.insert(&data[j].gene_bs);
			}
		}

		save_construction_data("ZT", n, N, L, timer.elapsed_nanoseconds(), 1, 0, true);

		timer.start();

		for (size_t _ = 0; _ < num_repetitions; ++_)
		{
			ZipTrie<Nucleotide, false, GeometricRank, 2> mi_zt(max_m, max_m);

			for (size_t j = 0; j < n; ++j)
			{
				mi_zt.insert(&data[j].gene_bs);
			}
		}

		save_construction_data("MI-ZT", n, N, L, timer.elapsed_nanoseconds(), 1, 0, true);

		timer.start();

		for (size_t _ = 0; _ < num_repetitions; ++_)
		{
			ParallelZipTrie<Nucleotide, true, GeometricRank, 2> pzt(max_m, max_m);

			for (size_t j = 0; j < n; ++j)
			{
				pzt.insert(&data[j].gene_bs);
			}
		}

		save_construction_data("PZT", n, N, L, timer.elapsed_nanoseconds(), 1, MIN_PAR_COMPARE_WORD_SIZE, true);

		timer.start();

		for (size_t _ = 0; _ < num_repetitions; ++_)
		{
			ParallelZipTrie<Nucleotide, false, GeometricRank, 2> mi_pzt(max_m, max_m);

			for (size_t j = 0; j < n; ++j)
			{
				mi_pzt.insert(&data[j].gene_bs);
			}
		}

		save_construction_data("MI-PZT", n, N, L, timer.elapsed_nanoseconds(), 1, MIN_PAR_COMPARE_WORD_SIZE, true);
	}
}

/**
 * @brief Runs search benchmarks on various trie implementations using genetic data.
 *
 * @details Measures the time taken to search for genes in different trie data structures
 * (C-Trie++, ZipTrie, memory-intensive ZipTrie, ParallelZipTrie, memory-intensive ParallelZipTrie)
 * with the same genetic data. For each implementation, the benchmark is repeated multiple times
 * with different search patterns, and the results are saved to CSV files.
 *
 * @param[in,out] data The vector of DataPair objects containing the genetic data.
 * @param n The number of genes to use in the benchmark.
 * @param num_trials The number of trials to run (each with a different search pattern).
 * @param num_repetitions The number of repetitions for each trial (default: 1000).
 *
 * @see save_search_data
 */
void run_contains_true_benchmark(std::vector<DataPair>& data, size_t n, size_t num_trials, size_t num_repetitions = 1000)
{
	shuffle_data(data, n);

	size_t max_m = 0;
	for (size_t j = 0; j < n; ++j)
	{
		const auto& gene = data[j].gene_bs;
		max_m = std::max(max_m, gene.size());
	}

	CTriePP<bool, false> ctriepp;
	ZipTrie<Nucleotide, true, GeometricRank, 2> zt(max_m, max_m);
	ZipTrie<Nucleotide, false, GeometricRank, 2> mi_zt(max_m, max_m);
	ParallelZipTrie<Nucleotide, true, GeometricRank, 2> pzt(max_m, max_m);
	ParallelZipTrie<Nucleotide, false, GeometricRank, 2> mi_pzt(max_m, max_m);

	for (size_t j = 0; j < n; ++j)
	{
		ctriepp.insert(data[j].gene_ls, true);
		zt.insert(&data[j].gene_bs);
		mi_zt.insert(&data[j].gene_bs);
		pzt.insert(&data[j].gene_bs);
		mi_pzt.insert(&data[j].gene_bs);
	}

	CPUTimer timer;

	for (size_t i = 0; i < num_trials; ++i)
	{
		auto m = data[i].gene_bs.size();
		auto l = m;

		timer.start();

		for (size_t _ = 0; _ < num_repetitions; ++_)
		{
			ctriepp.contains(data[i].gene_ls);
		}

		save_search_data("c-trie++", n, m, l, timer.elapsed_nanoseconds(), num_repetitions, 0, true);

		timer.start();

		for (size_t _ = 0; _ < num_repetitions; ++_)
		{
			zt.contains(&data[i].gene_bs);
		}

		save_search_data("ZT", n, m, l, timer.elapsed_nanoseconds(), num_repetitions, 0, true);

		timer.start();

		for (size_t _ = 0; _ < num_repetitions; ++_)
		{
			mi_zt.contains(&data[i].gene_bs);
		}

		save_search_data("MI-ZT", n, m, l, timer.elapsed_nanoseconds(), num_repetitions, 0, true);

		timer.start();

		for (size_t _ = 0; _ < num_repetitions; ++_)
		{
			pzt.contains(&data[i].gene_bs);
		}

		save_search_data("PZT", n, m, l, timer.elapsed_nanoseconds(), num_repetitions, MIN_PAR_COMPARE_WORD_SIZE, true);

		timer.start();

		for (size_t _ = 0; _ < num_repetitions; ++_)
		{
			mi_pzt.contains(&data[i].gene_bs);
		}

		save_search_data("MI-PZT", n, m, l, timer.elapsed_nanoseconds(), num_repetitions, MIN_PAR_COMPARE_WORD_SIZE, true);
	}
}

/**
 * @brief Main function that runs the genetic benchmarks.
 *
 * @details Loads genetic data from the ABC-HuMi dataset and runs construction and search
 * benchmarks on various trie implementations. The benchmark parameters (number of trials
 * and number of simulations) can be specified as command-line arguments.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line argument strings.
 * @return int Returns 0 on successful execution, 1 on invalid arguments.
 */
int main(int argc, char* argv[])
{
	setlocale(LC_ALL, "en_US.UTF-8");

	if (argc != 3)
	{
		std::cerr << "Usage: " << argv[0] << " <num_trials> <num_simulations>" << std::endl;
		return 1;
	}

	size_t num_trials = std::stoul(argv[1]);
	size_t num_simulations = std::stoul(argv[2]);

	CPUTimer timer;

	timer.start("Loading ABC-HuMi data");

	auto data = load_abchumi_data();

	timer.print();

	for (size_t n = 1; n <= data.size(); n *= 2)
	{
		for (size_t i = 0; i < num_simulations; ++i)
		{
			timer.start("Running benchmarks for n = " + std::to_string(n) + "(" + std::to_string(i + 1) + "/" + std::to_string(num_simulations) + ")");
			run_construction_benchmark(data, n, num_trials);
			run_contains_true_benchmark(data, n, num_trials);
			timer.print();
		}
	}

	return 0;
}

