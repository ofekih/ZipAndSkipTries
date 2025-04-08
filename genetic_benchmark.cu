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

struct DataPair
{
	Gene gene_bs; // gene represented as BitStrings
	LongString gene_ls; // gene represented as a LongString
};

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

std::vector<size_t> generate_indices(size_t n)
{
	std::vector<size_t> indices(n);
	std::iota(indices.begin(), indices.end(), 0);
	return indices;
}

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

		save_construction_data("c-trie++", n, N, L, timer.elapsed_nanoseconds());

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

		save_construction_data("ZT", n, N, L, timer.elapsed_nanoseconds());

		timer.start();

		for (size_t _ = 0; _ < num_repetitions; ++_)
		{
			ZipTrie<Nucleotide, false, GeometricRank, 2> mi_zt(max_m, max_m);

			for (size_t j = 0; j < n; ++j)
			{
				mi_zt.insert(&data[j].gene_bs);
			}
		}

		save_construction_data("MI-ZT", n, N, L, timer.elapsed_nanoseconds());

		timer.start();

		for (size_t _ = 0; _ < num_repetitions; ++_)
		{
			ParallelZipTrie<Nucleotide, true, GeometricRank, 2> pzt(max_m, max_m);

			for (size_t j = 0; j < n; ++j)
			{
				pzt.insert(&data[j].gene_bs);
			}
		}

		save_construction_data("PZT", n, N, L, timer.elapsed_nanoseconds(), MIN_PAR_COMPARE_WORD_SIZE);

		timer.start();

		for (size_t _ = 0; _ < num_repetitions; ++_)
		{
			ParallelZipTrie<Nucleotide, false, GeometricRank, 2> mi_pzt(max_m, max_m);

			for (size_t j = 0; j < n; ++j)
			{
				mi_pzt.insert(&data[j].gene_bs);
			}
		}

		save_construction_data("MI-PZT", n, N, L, timer.elapsed_nanoseconds(), MIN_PAR_COMPARE_WORD_SIZE);
	}
}

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

		save_search_data("c-trie++", n, m, l, timer.elapsed_nanoseconds(), num_repetitions);

		timer.start();

		for (size_t _ = 0; _ < num_repetitions; ++_)
		{
			zt.contains(&data[i].gene_bs);
		}

		save_search_data("ZT", n, m, l, timer.elapsed_nanoseconds(), num_repetitions);

		timer.start();

		for (size_t _ = 0; _ < num_repetitions; ++_)
		{
			mi_zt.contains(&data[i].gene_bs);
		}

		save_search_data("MI-ZT", n, m, l, timer.elapsed_nanoseconds(), num_repetitions);

		timer.start();

		for (size_t _ = 0; _ < num_repetitions; ++_)
		{
			pzt.contains(&data[i].gene_bs);
		}

		save_search_data("PZT", n, m, l, timer.elapsed_nanoseconds(), num_repetitions, MIN_PAR_COMPARE_WORD_SIZE);

		timer.start();

		for (size_t _ = 0; _ < num_repetitions; ++_)
		{
			mi_pzt.contains(&data[i].gene_bs);
		}

		save_search_data("MI-PZT", n, m, l, timer.elapsed_nanoseconds(), num_repetitions, MIN_PAR_COMPARE_WORD_SIZE);
	}
}

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

