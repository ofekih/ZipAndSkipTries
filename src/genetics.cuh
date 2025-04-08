#pragma once

#include "BitString.cuh"
#include "nucleotide.hpp"

#include <filesystem>
#include <istream>
#include <string>
#include <vector>

typedef BitString<Nucleotide, 2> Gene;

namespace genetics
{
	static const std::string GENETIC_DIRECTORY = "genetic/";
	static const std::string VIRUS_DIRECTORY = GENETIC_DIRECTORY + "virus/";
	static const std::string ABC_HUMI_DIRECTORY = GENETIC_DIRECTORY + "ABC-HuMi/";
	static const std::string DNA_EXTENSION = ".dna";
	static const std::string DNA_BIN_EXTENSION = ".dna.bin";
	static const std::string DNA_BIN_COMBINED_EXTENSION = ".dna.bin.combined";
	static const std::string DNA_BIN_SIZES_EXTENSION = ".dna.bin.sizes";
	static const std::string GBK_EXTENSION = ".gbk";
	static const std::string FNA_EXTENSION = ".fna";

	enum class FileType
	{
		GBK,
		FNA,
		DNA,
		DNA_BIN,
		DNA_BIN_COMBINED,
		UNKNOWN
	};

	Gene from_gbk(const std::string& filename);
	Gene from_fna(const std::string& filename);
	Gene from_dna(const std::string& filename);
	Gene from_file(const std::string& filename);
	Gene from_file(const std::filesystem::path& path);
	Gene from_file(const std::filesystem::directory_entry& entry);
	void to_dna(const Gene& gene, const std::string& filename);
	std::vector<Gene> from_directory(const std::string& dirname);
	std::vector<Gene> from_directory(const std::string& dirname, FileType type);
	std::vector<std::string> get_dna_strings_from_directory(const std::string& dirname);
	FileType get_file_type(const std::string& filename);
	FileType get_file_type(const std::filesystem::path& path);
	FileType get_file_type(const std::filesystem::directory_entry& entry);
	Gene from_stream(std::istream& stream);
	std::vector<std::filesystem::path> get_virus_nucleotide_paths();
	void combine_dna_bin_files(const std::string& dirname, const std::string& combined_filename);

	struct GeneManager
	{
	public:
		GeneManager(const std::string& path);
		~GeneManager();

		std::vector<Gene> all_genes() const;
		std::vector<Gene> get_random_genes(size_t count) const;
	private:
		struct GeneInfo
		{
			std::string name;
			std::streampos position;
			size_t size;
		};

		std::string m_path;
		mutable std::ifstream m_file;
		std::vector<GeneInfo> m_geneinfos;

		void save_sizes();
		void load_sizes();
	};
};
