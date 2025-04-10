/**
 * @file genetics.cu
 * @brief Implementation of genetics utility functions.
 * @details This file implements the functions declared in genetics.cuh for loading,
 * saving, and manipulating genetic data from various file formats. It provides
 * parsing functionality for different genetic data file formats and utilities for
 * working with directories of genetic data.
 *
 * @see genetics.cuh
 * @see nucleotide.hpp
 */

#include "genetics.cuh"
#include "nucleotide.hpp"
#include "data.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <unordered_set>
#include <vector>

using namespace genetics;

/**
 * @brief Implementation of from_stream function that parses genetic sequences from a stream.
 * @details Reads the stream line by line, extracting valid nucleotides and adding them to 
 * the gene sequence. Invalid characters are skipped. Aborts if a FASTA header marker ('>') 
 * is encountered, as this indicates incorrect parsing logic.
 * 
 * @param stream Input stream containing genetic data.
 * @return Gene The parsed genetic sequence.
 */
Gene genetics::from_stream(std::istream& stream)
{
	Gene gene;

	std::string line;
	while (std::getline(stream, line))
	{
		for (char c : line)
		{
			if (is_valid_nucleotide(c))
			{
				gene.push_back(char_to_nucleotide[c]);
			}
			if (c == '>')
			{
				printf("NOT GOOD\n");
				exit(0);
			}
		}
	}

	return gene;
}

/**
 * @brief Implementation of from_gbk function for parsing GenBank format files.
 * @details GenBank files have a specific format with a header section followed by
 * a section that begins with "ORIGIN" followed by the nucleotide sequence. This
 * function skips the header section and then passes the file stream to from_stream
 * to extract the actual nucleotide sequence.
 * 
 * @param filename Path to the GenBank file.
 * @return Gene The parsed genetic sequence.
 * @throws std::invalid_argument if the file cannot be opened.
 */
Gene genetics::from_gbk(const std::string& filename)
{
	std::ifstream file(filename);
	if (!file.is_open())
	{
		throw std::invalid_argument("File not found");
	}

	// step 1: skip header
	// keep reading until the line just read is "ORIGIN"

	std::string line;
	while (std::getline(file, line))
	{
		if (line == "ORIGIN")
		{
			break;
		}
	}

	// step 2: read nucleotides
	// read line by line, character by character, discarding any non-nucleotide characters

	return from_stream(file);
}

/**
 * @brief Implementation of from_fna function for parsing FASTA nucleotide format files.
 * @details FASTA files start with a header line (beginning with '>') followed by the 
 * nucleotide sequence. This function skips the first line (header) and then passes 
 * the file stream to from_stream to extract the actual nucleotide sequence.
 * 
 * @param filename Path to the FASTA nucleotide file.
 * @return Gene The parsed genetic sequence.
 * @throws std::invalid_argument if the file cannot be opened.
 */
Gene genetics::from_fna(const std::string& filename)
{
	std::ifstream file(filename);
	if (!file.is_open())
	{
		throw std::invalid_argument("File not found");
	}

	// step 1: skip header (just one line)

	std::string line;
	std::getline(file, line);

	return from_stream(file);
}

/**
 * @brief Implementation of to_dna function for saving genes to DNA format files.
 * @details Writes the gene to a text file in DNA format, using the stream output
 * operator of the Gene class to convert the gene to a string representation.
 * 
 * @param gene The genetic sequence to save.
 * @param filename Path where the DNA file will be saved.
 * @throws std::invalid_argument if the file cannot be opened for writing.
 */
void genetics::to_dna(const Gene& gene, const std::string& filename)
{
	std::ofstream file(filename);
	if (!file.is_open())
	{
		throw std::invalid_argument("File not found");
	}

	file << gene;
}

/**
 * @brief Implementation of from_dna function for loading genes from DNA format files.
 * @details Opens a DNA format file and passes the file stream to from_stream to
 * extract the nucleotide sequence. DNA format files contain only nucleotide characters.
 * 
 * @param filename Path to the DNA format file.
 * @return Gene The parsed genetic sequence.
 * @throws std::invalid_argument if the file cannot be opened.
 */
Gene genetics::from_dna(const std::string& filename)
{
	std::ifstream file(filename);
	if (!file.is_open())
	{
		throw std::invalid_argument("File not found");
	}

	Gene gene;

	std::string line;
	while (std::getline(file, line))
	{
		for (char c : line)
		{
			gene.push_back(char_to_nucleotide[c]);
		}
	}

	return gene;
}

/**
 * @brief Implementation of from_directory function that loads genes from all files in a directory.
 * @details Iterates through all files in the specified directory and attempts to load
 * each file as a genetic sequence using automatic format detection. This version handles
 * all supported file types.
 * 
 * @param dirname Path to the directory containing genetic data files.
 * @return std::vector<Gene> Vector of parsed genetic sequences.
 */
std::vector<Gene> genetics::from_directory(const std::string& dirname)
{
	std::vector<Gene> genes;

	for (const auto& entry : std::filesystem::directory_iterator(dirname))
	{
		genes.emplace_back(from_file(entry));
	}

	return genes;
}

std::vector<Gene> genetics::from_directory(const std::string& dirname, FileType type)
{
	std::vector<Gene> genes;

	// Iterate through directory and filter by the requested file type
	for (const auto& entry : std::filesystem::directory_iterator(dirname))
	{
		// Only process files matching the specified type
		if (get_file_type(entry) == type)
		{
			genes.emplace_back(from_file(entry));
		}
	}

	return genes;
}

std::vector<std::string> genetics::get_dna_strings_from_directory(const std::string& dirname)
{
	std::vector<std::string> dna_strings;

	// Iterate through all files in the directory
	for (const auto& entry : std::filesystem::directory_iterator(dirname))
	{
		// Process only DNA format files
		if (get_file_type(entry) == FileType::DNA)
		{
			std::ifstream file(entry.path());
			if (!file.is_open())
			{
				throw std::invalid_argument("File not found");
			}

			// Read just the first line from each file
			// Note: This assumes DNA files contain a single line with the entire sequence
			std::string dna_string;
			std::getline(file, dna_string);
			dna_strings.emplace_back(dna_string);
		}
	}

	return dna_strings;
}

namespace
{
	/**
	 * @brief Checks if a string ends with a specified suffix.
	 * @details Provides a simple utility to determine if a string has a particular ending,
	 * which is used for file extension checking in the get_file_type functions.
	 * 
	 * @param str The string to check.
	 * @param suffix The suffix to look for at the end of the string.
	 * @return bool True if the string ends with the specified suffix, false otherwise.
	 */
	bool string_ends_with(const std::string& str, const std::string& suffix)
	{
		return str.size() >= suffix.size() && str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
	}
}

FileType genetics::get_file_type(const std::string& filename)
{
	if (string_ends_with(filename, GBK_EXTENSION))
	{
		return FileType::GBK;
	}
	else if (string_ends_with(filename, FNA_EXTENSION))
	{
		return FileType::FNA;
	}
	else if (string_ends_with(filename, DNA_EXTENSION))
	{
		return FileType::DNA;
	}
	else if (string_ends_with(filename, DNA_BIN_EXTENSION))
	{
		return FileType::DNA_BIN;
	}
	else if (string_ends_with(filename, DNA_BIN_COMBINED_EXTENSION))
	{
		return FileType::DNA_BIN_COMBINED;
	}
	else
	{
		return FileType::UNKNOWN;
	}
}

Gene genetics::from_file(const std::string& filename)
{
	FileType type = get_file_type(filename);

	switch (type)
	{
	case FileType::GBK:
		return from_gbk(filename);
	case FileType::FNA:
		return from_fna(filename);
	case FileType::DNA:
		return from_dna(filename);
	case FileType::DNA_BIN:
		return Gene::from_file(filename);
	default:
		throw std::invalid_argument("Unknown file type");
	}
}

FileType genetics::get_file_type(const std::filesystem::path& path)
{
	return get_file_type(path.string());
}

FileType genetics::get_file_type(const std::filesystem::directory_entry& entry)
{
	return get_file_type(entry.path().string());
}

Gene genetics::from_file(const std::filesystem::path& path)
{
	return from_file(path.string());
}

Gene genetics::from_file(const std::filesystem::directory_entry& entry)
{
	return from_file(entry.path().string());
}

std::vector<std::filesystem::path> genetics::get_virus_nucleotide_paths()
{
	std::vector<std::filesystem::path> paths;

	for (const auto& entry : std::filesystem::directory_iterator(VIRUS_DIRECTORY))
	{
		paths.emplace_back(entry.path());
	}

	return paths;
}

void genetics::combine_dna_bin_files(const std::string& dirname, const std::string& combined_filename)
{
	// Phase 1: Collect paths of all DNA_BIN files in the directory
	WallTimer timer;
	timer.start("Adding filepaths to vector");
	std::vector<std::filesystem::path> filepaths;
	for (const auto& entry : std::filesystem::directory_iterator(dirname))
	{
		if (get_file_type(entry) == FileType::DNA_BIN)
		{
			filepaths.emplace_back(entry.path());
		}
	}
	timer.print();

	printf("Added %lu filepaths to vector\n", filepaths.size());

	// Phase 2: Create the combined file header with file count, names, and placeholders for positions
	timer.start("Writing names to binary file");
	std::ofstream file(dirname + combined_filename + DNA_BIN_COMBINED_EXTENSION, std::ios::binary);
	if (!file.is_open())
	{
		throw std::invalid_argument("File not found");
	}

	// Write the number of files first for easy parsing later
	size_t num_names = filepaths.size();
	file.write(reinterpret_cast<const char*>(&num_names), sizeof(num_names));

	// For each file, write its name and reserve space for its position in the file
	std::vector<std::streampos> streampos_positions;
	for (const auto& filepath : filepaths)
	{
		// Extract the base filename without the .dna.bin extension
		auto name = filepath.filename().string();
		name = name.substr(0, name.size() - DNA_BIN_EXTENSION.size());

		// Write the name with null terminator
		file.write(name.c_str(), name.size());
		file.write("\0", 1);
		
		// Store the current position to update later with the actual gene data position
		streampos_positions.push_back(file.tellp());
		
		// Write a placeholder for the position (will update later)
		file.write(reinterpret_cast<const char*>(&streampos_positions.back()), sizeof(streampos_positions.back()));
	}

	timer.print();

	// Phase 3: Write each gene to the file and update its position in the header
	timer.start();
	for (unsigned i = 0; i < filepaths.size(); ++i)
	{
		// Remember current position where gene data will start
		auto pos = file.tellp();
		
		// Go back to the position placeholder in the header and update it
		file.seekp(streampos_positions[i]);
		file.write(reinterpret_cast<const char*>(&pos), sizeof(pos));
		
		// Return to the end of the file to write the gene data
		file.seekp(pos);

		// Load and write the gene
		auto gene = Gene::from_file(filepaths[i].string());
		gene.to_file(file);

		// Progress reporting for long operations
		if (i % 1000 == 0)
		{
			printf("Wrote %u files, took %s\n", i, pretty_print(timer.elapsed_nanoseconds()).c_str());
		}
	}
}

GeneManager::GeneManager(const std::string& path) : m_path(path)
{
	// Open the combined binary file
	m_file.open(path + DNA_BIN_COMBINED_EXTENSION, std::ios::binary);
	if (!m_file.is_open())
	{
		throw std::invalid_argument("File not found");
	}

	WallTimer timer;
	timer.start("Loading names and positions");

	// Read the number of genes stored in the file
	size_t num_names;
	m_file.read(reinterpret_cast<char*>(&num_names), sizeof(num_names));

	// Reserve space for all gene metadata
	m_gene_infos.reserve(num_names);
	// Note: Previous implementation stored names and positions in separate vectors
	// m_names.reserve(num_names);
	// m_positions.reserve(num_names);

	// Parse the gene metadata header
	for (size_t i = 0; i < num_names; ++i)
	{
		// Read null-terminated gene name
		std::string name;
		char c;
		while (m_file.read(&c, 1) && c != '\0')
		{
			name.push_back(c);
		}

		// Read the position where this gene's data starts in the file
		std::streampos position;
		m_file.read(reinterpret_cast<char*>(&position), sizeof(position));
		
		// Store the metadata (size will be loaded or determined later)
		m_gene_infos.push_back({ name, position, 0 });
	}

	timer.print();

	// Load or create the size information for each gene
	if (std::filesystem::exists(path + DNA_BIN_SIZES_EXTENSION))
	{
		// If size cache already exists, load it
		load_sizes();
	}
	else
	{
		// Otherwise, calculate and save sizes
		save_sizes();
	}

	// Sort genes by their size for optimized access patterns
	timer.start("Sorting by size");
	std::sort(m_gene_infos.begin(), m_gene_infos.end(), [](const GeneInfo& a, const GeneInfo& b) { return a.size < b.size; });
	timer.print();
}

void GeneManager::save_sizes()
{
	WallTimer timer;
	timer.start("Loading & saving sizes");

	// Create a companion file to store size information for quick loading in future
	std::ofstream sizes_file(m_path + DNA_BIN_SIZES_EXTENSION, std::ios::binary);
	if (!sizes_file.is_open())
	{
		throw std::invalid_argument("File not found");
	}

	// For each gene, read its size from the combined file and save to the sizes file
	for (auto& info : m_gene_infos)
	{
		// Navigate to where this gene's data begins
		m_file.seekg(info.position);
		
		// Read the size value from the gene data header
		m_file.read(reinterpret_cast<char*>(&info.size), sizeof(info.size));
		
		// Write the size to the companion file for future quick loading
		sizes_file.write(reinterpret_cast<const char*>(&info.size), sizeof(info.size));
	}

	timer.print();
}

void GeneManager::load_sizes()
{
	WallTimer timer;
	timer.start("Loading sizes");

	// Open the pre-computed sizes file which contains size information for each gene
	std::ifstream sizes_file(m_path + DNA_BIN_SIZES_EXTENSION, std::ios::binary);
	if (!sizes_file.is_open())
	{
		throw std::invalid_argument("File not found");
	}

	// Read size for each gene in the same order as the gene infos
	for (auto& info : m_gene_infos)
	{
		// Load the pre-computed size directly into the gene info structure
		sizes_file.read(reinterpret_cast<char*>(&info.size), sizeof(info.size));
	}

	timer.print();
}

GeneManager::~GeneManager()
{
	m_file.close();
}

std::vector<Gene> GeneManager::all_genes() const
{
	std::vector<Gene> genes;
	// Iterate through all gene infos and load each gene from the combined file
	for (const auto& info : m_gene_infos)
	{
		// Position file pointer at the start of this gene's data
		m_file.seekg(info.position);
		
		// Load the gene from the current file position and add to result vector
		genes.emplace_back(Gene::from_file(m_file));
	}

	return genes;
}

std::vector<Gene> GeneManager::get_random_genes(size_t count) const
{
	// Initialize static random number generator for efficient repeated calls
	static std::random_device rd;
	static std::mt19937 gen(rd());
	static std::uniform_int_distribution<size_t> dist(0, m_gene_infos.size() - 1);

	// Use a set to ensure we don't select the same gene twice
	// Store positions rather than indices to directly access the file
	std::unordered_set<size_t> positions;

	// Select random genes until we have the requested count
	while (positions.size() < count)
	{
		positions.insert(m_gene_infos[dist(gen)].position);
	}

	std::vector<Gene> genes;

	// Load each selected gene from its file position
	for (const auto& position : positions)
	{
		// Position file pointer at the gene's data
		m_file.seekg(position);
		
		// Load the gene and add to result vector
		genes.emplace_back(Gene::from_file(m_file));
	}

	return genes;
}
