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

void genetics::to_dna(const Gene& gene, const std::string& filename)
{
	std::ofstream file(filename);
	if (!file.is_open())
	{
		throw std::invalid_argument("File not found");
	}

	file << gene;
}

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

	for (const auto& entry : std::filesystem::directory_iterator(dirname))
	{
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

	for (const auto& entry : std::filesystem::directory_iterator(dirname))
	{
		if (get_file_type(entry) == FileType::DNA)
		{
			std::ifstream file(entry.path());
			if (!file.is_open())
			{
				throw std::invalid_argument("File not found");
			}

			std::string dna_string;
			std::getline(file, dna_string);
			dna_strings.emplace_back(dna_string);
		}
	}

	return dna_strings;
}

namespace
{
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
	// load names of all files in the directory, adding them to a vector
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

	// print all to binary file, add one streampos of space between each, and then
	// store in another vector the streampos of that streampos
	timer.start("Writing names to binary file");
	std::ofstream file(dirname + combined_filename + DNA_BIN_COMBINED_EXTENSION, std::ios::binary);
	if (!file.is_open())
	{
		throw std::invalid_argument("File not found");
	}

	size_t num_names = filepaths.size();
	file.write(reinterpret_cast<const char*>(&num_names), sizeof(num_names));

	std::vector<std::streampos> streampos_positions;
	for (const auto& filepath : filepaths)
	{
		auto name = filepath.filename().string();
		name = name.substr(0, name.size() - DNA_BIN_EXTENSION.size());

		file.write(name.c_str(), name.size());
		file.write("\0", 1);
		streampos_positions.push_back(file.tellp());
		file.write(reinterpret_cast<const char*>(&streampos_positions.back()), sizeof(streampos_positions.back()));
	}

	timer.print();

	timer.start();
	for (unsigned i = 0; i < filepaths.size(); ++i)
	{
		auto pos = file.tellp();
		file.seekp(streampos_positions[i]);
		file.write(reinterpret_cast<const char*>(&pos), sizeof(pos));
		file.seekp(pos);

		auto gene = Gene::from_file(filepaths[i].string());
		gene.to_file(file);

		if (i % 1000 == 0)
		{
			printf("Wrote %u files, took %s\n", i, pretty_print(timer.elapsed_nanoseconds()).c_str());
		}
	}
}

GeneManager::GeneManager(const std::string& path) : m_path(path)
{
	m_file.open(path + DNA_BIN_COMBINED_EXTENSION, std::ios::binary);
	if (!m_file.is_open())
	{
		throw std::invalid_argument("File not found");
	}

	WallTimer timer;
	timer.start("Loading names and positions");

	size_t num_names;
	m_file.read(reinterpret_cast<char*>(&num_names), sizeof(num_names));

	m_geneinfos.reserve(num_names);
	// m_names.reserve(num_names);
	// m_positions.reserve(num_names);

	for (size_t i = 0; i < num_names; ++i)
	{
		std::string name;
		char c;
		while (m_file.read(&c, 1) && c != '\0')
		{
			name.push_back(c);
		}

		// m_names.emplace_back(name);

		std::streampos position;
		m_file.read(reinterpret_cast<char*>(&position), sizeof(position));
		// m_positions.push_back(position);
		m_geneinfos.push_back({ name, position, 0 });
	}

	timer.print();

	if (std::filesystem::exists(path + DNA_BIN_SIZES_EXTENSION))
	{
		load_sizes();
	}
	else
	{
		save_sizes();
	}

	timer.start("Sorting by size");
	std::sort(m_geneinfos.begin(), m_geneinfos.end(), [](const GeneInfo& a, const GeneInfo& b) { return a.size < b.size; });
	timer.print();
}

void GeneManager::save_sizes()
{
	WallTimer timer;
	timer.start("Loading & saving sizes");

	std::ofstream sizes_file(m_path + DNA_BIN_SIZES_EXTENSION, std::ios::binary);
	if (!sizes_file.is_open())
	{
		throw std::invalid_argument("File not found");
	}

	for (auto& info : m_geneinfos)
	{
		m_file.seekg(info.position);
		m_file.read(reinterpret_cast<char*>(&info.size), sizeof(info.size));
		sizes_file.write(reinterpret_cast<const char*>(&info.size), sizeof(info.size));
	}

	timer.print();
}

void GeneManager::load_sizes()
{
	WallTimer timer;
	timer.start("Loading sizes");

	std::ifstream sizes_file(m_path + DNA_BIN_SIZES_EXTENSION, std::ios::binary);
	if (!sizes_file.is_open())
	{
		throw std::invalid_argument("File not found");
	}

	for (auto& info : m_geneinfos)
	{
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
	for (const auto& info : m_geneinfos)
	{
		m_file.seekg(info.position);
		genes.emplace_back(Gene::from_file(m_file));
	}

	return genes;
}

std::vector<Gene> GeneManager::get_random_genes(size_t count) const
{
	static std::random_device rd;
	static std::mt19937 gen(rd());
	static std::uniform_int_distribution<size_t> dist(0, m_geneinfos.size() - 1);

	std::unordered_set<size_t> positions;

	while (positions.size() < count)
	{
		positions.insert(m_geneinfos[dist(gen)].position);
	}

	std::vector<Gene> genes;

	for (const auto& position : positions)
	{
		m_file.seekg(position);
		genes.emplace_back(Gene::from_file(m_file));
	}

	return genes;
}
