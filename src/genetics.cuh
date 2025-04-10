/**
 * @file genetics.cuh
 * @brief Provides utilities for working with genetic data in various formats.
 * @details This file defines types and functions for loading, saving, and manipulating 
 * genetic data from various file formats such as GBK, FNA, and DNA. It uses the BitString 
 * template specialized for nucleotides to efficiently represent genetic sequences.
 *
 * @see BitString.cuh
 * @see nucleotide.hpp
 */

#pragma once

#include "BitString.cuh"
#include "nucleotide.hpp"

#include <filesystem>
#include <istream>
#include <string>
#include <vector>

/**
 * @typedef Gene
 * @brief A specialized BitString for representing genetic sequences.
 * @details Represents DNA sequences as bit strings with 2 bits per nucleotide,
 * allowing for compact representation and efficient operations on genetic data.
 */
typedef BitString<Nucleotide, 2> Gene;

/**
 * @namespace genetics
 * @brief Contains utilities for working with genetic data and file formats.
 * @details Provides constants, enums, and functions for loading, saving, and manipulating
 * genetic data from various file formats. The namespace encapsulates all functionality
 * related to genetic sequence processing.
 */
namespace genetics
{
	/** @brief Path to the genetic data directory */
	static const std::string GENETIC_DIRECTORY = "genetic/";
	
	/** @brief Path to the virus genome data directory */
	static const std::string VIRUS_DIRECTORY = GENETIC_DIRECTORY + "virus/";
	
	/** @brief Path to the ABC-HuMi dataset directory */
	static const std::string ABC_HUMI_DIRECTORY = GENETIC_DIRECTORY + "ABC-HuMi/";
	
	/** @brief File extension for plain text DNA files */
	static const std::string DNA_EXTENSION = ".dna";
	
	/** @brief File extension for binary DNA files */
	static const std::string DNA_BIN_EXTENSION = ".dna.bin";
	
	/** @brief File extension for combined binary DNA files */
	static const std::string DNA_BIN_COMBINED_EXTENSION = ".dna.bin.combined";
	
	/** @brief File extension for binary DNA file size information */
	static const std::string DNA_BIN_SIZES_EXTENSION = ".dna.bin.sizes";
	
	/** @brief File extension for GenBank format files */
	static const std::string GBK_EXTENSION = ".gbk";
	
	/** @brief File extension for FASTA nucleotide format files */
	static const std::string FNA_EXTENSION = ".fna";

	/**
	 * @enum FileType
	 * @brief Enumerates the supported genetic file formats.
	 * @details This enumeration identifies the different file formats that can be
	 * processed by the genetics utility functions. Each format has specific parsing
	 * rules implemented in the corresponding from_* functions.
	 */
	enum class FileType
	{
		GBK,            ///< GenBank format (.gbk)
		FNA,            ///< FASTA nucleotide format (.fna)
		DNA,            ///< Plain text DNA format (.dna)
		DNA_BIN,        ///< Binary DNA format (.dna.bin)
		DNA_BIN_COMBINED, ///< Combined binary DNA format (.dna.bin.combined)
		UNKNOWN         ///< Unknown or unsupported file format
	};

	/**
	 * @brief Loads a gene from a GenBank format file.
	 * @param filename Path to the GenBank (.gbk) file.
	 * @return Gene The parsed genetic sequence.
	 * @throws std::invalid_argument if the file cannot be opened.
	 * @see from_file
	 */
	Gene from_gbk(const std::string& filename);
	
	/**
	 * @brief Loads a gene from a FASTA nucleotide format file.
	 * @param filename Path to the FASTA nucleotide (.fna) file.
	 * @return Gene The parsed genetic sequence.
	 * @throws std::invalid_argument if the file cannot be opened.
	 * @see from_file
	 */
	Gene from_fna(const std::string& filename);
	
	/**
	 * @brief Loads a gene from a plain text DNA format file.
	 * @param filename Path to the DNA (.dna) file.
	 * @return Gene The parsed genetic sequence.
	 * @throws std::invalid_argument if the file cannot be opened.
	 * @see from_file
	 */
	Gene from_dna(const std::string& filename);
	
	/**
	 * @brief Loads a gene from a file with automatic format detection.
	 * @details Determines the file type based on its extension and calls
	 * the appropriate parsing function.
	 * @param filename Path to the genetic data file.
	 * @return Gene The parsed genetic sequence.
	 * @throws std::invalid_argument if the file cannot be opened or has an unsupported format.
	 * @see get_file_type
	 */
	Gene from_file(const std::string& filename);
	
	/**
	 * @brief Loads a gene from a file path with automatic format detection.
	 * @param path A filesystem path to the genetic data file.
	 * @return Gene The parsed genetic sequence.
	 * @throws std::invalid_argument if the file cannot be opened or has an unsupported format.
	 * @see from_file(const std::string&)
	 */
	Gene from_file(const std::filesystem::path& path);
	
	/**
	 * @brief Loads a gene from a directory entry with automatic format detection.
	 * @param entry A filesystem directory entry pointing to the genetic data file.
	 * @return Gene The parsed genetic sequence.
	 * @throws std::invalid_argument if the file cannot be opened or has an unsupported format.
	 * @see from_file(const std::filesystem::path&)
	 */
	Gene from_file(const std::filesystem::directory_entry& entry);
	
	/**
	 * @brief Saves a gene to a plain text DNA format file.
	 * @param gene The genetic sequence to save.
	 * @param filename Path where the DNA file will be saved.
	 * @throws std::runtime_error if the file cannot be created or written to.
	 */
	void to_dna(const Gene& gene, const std::string& filename);
	
	/**
	 * @brief Loads all genes from files in a directory with automatic format detection.
	 * @param dirname Path to the directory containing genetic data files.
	 * @return std::vector<Gene> Vector of parsed genetic sequences.
	 * @see from_directory(const std::string&, FileType)
	 */
	std::vector<Gene> from_directory(const std::string& dirname);
	
	/**
	 * @brief Loads all genes from files of a specific format in a directory.
	 * @param dirname Path to the directory containing genetic data files.
	 * @param type The specific file type to load.
	 * @return std::vector<Gene> Vector of parsed genetic sequences.
	 * @see from_file
	 */
	std::vector<Gene> from_directory(const std::string& dirname, FileType type);
	
	/**
	 * @brief Extracts DNA strings from all supported files in a directory.
	 * @param dirname Path to the directory containing genetic data files.
	 * @return std::vector<std::string> Vector of DNA strings.
	 */
	std::vector<std::string> get_dna_strings_from_directory(const std::string& dirname);
	
	/**
	 * @brief Determines the file type based on its filename extension.
	 * @param filename The filename to analyze.
	 * @return FileType The determined file type, or FileType::UNKNOWN if not recognized.
	 */
	FileType get_file_type(const std::string& filename);
	
	/**
	 * @brief Determines the file type based on its path's extension.
	 * @param path A filesystem path to analyze.
	 * @return FileType The determined file type, or FileType::UNKNOWN if not recognized.
	 * @see get_file_type(const std::string&)
	 */
	FileType get_file_type(const std::filesystem::path& path);
	
	/**
	 * @brief Determines the file type based on a directory entry's path extension.
	 * @param entry A filesystem directory entry to analyze.
	 * @return FileType The determined file type, or FileType::UNKNOWN if not recognized.
	 * @see get_file_type(const std::filesystem::path&)
	 */
	FileType get_file_type(const std::filesystem::directory_entry& entry);
	
	/**
	 * @brief Loads a gene from an input stream.
	 * @details Reads valid nucleotides from the stream until EOF, skipping invalid characters.
	 * @param stream The input stream to read from.
	 * @return Gene The parsed genetic sequence.
	 */
	Gene from_stream(std::istream& stream);
	
	/**
	 * @brief Gets paths to all virus nucleotide files in the virus directory.
	 * @return std::vector<std::filesystem::path> Vector of paths to virus nucleotide files.
	 */
	std::vector<std::filesystem::path> get_virus_nucleotide_paths();
	
	/**
	 * @brief Combines multiple binary DNA files into a single combined file.
	 * @details Creates a combined binary file and a separate size file that tracks
	 * the size of each original file for later extraction.
	 * @param dirname Directory containing binary DNA files to combine.
	 * @param combined_filename Name of the output combined file.
	 */
	void combine_dna_bin_files(const std::string& dirname, const std::string& combined_filename);

	/**
	 * @brief Manages efficient loading of multiple genetic sequences from a combined file.
	 * @details This class provides utilities for loading genes from a combined binary file
	 * created by combine_dna_bin_files(). It maintains an index of genes in the file for
	 * fast random access without loading the entire dataset into memory at once.
	 * 
	 * @see combine_dna_bin_files
	 */
	struct GeneManager
	{
	public:
		/**
		 * @brief Constructs a GeneManager for accessing genes in the specified path.
		 * @param path Path to the combined binary file containing genetic sequences.
		 * @throws std::runtime_error If the file cannot be opened or is improperly formatted.
		 */
		GeneManager(const std::string& path);
		
		/**
		 * @brief Destructor for GeneManager.
		 * @details Closes any open file streams.
		 */
		~GeneManager();

		/**
		 * @brief Retrieves all genes from the combined file.
		 * @return std::vector<Gene> A vector containing all genetic sequences from the file.
		 */
		std::vector<Gene> all_genes() const;
		
		/**
		 * @brief Retrieves a random subset of genes from the combined file.
		 * @param count The number of random genes to retrieve.
		 * @return std::vector<Gene> A vector containing the randomly selected genetic sequences.
		 */
		std::vector<Gene> get_random_genes(size_t count) const;
	private:
		/**
		 * @brief Structure to track information about a gene in the combined file.
		 * @details Stores metadata about each gene including its name, position in the file,
		 * and size, allowing for efficient direct access to individual genes.
		 */
		struct GeneInfo
		{
			std::string name;     ///< Name or identifier of the gene
			std::streampos position; ///< File position where the gene data begins
			size_t size;          ///< Size of the gene data in bytes
		};

		std::string m_path;       ///< Path to the combined binary file
		mutable std::ifstream m_file; ///< File stream for reading gene data (mutable to allow const methods to read)
		std::vector<GeneInfo> m_gene_infos; ///< Index of all genes in the combined file

		/**
		 * @brief Saves the size information of genes to a companion file.
		 * @details Creates a .sizes file alongside the combined binary file that
		 * stores metadata about each gene for efficient future loading.
		 */
		void save_sizes();
		
		/**
		 * @brief Loads gene size information from the companion file.
		 * @details Reads the .sizes file to populate the m_gene_infos vector
		 * with metadata about each gene in the combined binary file.
		 * @throws std::runtime_error If the sizes file cannot be opened or is corrupted.
		 */
		void load_sizes();
	};
};
