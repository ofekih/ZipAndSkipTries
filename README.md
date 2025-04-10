# Zip and Skip Tries

This repository contains the implementation and evaluation of several trie data structures with a focus on performance and memory efficiency, as described in the accompanying academic paper.

## Overview

The project implements and evaluates four main trie variants:

- **SkipTrie**: A sequential trie implementation that uses skipping techniques to improve traversal speed
- **ZipTrie**: A memory-efficient trie that compresses paths with single children
- **ParallelSkipTrie**: A GPU-accelerated version of SkipTrie using CUDA
- **ParallelZipTrie**: A GPU-accelerated version of ZipTrie using CUDA

Each implementation is designed to efficiently store and search for strings, with particular optimizations for genetic sequence data and other applications with long common prefixes.

## Documentation

You can find the full generated documentation [here](https://ofekih.github.io/ZipAndSkipTries/).

To generate the documentation locally:

```bash
make docs
```

The generated documentation will be available in the `docs/html` directory. Open `docs/html/index.html` in a web browser to view it.

## Requirements

- CUDA Toolkit (for parallel implementations)
- C++17 compatible compiler
- Python 3.x with matplotlib (for plotting)

## Building

The project uses a Makefile build system. To build the various components:

```bash
# Build the verification program
make verify

# Build the genetic data benchmark
make genetic_benchmark

# Build the synthetic data benchmark
make synthetic_benchmark
```

## Verification

To verify the correctness of the trie implementations:

```bash
./verify
```

This program tests all trie implementations with various inputs to ensure they produce correct results.

## Benchmarking

The project includes two benchmark programs:

### Genetic Benchmark

This benchmark evaluates performance on genetic data from the ABC-HuMi dataset:

```bash
./genetic_benchmark <num_trials> <num_simulations>
```

Example:
```bash
./genetic_benchmark 5 10
```

This will run 5 trials with 10 simulations each and save the results to the `data-genetic` directory.

### Synthetic Benchmark

This benchmark evaluates performance on synthetic data with controlled properties. The recommended way to run the synthetic benchmark is using the provided bash script:

```bash
bash run_synthetic_benchmark.sh
```

This script runs a series of benchmarks with appropriate parameters and ensures consistent testing across all trie implementations.

Alternatively, you can run the benchmark directly with custom parameters:

```bash
./synthetic_benchmark <num_words> <word_length> <mean_lcp> <num_repetitions>
```

Example:
```bash
./synthetic_benchmark 1000000 1024 512 5
```

This will generate 1,000,000 strings of length 1024 with a mean longest common prefix (LCP) length of 512, run 5 repetitions, and save the results to the `data-synthetic` directory.

## Plotting Results

The project includes Python scripts to visualize benchmark results:

### Genetic Data Plots

```bash
python plot_genetic.py
```

This generates plots comparing construction and search times for different trie implementations on genetic data.

### Synthetic Data Plots

```bash
python plot_synthetic.py
```

This generates plots comparing construction and search times for different trie implementations on synthetic data with controlled properties.

The plots are saved to the `figures` directory.

## Data Organization

The project organizes benchmark data into two directories:

- `data-genetic/`: Contains benchmark results from genetic data
- `data-synthetic/`: Contains benchmark results from synthetic data

## Key Features

- **Memory Efficiency**: ZipTrie variants compress paths to reduce memory usage
- **GPU Acceleration**: Parallel implementations leverage CUDA for high-performance operations
- **Genetic Data Support**: Specialized support for nucleotide sequences
- **Synthetic Data Generation**: Tools for creating controlled test data with specific properties
- **Comprehensive Benchmarking**: Detailed performance analysis across different parameters
