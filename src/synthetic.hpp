#pragma once

#include <string>
#include <vector>

std::string get_random_word(size_t length) noexcept;

std::vector<std::string> get_random_words(size_t length, size_t num_words, double mean_lcp_length) noexcept;

