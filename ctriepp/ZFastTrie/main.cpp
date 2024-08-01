#include <unistd.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>

#include <city.h>
#include "ZFastTrie.hpp"

struct str_hash {
    std::size_t operator()(const char* key, std::size_t key_size) const {
        return CityHash64(key, key_size);
    }

    std::size_t operator()(const std::string& key) const {
        return CityHash64(key.c_str(), key.size());
    }
};

using value_type = int8_t;

using map_type = ZFastTrie<char>;

map_type construct() { return map_type(); }

void insert(map_type& map, const std::string& str, value_type value) {
    map.insert(str, value);
}

bool find(const map_type& map, const std::string& str) {
    return map.containsPrefix(str);
}

void destroy(map_type& map) {}

std::size_t get_memory_usage_bytes() {
    std::ifstream file("/proc/self/statm");
    std::size_t memory;
    file >> memory;  // Ignore first
    file >> memory;

    return memory * getpagesize();
}

size_t get_process_size() {
#ifdef __APPLE__
    struct task_basic_info t_info;
    mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
    task_info(current_task(), TASK_BASIC_INFO,
              reinterpret_cast<task_info_t>(&t_info), &t_info_count);
    return t_info.resident_size;
#else
    FILE* fp = std::fopen("/proc/self/statm", "r");
    size_t dummy(0), vm(0);
    std::fscanf(fp, "%ld %ld ", &dummy, &vm);  // get resident (see procfs)
    std::fclose(fp);
    return vm * ::getpagesize();
#endif
}

template <typename T>
void bench_insert(T& map, const std::string& file_words_insert) {
    std::ifstream file(file_words_insert);
    if (!file.is_open()) {
        throw std::runtime_error("Couldn't read " + file_words_insert + ".");
    }

    std::vector<std::string>* _lines = new std::vector<std::string>();
    std::vector<std::string>& lines = *_lines;

#ifdef SEP_SPACE
    // std::copy(std::istream_iterator<std::string>(file),
    //           std::istream_iterator<std::string>(),
    //           std::back_inserter(lines));
#else
    std::string line;
    while (file >> line) {
        if (line.empty()) {
            continue;
        }
        lines.push_back(line);
    }
#endif

    const size_t mem_start = get_process_size();
    const auto chrono_start = std::chrono::high_resolution_clock::now();

#ifdef BENCH_RESERVE
    reserve(map, lines.size());
#endif

    int i = 0;
    for (auto& line : lines) {
        insert(map, line, i);
        i++;
    }

    const auto chrono_end = std::chrono::high_resolution_clock::now();
    const size_t mem_end = get_process_size();

    std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(
                     chrono_end - chrono_start)
                         .count() /
                     double(lines.size())
              << " ns/key" << std::endl;
    std::cout << double(mem_end - mem_start) / double(1024 * 1024) << " MiB ("
              << (mem_end - mem_start) << " bytes)." << std::endl;
}

template <typename T>
void bench_read(T& map, const std::string& file_words_read) {
    std::ifstream file(file_words_read);
    if (!file.is_open()) {
        throw std::runtime_error("Couldn't read " + file_words_read + ".");
    }

    std::vector<std::string> lines;

#ifdef SEP_SPACE
    std::copy(std::istream_iterator<std::string>(file),
              std::istream_iterator<std::string>(), std::back_inserter(lines));
#else
    std::string line;
    while (file >> line) {
        if (line.empty()) {
            continue;
        }
        lines.push_back(line);
    }
#endif

    const auto chrono_start = std::chrono::high_resolution_clock::now();

    std::size_t nb_found = 0;
    for (auto& line : lines) {
        if (find(map, line)) {
            nb_found++;
        } else {
            exit(1);
        }
    }

    const auto chrono_end = std::chrono::high_resolution_clock::now();

    std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(
                     chrono_end - chrono_start)
                         .count() /
                     double(lines.size())
              << " ns/key" << std::endl;

    std::cout << "Found " << nb_found << " of " << lines.size() << " elements."
              << std::endl;
}

void bench(const std::string& file_words_insert,
           const std::string& file_words_read) {
#ifdef NO_CONSTRUCT_METHOD
    map_type map;
#else
    auto map = construct();
#endif

    bench_insert(map, file_words_insert);
    bench_read(map, file_words_read);

    destroy(map);
}

int main(int argc, char** argv) {
    if (argc == 2) {
        // insert_no_measure(argv[1]);
    } else if (argc == 3) {
        bench(argv[1], argv[2]);
    } else {
        std::cout << "Usage: " << argv[0] << " words_file_insert" << std::endl;
        std::cout << "Usage: " << argv[0]
                  << " words_file_insert words_file_read" << std::endl;
        exit(1);
    }
}
