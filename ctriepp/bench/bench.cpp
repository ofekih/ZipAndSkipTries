#include <unistd.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>

#ifdef __APPLE__
#include <mach/mach.h>
#endif

#include "city.h"

struct str_hash {
    std::size_t operator()(const char* key, std::size_t key_size) const {
        return CityHash64(key, key_size);
    }

    std::size_t operator()(const std::string& key) const {
        return CityHash64(key.c_str(), key.size());
    }
};

using value_type = uint8_t;

// #define CTriePP
#ifdef C_TRIE_PP

#include "../ctriepp/CTriePP.hpp"
using map_type = ctriepp::CTriePP<value_type>;

map_type construct() { return map_type(); }

void insert(map_type& map, const std::string& str, value_type value) {
    map.insert(&str, value);
}

bool prefixLocate(const map_type& map, const std::string& str) {
    return map.containsPrefix(str);
}

bool locate(const map_type& map, const std::string& str) {
    return map.contains(str);
}

void deleteKey(map_type& map, const std::string& str) { map.erase(str); }

void destroy(map_type& map) {}

#elif C_TRIE

#include "packed-ctrie/naive/naive_ctrie.cc"
using map_type = CompactTrie<value_type>;
map_type construct() { return map_type(); }
void insert(map_type& map, Edge* str, value_type value) { map.Update(str); }
bool prefixLocate(map_type& map, Edge* str) { return map.Member(str); }
bool locate(map_type& map, Edge* str) { return map.Member(str); }
void destroy(map_type& map) {}

#elif PACKED_C_TRIE_XOR

#include "packed-ctrie/xor/xor_packed_ctrie.cc"
using map_type = CompactTrie<value_type>;
map_type construct() { return map_type(); }
void insert(map_type& map, Edge* str, value_type value) { map.Update(str); }
bool prefixLocate(map_type& map, Edge* str) { return map.Member(str); }
bool locate(map_type& map, Edge* str) { return map.Member(str); }
void destroy(map_type& map) {}

#elif PACKED_C_TRIE_HASH

#include "packed-ctrie/hash/hash_packed_ctrie.cc"
using map_type = CompactTrie<value_type>;
map_type construct() { return map_type(); }
void insert(map_type& map, Edge* str, value_type value) { map.Update(str); }
bool prefixLocate(map_type& map, Edge* str) { return map.Member(str); }
bool locate(map_type& map, Edge* str) { return map.Member(str); }
void destroy(map_type& map) {}

#elif (BENCH_TSL_HTRIE_MAP || BENCH_TSL_HTRIE_MAP_LF_4 || \
       BENCH_TSL_HTRIE_MAP_LF_2 || BENCH_TSL_HTRIE_MAP_LF_1)

#include <tsl/htrie_map.h>

#if BENCH_KEY_SIZE_32
using map_type = tsl::htrie_map<char, value_type, str_hash, uint32_t>;
#else
using map_type = tsl::htrie_map<char, value_type, str_hash>;
#endif

map_type construct() {
    map_type map;

#ifdef BENCH_TSL_HTRIE_MAP_LF_4
    map.max_load_factor(4.0f);
#elif BENCH_TSL_HTRIE_MAP_LF_2
    map.max_load_factor(2.0f);
#elif BENCH_TSL_HTRIE_MAP_LF_1
    map.max_load_factor(1.0f);
#endif

    return map;
}

void insert(map_type& map, const std::string& str, value_type value) {
    map.insert(str, value);
}

bool prefixLocate(const map_type& map, const std::string& str) {
    auto iters = map.equal_prefix_range(str);
    // for (auto iter = iters.first; iter != iters.second; ++iter) {
    //     std::cout << iter.key() << std::endl;
    // }
    return iters.first != iters.second;
}

bool locate(const map_type& map, const std::string& str) {
    return map.find(str) != map.end();
}

void deleteKey(map_type& map, const std::string& str) { map.erase(str); }

void destroy(map_type& /*map*/) {}

#elif BENCH_CEDARPP

#include <cedar/cedarpp.h>
using map_type = cedar::da<value_type, -1, -2, true>;
static const size_t NUM_RESULT = 1;
map_type::result_triple_type result_triple[NUM_RESULT];

#define NO_CONSTRUCT_METHOD

void insert(map_type& map, const std::string& str, value_type value) {
    map.update(str.c_str(), str.size(), value);
}

bool prefixLocate(map_type& map, const std::string& str) {
    return map.commonPrefixPredict(str.c_str(), result_triple, NUM_RESULT) != 0;
}

bool locate(map_type& map, const std::string& str) {
    return map.exactMatchSearch<value_type>(str.c_str(), str.size()) !=
           value_type(cedar::NaN<int>::N1);
}

void deleteKey(map_type& map, const std::string& str) {
    map.erase(str.c_str(), str.size());
}

void destroy(map_type& /*map*/) {}

#elif BENCH_CEDARPP_UNORDERED

#include <cedar/cedarpp.h>
using map_type = cedar::da<value_type, -1, -2, false>;
static const size_t NUM_RESULT = 1;
map_type::result_triple_type result_triple[NUM_RESULT];

#define NO_CONSTRUCT_METHOD

void insert(map_type& map, const std::string& str, value_type value) {
    map.update(str.c_str(), str.size(), value);
}

bool prefixLocate(map_type& map, const std::string& str) {
    return map.commonPrefixPredict(str.c_str(), result_triple, NUM_RESULT) != 0;
}

bool locate(map_type& map, const std::string& str) {
    return map.exactMatchSearch<value_type>(str.c_str(), str.size()) != -1;
}

void destroy(map_type& /*map*/) {}

#elif Z_FAST_TRIE

#include "../ZFastTrie/ZFastTrie.hpp"

using map_type = ZFastTrie<char>;

map_type construct() { return map_type(); }

void insert(map_type& map, const std::string& str, value_type value) {
    map.insert(str, value);
}

bool prefixLocate(const map_type& map, const std::string& str) {
    return map.containsPrefix(str);
}

bool locate(const map_type& map, const std::string& str) {
    return map.contains(str);
}

void deleteKey(map_type& map, const std::string& str) { map.erase(str); }

void destroy(map_type& map) {}

#endif

#if (C_TRIE || PACKED_C_TRIE_HASH || PACKED_C_TRIE_XOR)
inline Edge* newEdge(const std::string& str) {
    Edge* edge = new Edge();
    edge->clear();
    for (std::size_t i = 0; i < str.size(); i++) {  // char to bit
        uint8_t c = str.c_str()[i];
        for (int j = 0; j < 8; j++) {
            ++textsize;

            if (c >= 128) {
                edge->push_back(true);
            } else {
                edge->push_back(false);
            }
            c = c << 1;
        }
    }
    return edge;
}
#endif

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

double insert_time_per_query = 0.;
double map_mem = 0.;

template <typename T>
void bench_insert(T& map, const std::vector<std::string>& lines) {
#if (C_TRIE || PACKED_C_TRIE_HASH || PACKED_C_TRIE_XOR)
    std::vector<Edge*> texts;
    for (auto& line : lines) {
        texts.push_back(newEdge(line));
    }
#endif

    const size_t mem_start = get_process_size();
    const auto chrono_start = std::chrono::high_resolution_clock::now();

#ifdef BENCH_RESERVE
    reserve(map, lines.size());
#endif

    int i = 0;
#if (C_TRIE || PACKED_C_TRIE_HASH || PACKED_C_TRIE_XOR)
    for (auto& line : texts) {
#else
    for (auto& line : lines) {
#endif
        insert(map, line, 1);
        i++;
    }

    const auto chrono_end = std::chrono::high_resolution_clock::now();
    const size_t mem_end = get_process_size();

    // std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(
    //                  chrono_end - chrono_start)
    //                      .count() /
    //                  double(lines.size())
    //           << " ns/key" << std::endl;
    // std::cout << double(mem_end - mem_start) / double(1024 * 1024) << " MiB
    // ("
    //           << (mem_end - mem_start) << " bytes)." << std::endl;
    insert_time_per_query =
        std::chrono::duration_cast<std::chrono::nanoseconds>(chrono_end -
                                                             chrono_start)
            .count() /
        double(lines.size());
    map_mem = double(mem_end - mem_start) / double(1024 * 1024);
}

double prefix_time_per_query = 0.;
template <typename T>
void bench_prefix(T& map, const std::vector<std::string>& lines) {
#if (C_TRIE || PACKED_C_TRIE_HASH || PACKED_C_TRIE_XOR)
    std::vector<Edge*> texts;
    for (auto& line : lines) {
        texts.push_back(newEdge(line));
    }
#endif

    const auto chrono_start = std::chrono::high_resolution_clock::now();

    std::size_t nb_found = 0;
#if (C_TRIE || PACKED_C_TRIE_HASH || PACKED_C_TRIE_XOR)
    for (auto& line : texts) {
#else
    for (auto& line : lines) {
#endif
        if (prefixLocate(map, line)) {
            nb_found++;
        } else {
            std::cerr << "Not Found query!" << std::endl;
            exit(1);
        }
    }

    const auto chrono_end = std::chrono::high_resolution_clock::now();

    // std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(
    //                  chrono_end - chrono_start)
    //                      .count() /
    //                  double(lines.size())
    //           << " ns/key" << std::endl;

    // std::cout << "Found " << nb_found << " of " << lines.size() << "
    prefix_time_per_query =
        std::chrono::duration_cast<std::chrono::nanoseconds>(chrono_end -
                                                             chrono_start)
            .count() /
        double(lines.size());
}

double locate_time_per_query = 0.;
template <typename T>
void bench_locate(T& map, const std::vector<std::string>& lines) {
#if (C_TRIE || PACKED_C_TRIE_HASH || PACKED_C_TRIE_XOR)
    std::vector<Edge*> texts;
    for (auto& line : lines) {
        texts.push_back(newEdge(line));
    }
#endif

    const auto chrono_start = std::chrono::high_resolution_clock::now();

    std::size_t nb_found = 0;
#if (C_TRIE || PACKED_C_TRIE_HASH || PACKED_C_TRIE_XOR)
    for (auto& line : texts) {
#else
    for (auto& line : lines) {
#endif
        if (locate(map, line)) {
            nb_found++;
        } else {
            std::cerr << "Not Found query!" << std::endl;
            std::cerr << line << std::endl;
            exit(1);
        }
    }

    const auto chrono_end = std::chrono::high_resolution_clock::now();

    // std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(
    //                  chrono_end - chrono_start)
    //                      .count() /
    //                  double(lines.size())
    //           << " ns/key" << std::endl;

    // std::cout << "Found " << nb_found << " of " << lines.size() << "
    locate_time_per_query =
        std::chrono::duration_cast<std::chrono::nanoseconds>(chrono_end -
                                                             chrono_start)
            .count() /
        double(lines.size());
    if (nb_found != lines.size()) {
        std::cerr << "Not Found query!" << std::endl;
    }
}

double delete_time_per_query = 0.;
template <typename T>
void bench_delete(T& map, const std::vector<std::string>& lines) {

    const auto chrono_start = std::chrono::high_resolution_clock::now();

    for (auto& line : lines) {
        deleteKey(map, line);
    }

    const auto chrono_end = std::chrono::high_resolution_clock::now();

    // std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(
    //                  chrono_end - chrono_start)
    //                      .count() /
    //                  double(lines.size())
    //           << " ns/key" << std::endl;

    // std::cout << "Found " << nb_found << " of " << lines.size() << "
    delete_time_per_query =
        std::chrono::duration_cast<std::chrono::nanoseconds>(chrono_end -
                                                             chrono_start)
            .count() /
        double(lines.size());

    std::size_t nb_found = 0;
    for (auto& line : lines) {
        if (locate(map, line)) {
            ++nb_found;
        }
    }
    if (nb_found != 0) {
        std::cerr << "Found deleted text! : " << nb_found << std::endl;
    }
}

void bench_insert_prefix(const std::string& file_words_insert,
                         const std::string& file_words_read) {
#ifdef NO_CONSTRUCT_METHOD
    map_type map;
#else
    auto map = construct();
#endif

    std::ifstream file(file_words_insert);
    if (!file.is_open()) {
        throw std::runtime_error("Couldn't read " + file_words_insert + ".");
    }

    std::vector<std::string> texts = std::vector<std::string>();

    std::string text;
    while (std::getline(file, text)) {
        if (text.empty()) {
            continue;
        }
        texts.push_back(text);
    }

    bench_insert(map, texts);

    std::ifstream qfile(file_words_read);
    if (!qfile.is_open()) {
        throw std::runtime_error("Couldn't read " + file_words_read + ".");
    }

    std::vector<std::string> querys;

    std::string query;
    double query_total_length = 0;
    while (std::getline(qfile, query)) {
        if (query.empty()) {
            continue;
        }
        querys.push_back(query);
        query_total_length += query.size();
    }
    double average_query_length = query_total_length / (double)querys.size();

    bench_prefix(map, querys);
    std::cout << average_query_length << " " << prefix_time_per_query
              << std::endl;
    destroy(map);
}

void bench_insert_locate(const std::string& file_words_insert,
                         const std::string& file_words_read) {
#ifdef NO_CONSTRUCT_METHOD
    map_type map;
#else
    auto map = construct();
#endif

    std::ifstream file(file_words_insert);
    if (!file.is_open()) {
        throw std::runtime_error("Couldn't read " + file_words_insert + ".");
    }

    std::vector<std::string> texts = std::vector<std::string>();

    std::string text;
    long texts_total_length = 0;
    while (std::getline(file, text)) {
        if (text.empty()) {
            continue;
        }
        texts.push_back(text);
        texts_total_length += text.size();
    }
    double texts_average_length =
        (double)texts_total_length / (double)texts.size();

    bench_insert(map, texts);

// 25#define TRACE_NODE
#ifdef TRACE_NODE
#if C_TRIE_PP
    TRACE(map_type::nodeFactory_.size());
    TRACE(AlphabetAwareZFastTrie<uint32_t>::nodeFactory_.size());
    TRACE(map_type::microTrieFactory_.size());
    long long macro_hash_size = 0;
    for (size_t i = 0; i < map_type::nodeFactory_.size(); i++) {
        auto& hash = map_type::nodeFactory_.at(i);
        if (hash.children_ != nullptr) {
            macro_hash_size += hash.children_->capacity();
        }
    }
    TRACE(macro_hash_size);
    long long micro_hash_size = 0;
    for (size_t i = 0;
         i < AlphabetAwareZFastTrie<uint32_t>::nodeFactory_.size(); i++) {
        auto& hash = AlphabetAwareZFastTrie<uint32_t>::nodeFactory_.at(i);
        micro_hash_size += hash.children_.capacity();
    }
    TRACE(micro_hash_size);
    long long h2n_hash_size = 0;
    long long h2n_hash_capacity = 0;
    for (size_t i = 0; i < map_type::microTrieFactory_.size(); i++) {
        auto& hash = map_type::microTrieFactory_.at(i);
        h2n_hash_size += hash.handle2NodeMap_.size();
        h2n_hash_capacity += hash.handle2NodeMap_.capacity();
        if (hash.handle2NodeMap_.capacity() /
                (hash.handle2NodeMap_.size() + 1) >
            2) {
            // INFO(hash.handle2NodeMap_.capacity()
            //      << " " << hash.handle2NodeMap_.size() << " "
            //      << hash.handle2NodeMap_.capacity() /
            //             (hash.handle2NodeMap_.size() + 1));
        }
    }
    TRACE(h2n_hash_size);
    TRACE(h2n_hash_capacity);
    TRACE(NUM_CUCKOO_HASH_FUNCTIONS);
#endif
#endif

    // #ifdef C_TRIE
    //     map.traversal();
    //     std::cout << result_nodenum << std::endl;
    // #endif
    // #ifdef C_TRIE_PP
    // long total = 0;
    //     std::map<long, long> size2count;
    //     // auto factory = *AlphabetAwareZFastTrie<Uint>::m_node_factory;
    //     // for (size_t i = 0; i < factory.size_; i++) {
    //     //     auto node = factory.at(i);
    //     //     long size = node.m_children.size();
    //     auto factory = map.nodeFactory_;
    //     for (size_t i = 0; i < factory.size_; i++) {
    //         auto node = map.node(i);
    //         if(node.microTrieIndex_ != -1){
    //             auto microTrie = map.microTrie(node.microTrieIndex_);
    //             long size = microTrie.m_handle2node_map.size();
    //             if(i == 0){
    //                 TRACE(size);
    //             }
    //             total += size;
    //             if(size2count.count(size) == 0){
    //                 size2count[size] = 1;
    //             }else{
    //                 ++size2count[size];
    //             }
    //         }
    //     }

    //     TRACE(total);
    //     TRACE(AlphabetAwareZFastTrie<Uint>::m_node_factory->size_);
    //     // for (auto i = count.begin(); i != count.end(); i++) {
    //     //         std::cout << i->first<<" "<<i->second << std::endl;
    //     // }

    //     std::vector<long> histgram(22, 0);
    //     for (auto i = size2count.begin(); i != size2count.end(); i++) {
    //         if(i->first == 0){
    //             ++histgram[0];
    //         }else{
    //             int pow = std::ceil(std::log2(i->first));
    //             histgram[pow + 1] += i->second;
    //         }

    //     }

    //     for (int i = 0; i < histgram.size(); i++) {
    //             std::cout << (i == 0 ? 0 : (1 << (i - 1))) << " " <<
    //             histgram[i] << std::endl;
    //     }

    // #endif

    std::ifstream qfile(file_words_read);
    std::vector<std::string> querys;
    std::string query;
    double query_total_length = 0;
    while (std::getline(qfile, query)) {
        if (query.empty()) {
            continue;
        }
        querys.push_back(query);
        query_total_length += query.size();
    }
    double average_query_length = query_total_length / (double)querys.size();

    bench_locate(map, querys);
// #if !(C_TRIE || PACKED_C_TRIE_HASH || PACKED_C_TRIE_XOR)
#if (C_TRIE_PP || BENCH_TSL_HTRIE_MAP || BENCH_TSL_HTRIE_MAP_LF_4 || \
     BENCH_TSL_HTRIE_MAP_LF_2 || BENCH_TSL_HTRIE_MAP_LF_1 ||         \
     Z_FAST_TRIE || BENCH_CEDARPP)
    bench_delete(map, texts);
#endif

    std::cout << texts_average_length << " " << map_mem << " "
              << insert_time_per_query << " " << locate_time_per_query << " "
              << delete_time_per_query << std::endl;
    destroy(map);
}

int main(int argc, char** argv) {
    if (argc == 3) {
        bench_insert_prefix(argv[1], argv[2]);
    } else if (argc == 4) {
        bench_insert_locate(argv[1], argv[2]);
    } else {
        std::cout << "Usage: " << argv[0] << " words_file_insert" << std::endl;
        std::cout << "Usage: " << argv[0]
                  << " words_file_insert words_file_read" << std::endl;
        exit(1);
    }
}
