#include <unistd.h>
#include <chrono>
#include <fstream>
#include <iostream>

#include "CTriePP.hpp"
#include "LongString.hpp"

int _main(int argc, char** argv) {
    // Edge *input_text = new Edge();

    if (argc < 2) {
        FATAL("");
    }

    std::cout << "text file : " << argv[1] << std::endl;

    std::ifstream ifs(argv[1]);
    if (ifs.fail()) {
        FATAL("file open failed");
        return -1;
    }

    std::vector<std::string> texts;
    // std::copy(std::istream_iterator<std::string>(ifs),
    //           std::istream_iterator<std::string>(),
    //           std::back_inserter(texts));
    std::string line;
    // while (ifs >> line) {
    while (std::getline(ifs, line)) {
        if (line.empty()) {
            continue;
        }
        texts.push_back(line);
    }

    std::cout << "texts size : " << texts.size() << std::endl;

    // start time
    std::chrono::system_clock::time_point start;
    // end time
    std::chrono::system_clock::time_point end;
    // duration = end - start;
    // std::chrono::system_clock::duration duration;
    Long msec;

    start = std::chrono::system_clock::now();
    ctriepp::CTriePP<bool> trie;
    for (size_t i = 0; i < texts.size(); ++i) {
        trie.insert(&(texts[i]), false);
    }

    end = std::chrono::system_clock::now();
    msec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
               .count();
    std::cout << "construction time " << msec << " [millisec]\n";

    std::vector<std::string>& patterns = texts;
    // std::ifstream ifs1(argv[1]);
    // if (ifs1.fail()) {
    //     FATAL("file open failed");
    //     return -1;
    // }

    // std::string pattern;
    // while (getline(ifs1, pattern)) {
    //     if (pattern.empty()) {
    //         continue;
    //     }
    //     patterns.push_back(new LongString(pattern));
    //     // patterns.push_back(
    //     //     new LongString(pattern.substr(0, pattern.size() - 1)));
    // }

    size_t NUM_BATCH = 1;
    Int NUM_QUERY = patterns.size();
    // mtries_msec = 0;
    start = std::chrono::system_clock::now();
    // assert(trie.contains(texts[19]));

    for (auto& pattern : patterns) {
        if (!trie.contains(&pattern)) {
            ERROR(pattern);
            exit(1);
            ERROR("NG");
        }
        if (!trie.containsPrefix(&pattern)) {
            ERROR(pattern);
            exit(1);
            ERROR("NG");
        }
    }

    for (auto& pattern : patterns) {
        trie.erase(pattern);
        if (trie.contains(&pattern)) {
            ERROR(pattern);
            exit(1);
            ERROR("NG");
        }
    }

    end = std::chrono::system_clock::now();
    msec = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
               .count();
    std::cout << "total query time " << msec << " [nanosec]\n";
    std::cout << "average query time " << msec / NUM_QUERY / NUM_BATCH
              << " [nanosec]\n";

    return 0;
}
