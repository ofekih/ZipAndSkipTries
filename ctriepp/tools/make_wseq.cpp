#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

#define LOG(message) std::cerr << "[DEBUG] " << message << std::endl
#define TRACE(message)                                     \
    std::cerr << "[TRACE] " << __FILE__ << "/" << __LINE__ \
              << "/ " #message " : " << message << std::endl
#define COUT(message) std::cout << message << std::endl

std::vector<std::string> split(const std::string& s, char delim) {
    std::vector<std::string> elems;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        if (!item.empty()) {
            elems.push_back(item);
        }
    }
    return elems;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        COUT("argument error");
        exit(1);
    }

    std::string text_path = argv[1];
    std::string out_text_path = argv[2];
    COUT(text_path);
    COUT(out_text_path);

    std::ifstream file(text_path);
    if (!file.is_open()) {
        throw std::runtime_error("Couldn't read " + text_path + ".");
    }

    std::string prev;
    std::string word;
    std::string seq;
    std::map<std::string, uint> string_set;
    while (file >> word) {
        if (prev.empty() && word.empty()) {
            continue;
        }
        seq = prev + " " + word;
        if (string_set.count(seq) == 0) {
            string_set[seq] = 1;
        } else {
            ++string_set[seq];
        }
        prev = word;
    }

    std::ofstream ofile(out_text_path);
    for (auto iter = string_set.begin(); iter != string_set.end(); ++iter) {
        // if (iter->second < 9) {
        //     continue;
        // }
        ofile << iter->first << std::endl;
    }
}
