#include <fstream>
#include <iostream>
#include <set>
#include <vector>

#define LOG(message) std::cerr << "[DEBUG] " << message << std::endl
#define TRACE(message)                                     \
    std::cerr << "[TRACE] " << __FILE__ << "/" << __LINE__ \
              << "/ " #message " : " << message << std::endl
#define COUT(message) std::cout << message << std::endl

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

    std::string line;
    std::set<std::string> string_set;
    // while (file >> line) {
    while (std::getline(file, line)) {
        string_set.insert(line);
    }

    std::ofstream ofile(out_text_path);
    for (auto iter = string_set.begin(); iter != string_set.end(); ++iter) {
        // COUT(*iter);
        if (iter->empty()) {
            continue;
        }
        ofile << *iter << std::endl;
    }
}
