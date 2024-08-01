#include <unistd.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

typedef uint32_t Uint;

Uint hoge(Uint size) {
    static const float SHORT_RATE = 0.6f;
    static const float LONG_RATE = 1.0f;

    Uint short_size = size * SHORT_RATE;
    Uint long_size = size * LONG_RATE;

    return short_size + rand() % (long_size - short_size);
}

void make_query(const std::string& filePath, const std::string& out_filePath) {
    std::ifstream file(filePath);
    if (!file.is_open()) {
        throw std::runtime_error("Couldn't read " + filePath + ".");
    }

    std::ofstream out_file(out_filePath);
    if (!out_file.is_open()) {
        throw std::runtime_error("Couldn't read " + out_filePath + ".");
    }

    srand(0);

    // std::vector<std::string> texts;
    // std::copy(std::istream_iterator<std::string>(file),
    //           std::istream_iterator<std::string>(),
    //           std::back_inserter(texts));
    std::string text;
    while (file >> text) {
        if (text.empty()) {
            continue;
        }
        out_file << text.substr(0, hoge(text.size())) << std::endl;
    }
}

int main(int argc, char** argv) {
    if (argc == 3) {
        make_query(argv[1], argv[2]);
    } else {
        std::cout << "Usage: " << argv[0] << " in_file out_file" << std::endl;
        exit(1);
    }
}
