#include <unistd.h>
#include <fstream>
#include <iostream>
#include <string>
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
    COUT(text_path);

    std::string rate[] = {"06", "07", "08", "09", "10"};

    for (size_t i = 0; i < 5; i++) {
        COUT(rate[i]);
        std::string out_text_path = argv[2];
        out_text_path += "." + rate[i];
        COUT(out_text_path);
        std::ifstream file(text_path);
        if (!file.is_open()) {
            throw std::runtime_error("Couldn't read " + text_path + ".");
        }

        std::ofstream out_file(out_text_path);
        std::string line;
        std::string out_line;
        while (std::getline(file, line)) {
            out_line =
                line.substr(0, line.size() * std::atoi(rate[i].c_str()) / 10);
            if (line.empty()) {
                continue;
            }
            out_file << out_line << std::endl;
        }
    }
}
